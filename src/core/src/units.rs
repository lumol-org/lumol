// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

//! This module allow to convert from and to the internal unit system.
//!
//! Internal units are:
//!
//!  - Angstrom (A) for distances;
//!  - femtosecond (fs) for time;
//!  - Unified atomic mass unit (u or Da) for mass;
//!  - Kelvin (K) for temperature;
//!  - Number of particles for quantity of matter;
//!  - radian (rad) for angles;
//!
//! Other units are derived from these primitives units. For examples, the
//! internal unit for energy is 1e-4 kJ/mol.

use std::error::Error;
use std::fmt;
use std::num;

// Using a separated module because lazy_static does not support pub(crate)
mod detail {
    use consts::{BOHR_RADIUS, NA};
    use std::collections::BTreeMap;
    use std::f64::consts::PI;

    // Atomic mass unit in kg
    const U_IN_KG: f64 = 1.660538782e-27;

    /// Get the conversion factors from a string unit to the internal units.
    lazy_static!{
        pub static ref FACTORS: BTreeMap<&'static str, f64> = {
            let mut map = BTreeMap::new();
            // Distances units.
            assert!(map.insert("A", 1.0).is_none());
            assert!(map.insert("nm", 10.0).is_none());
            assert!(map.insert("pm", 1e-2).is_none());
            assert!(map.insert("fm", 1e-5).is_none());
            assert!(map.insert("m", 1e10).is_none());
            assert!(map.insert("bohr", BOHR_RADIUS).is_none());

            // Time units.
            assert!(map.insert("fs", 1.0).is_none());
            assert!(map.insert("ps", 1e3).is_none());
            assert!(map.insert("ns", 1e6).is_none());

            // Mass units.
            assert!(map.insert("u", 1.0).is_none());
            assert!(map.insert("Da", 1.0).is_none());
            assert!(map.insert("kDa", 1.0).is_none());
            assert!(map.insert("g", 1e-3 / U_IN_KG).is_none());
            assert!(map.insert("kg", 1.0 / U_IN_KG).is_none());

            // Temperature units.
            assert!(map.insert("K", 1.0).is_none());
            // Quantity of matter units.
            assert!(map.insert("mol", NA).is_none());

            // Angle units.
            assert!(map.insert("rad", 1.0).is_none());
            assert!(map.insert("deg", PI / 180.0).is_none());

            // Energy units.
            assert!(map.insert("J", 1e-10 / U_IN_KG).is_none());
            assert!(map.insert("kJ", 1e-7 / U_IN_KG).is_none());
            assert!(map.insert("kcal", 4.184 * 1e-7 / U_IN_KG).is_none());
            assert!(map.insert("eV", 1.60217653e-19 * 1e-10 / U_IN_KG).is_none());
            assert!(map.insert("H", 4.35974417e-18 * 1e-10 / U_IN_KG).is_none());
            assert!(map.insert("Ry", 4.35974417e-18 / 2.0 * 1e-10 / U_IN_KG).is_none());

            // Force unit.
            assert!(map.insert("N", 1e-20 / U_IN_KG).is_none());

            // Pressure units.
            assert!(map.insert("Pa", 1e-40 / U_IN_KG).is_none());
            assert!(map.insert("kPa", 1e-37 / U_IN_KG).is_none());
            assert!(map.insert("MPa", 1e-34 / U_IN_KG).is_none());
            assert!(map.insert("bar", 1e-35 / U_IN_KG).is_none());
            assert!(map.insert("atm", 101325.0 * 1e-40 / U_IN_KG).is_none());

            return map;
        };
    }
}

pub(crate) use self::detail::FACTORS;

/// Possible error causes when parsing an unit string.
#[derive(Debug)]
pub enum ParseError {
    /// Error while parsing a power in `x^y` expressions
    Power(num::ParseIntError),
    /// Error while parsing the value part of an unit string
    Value(num::ParseFloatError),
    /// Parentheses are not balanced in this unit
    ParenthesesMismatch,
    /// This unit was not found
    NotFound {
        /// The unit that created this error
        unit: String,
    },
    /// Any other error
    MalformedExpr(String),
}

impl From<num::ParseIntError> for ParseError {
    fn from(err: num::ParseIntError) -> ParseError {
        ParseError::Power(err)
    }
}

impl From<num::ParseFloatError> for ParseError {
    fn from(err: num::ParseFloatError) -> ParseError {
        ParseError::Value(err)
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            ParseError::Power(ref err) => err.fmt(f),
            ParseError::Value(ref err) => err.fmt(f),
            ParseError::ParenthesesMismatch => write!(f, "Parentheses are not equilibrated."),
            ParseError::NotFound { ref unit } => write!(f, "Unit '{}' not found.", unit),
            ParseError::MalformedExpr(ref err) => write!(f, "Malformed expression: {}", err),
        }
    }
}

impl Error for ParseError {
    fn description(&self) -> &str {
        match *self {
            ParseError::Power(ref err) => err.description(),
            ParseError::Value(ref err) => err.description(),
            ParseError::ParenthesesMismatch => "Parentheses are not equilibrated.",
            ParseError::NotFound { .. } => "Unit not found.",
            ParseError::MalformedExpr(..) => "Malformed expression",
        }
    }
}

/// Possible tokens in unit strings
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum Token {
    /// Left parentheses
    LParen,
    /// Right parentheses
    RParen,
    /// '*' token
    Mul,
    /// '/' token
    Div,
    /// '^' token
    Pow,
    /// Any other whitespaces separated value
    Value(String),
}

impl Token {
    /// What is the precedence of a specific token
    fn precedence(&self) -> usize {
        match *self {
            Token::LParen | Token::RParen => 0,
            Token::Div | Token::Mul => 10,
            Token::Pow => 20,
            Token::Value(..) => internal_error!("invalid call to UnitTok::precedence for values"),
        }
    }

    /// Get the string used to build this token in tokenize
    fn as_str(&self) -> &str {
        match *self {
            Token::LParen => "(",
            Token::RParen => ")",
            Token::Div => "/",
            Token::Mul => "*",
            Token::Pow => "^",
            Token::Value(ref value) => value,
        }
    }
}

/// Transform a string to a stream of tokens
fn tokenize(unit: &str) -> Vec<Token> {
    let mut tokens = Vec::new();
    let mut token = String::new();
    for c in unit.chars() {
        match c {
            '*' | '/' | '^' | '(' | ')' => {
                if !token.is_empty() {
                    tokens.push(Token::Value(token.clone()));
                    token.clear();
                }
                match c {
                    '*' => tokens.push(Token::Mul),
                    '/' => tokens.push(Token::Div),
                    '^' => tokens.push(Token::Pow),
                    '(' => tokens.push(Token::LParen),
                    ')' => tokens.push(Token::RParen),
                    _ => internal_error!("invalid unit operator"),
                }
            }
            other if !other.is_whitespace() => {
                token.push(other);
            }
            _ => assert!(c.is_whitespace()),
        }
    }
    // Last token
    if !token.is_empty() {
        tokens.push(Token::Value(token));
    }
    return tokens;
}

static MISSING_OPERATOR: &'static str = "Oops, sorry explorator, but you felt \
                                         in a space-time hole. We are missing an operator here";

/// Create the AST for unit expression using the Shunting-Yard algorithm.
///
/// See /// https://en.wikipedia.org/wiki/Shunting-yard_algorithm for a
/// description of the algorithm.
#[allow(trivial_casts)]
fn shunting_yard(tokens: Vec<Token>) -> Result<Vec<Token>, ParseError> {
    let mut operators = Vec::new();
    let mut output = Vec::new();
    for token in tokens {
        match token {
            Token::Value(..) => output.push(token),
            Token::Mul | Token::Div | Token::Pow => {
                while !operators.is_empty() {
                    // The cast is useless here, but rustc can't figure out
                    // the type of the expression after the call to `expect`
                    let top_operator =
                        (operators.last().expect(MISSING_OPERATOR) as &Token).clone();
                    // All the operators are left-associative
                    if token.precedence() <= top_operator.precedence() {
                        output.push(operators.pop().expect(MISSING_OPERATOR));
                    } else {
                        break;
                    }
                }
                operators.push(token);
            }
            Token::LParen => operators.push(token),
            Token::RParen => {
                while !operators.is_empty() && operators.last() != Some(&Token::LParen) {
                    output.push(operators.pop().expect(MISSING_OPERATOR))
                }
                if operators.is_empty() || operators.last() != Some(&Token::LParen) {
                    return Err(ParseError::ParenthesesMismatch);
                } else {
                    let _ = operators.pop();
                }
            }
        }
    }

    while !operators.is_empty() {
        match *operators.last().expect(MISSING_OPERATOR) {
            Token::LParen | Token::RParen => return Err(ParseError::ParenthesesMismatch),
            _ => output.push(operators.pop().expect(MISSING_OPERATOR)),
        }
    }

    return Ok(output);
}


/// Possible members in unit expressions
#[derive(Debug, PartialEq)]
enum UnitExpr {
    /// A single value
    Val(f64),
    /// Multiplication of left-hand side by right-hand side
    Mul(Box<UnitExpr>, Box<UnitExpr>),
    /// Division of left-hand side by right-hand side
    Div(Box<UnitExpr>, Box<UnitExpr>),
    /// Take the power of the expr by the `i32` value
    Pow(Box<UnitExpr>, i32),
}

impl UnitExpr {
    /// Recursively evaluate an unit expression
    fn eval(&self) -> f64 {
        match *self {
            UnitExpr::Val(v) => v,
            UnitExpr::Mul(ref lhs, ref rhs) => lhs.eval() * rhs.eval(),
            UnitExpr::Div(ref lhs, ref rhs) => lhs.eval() / rhs.eval(),
            UnitExpr::Pow(ref expr, pow) => expr.eval().powi(pow),
        }
    }

    /// Parse a string, and generate the corresponding unit expression
    fn parse(unit: &str) -> Result<UnitExpr, ParseError> {
        let tokens = tokenize(unit);
        let mut stream = try!(shunting_yard(tokens));
        let ast = try!(read_expr(&mut stream));
        if stream.is_empty() {
            Ok(ast)
        } else {
            let remaining = stream.iter().map(|t| t.as_str()).collect::<Vec<_>>().join(" ");
            return Err(ParseError::MalformedExpr(
                format!("remaining values after the end of the unit: {}", remaining),
            ));
        }
    }
}

/// Read and pop (recursively) a single expression from the `stream`.
/// The `stream` must be in reverse polish notation.
fn read_expr(stream: &mut Vec<Token>) -> Result<UnitExpr, ParseError> {
    if let Some(token) = stream.pop() {
        match token {
            Token::Value(unit) => {
                match FACTORS.get(&*unit) {
                    Some(&value) => Ok(UnitExpr::Val(value)),
                    None => Err(ParseError::NotFound { unit: unit }),
                }
            }
            Token::Mul => {
                let rhs = try!(read_expr(stream).map_err(|err| {
                    ParseError::MalformedExpr(format!("Error in unit at the right of '*': {}", err))
                }));
                let lhs = try!(read_expr(stream).map_err(|err| {
                    ParseError::MalformedExpr(format!("Error in unit at the left of '*': {}", err))
                }));
                Ok(UnitExpr::Mul(Box::new(lhs), Box::new(rhs)))
            }
            Token::Div => {
                let rhs = try!(read_expr(stream).map_err(|err| {
                    ParseError::MalformedExpr(format!("Error in unit at the right of '/': {}", err))
                }));
                let lhs = try!(read_expr(stream).map_err(|err| {
                    ParseError::MalformedExpr(format!("Error in unit at the left of '/': {}", err))
                }));
                Ok(UnitExpr::Div(Box::new(lhs), Box::new(rhs)))
            }
            Token::Pow => {
                let pow = match stream.pop() {
                    Some(pow) => {
                        match pow {
                            Token::Value(value) => try!(value.parse()),
                            _ => {
                                return Err(ParseError::MalformedExpr(
                                    format!("Invalid value after ^: {}", pow.as_str()),
                                ))
                            }
                        }
                    }
                    None => {
                        return Err(
                            ParseError::MalformedExpr(String::from("Missing value after '^'")),
                        )
                    }
                };
                let expr = try!(read_expr(stream).map_err(|err| {
                    ParseError::MalformedExpr(format!("Error in unit at the left of '*': {}", err))
                }));
                Ok(UnitExpr::Pow(Box::new(expr), pow))
            }
            Token::LParen | Token::RParen => {
                internal_error!("there should not be any parenthese here")
            }
        }
    } else {
        Err(ParseError::MalformedExpr(String::from("missing a value")))
    }
}

/// Convert the numeric value `val` from the unit `unit` to the internal unit.
///
/// ```
/// use lumol_core::units;
/// let internal = units::from(10.0, "A").unwrap();
/// assert!(internal == 10.0);
/// ```
pub fn from(value: f64, unit: &str) -> Result<f64, ParseError> {
    let unit = try!(UnitExpr::parse(unit));
    return Ok(unit.eval() * value);
}

/// Parse the string `val` and convert it to the corresponding internal unit
///
/// ```
/// use lumol_core::units;
/// let internal = units::from_str("10 A").unwrap();
/// assert!(internal == 10.0);
/// ```
pub fn from_str(value: &str) -> Result<f64, ParseError> {
    let unit = value.split_whitespace().skip(1).collect::<Vec<&str>>().join(" ");
    let unit = if unit.is_empty() {
        UnitExpr::Val(1.0)
    } else {
        try!(UnitExpr::parse(&unit))
    };
    let value = value.split_whitespace().take(1).collect::<Vec<&str>>()[0];
    let value = try!(value.parse::<f64>());
    return Ok(unit.eval() * value);
}

/// Convert the numeric value `val` (in internal units) to the unit `unit`.
///
/// ```
/// use lumol_core::units;
/// let real = units::to(10.0, "A").unwrap();
/// assert!(real == 10.0);
/// ```
pub fn to(value: f64, unit: &str) -> Result<f64, ParseError> {
    let unit = try!(UnitExpr::parse(unit));
    return Ok(value / unit.eval());
}

#[cfg(test)]
mod test {
    use super::*;
    use super::{Token, UnitExpr};
    use super::{shunting_yard, tokenize};

    #[test]
    fn tokens() {
        assert_eq!(tokenize("(")[0], Token::LParen);
        assert_eq!(tokenize(")")[0], Token::RParen);
        assert_eq!(tokenize("*")[0], Token::Mul);
        assert_eq!(tokenize("/")[0], Token::Div);
        assert_eq!(tokenize("^")[0], Token::Pow);
        assert_eq!(tokenize("foo")[0], Token::Value(String::from("foo")));
        assert_eq!(tokenize("45")[0], Token::Value(String::from("45")));

        assert_eq!(tokenize("(bar/m").len(), 4);
        assert_eq!(tokenize(" ( bar\t/\n   m").len(), 4);
    }

    fn ast_str(unit: &str) -> Result<String, ParseError> {
        let tokens = tokenize(unit);
        let ast = try!(shunting_yard(tokens));
        return Ok(ast.iter().map(|t| t.as_str()).collect::<Vec<_>>().join(" "));
    }

    #[test]
    fn ast() {
        assert_eq!(ast_str("").unwrap(), "");
        assert_eq!(ast_str("()").unwrap(), "");
        assert_eq!(ast_str("foo").unwrap(), "foo");
        assert_eq!(ast_str("foo*bar").unwrap(), "foo bar *");
        assert_eq!(ast_str("foo / bar").unwrap(), "foo bar /");
        assert_eq!(ast_str("foo^4").unwrap(), "foo 4 ^");
        assert_eq!(ast_str("bar/foo ^ 4").unwrap(), "bar foo 4 ^ /");
        assert_eq!(ast_str("k*bar /foo^ 4").unwrap(), "k bar * foo 4 ^ /");
    }

    #[test]
    fn ast_errors() {
        assert!(ast_str("(").is_err());
        assert!(ast_str(")").is_err());
        assert!(ast_str("(bar/m").is_err());
        assert!(ast_str("m/K)").is_err());
    }

    #[test]
    fn eval() {
        assert_eq!(UnitExpr::parse("A").unwrap(), UnitExpr::Val(1.0));
        assert_eq!(UnitExpr::parse("nm").unwrap(), UnitExpr::Val(10.0));

        assert_eq!(UnitExpr::parse("bohr/fs").unwrap().eval(), 0.52917720859);
        assert_eq!(UnitExpr::parse("(Ry / rad^-3   )").unwrap().eval(), 0.13127498789124938);
        assert_eq!(UnitExpr::parse("bar/(m * fs^2)").unwrap().eval(), 6.0221417942167636e-19);
        assert_eq!(UnitExpr::parse("kJ/mol/deg^2").unwrap().eval(), 0.3282806352310398);

        assert_ulps_eq!(UnitExpr::parse("kcal/mol/A^2").unwrap().eval(), 4.184e-4, epsilon = 1e-9);
    }

    #[test]
    fn parsing_errrors() {
        assert!(UnitExpr::parse("m^4-8").is_err());
        assert!(UnitExpr::parse("foo ^ bar").is_err());
        assert!(UnitExpr::parse("m^z4").is_err());
        assert!(UnitExpr::parse("HJK").is_err());
    }

    #[test]
    fn unit_from_str() {
        assert_eq!(from_str("10.0 A").unwrap(), 10.0);
        assert_eq!(from_str("10 A").unwrap(), 10.0);
        assert_eq!(from_str("1e1 A").unwrap(), 10.0);
        assert_eq!(from_str("10").unwrap(), 10.0);

        assert!(from_str("10a.0 bar").is_err());
        assert!(from_str("h10").is_err());
    }

    #[test]
    fn unit_to() {
        assert_eq!(to(25.0, "m").unwrap(), 2.5e-9);
        assert_eq!(to(25.0, "bar").unwrap(), 4.1513469550000005e9);
        assert_eq!(to(25.0, "kJ/mol").unwrap(), 249999.99982494753);
    }
}
