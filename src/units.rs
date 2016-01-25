/* Cymbalum, Molecular Simulation in Rust - Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
 */

//! This module allow to convert from and to the internal unit system, using the
//! following units:
//!
//!  - Angstrom (A) for distances;
//!  - femtosecond (fs) for time;
//!  - Unified atomic mass unit (u or Da) for mass;
//!  - Kelvin (K) for temperature;
//!  - Number of particles for quantity of matter;
//!  - radian (rad) for angles;
//!
//!  Using this system, the internal unit for energy is 1e-4 kJ/mol.

use std::collections::HashMap;
use std::f64::consts::PI;
use std::num;
use std::fmt;
use std::error::Error;

use constants::{BOHR_RADIUS, NA};

// Atomic mass unit in kg
const U_IN_KG : f64 = 1.660538782e-27;

/// Get the conversion factors from a string unit to the internal units.
fn conversion_factors() -> HashMap<&'static str, f64> {
    let mut map = HashMap::new();
    // Distances units.
    map.insert("A", 1.0);
    map.insert("Å", 1.0);
    map.insert("nm", 10.0);
    map.insert("pm", 1e-2);
    map.insert("fm", 1e-5);
    map.insert("m", 1e10);
    map.insert("bohr", BOHR_RADIUS);

    // Time units.
    map.insert("fs", 1.0);
    map.insert("ps", 1e3);
    map.insert("ns", 1e6);

    // Mass units.
    map.insert("u", 1.0);
    map.insert("Da", 1.0);
    map.insert("kDa", 1.0);
    map.insert("g", 1e-3 / U_IN_KG);
    map.insert("kg", 1.0 / U_IN_KG);

    // Temperature units.
    map.insert("K", 1.0);
    // Quantity of matter units.
    map.insert("mol", NA);

    // Angle units.
    map.insert("rad", 1.0);
    map.insert("deg", PI / 180.0);

    // Energy units.
    map.insert("J", 1e-10 / U_IN_KG);
    map.insert("kJ", 1e-7 / U_IN_KG);
    map.insert("kcal", 4.184 * 1e-7 / U_IN_KG);
    map.insert("eV", 1.60217653e-19 * 1e-10 / U_IN_KG);
    map.insert("H", 4.35974417e-18 * 1e-10 / U_IN_KG);
    map.insert("Ry", 4.35974417e-18 / 2.0 * 1e-10 / U_IN_KG);

    // Force unit.
    map.insert("N", 1e-20 / U_IN_KG);

    // Pressure units.
    map.insert("Pa", 1e-40 / U_IN_KG);
    map.insert("kPa", 1e-37 / U_IN_KG);
    map.insert("MPa", 1e-34 / U_IN_KG);
    map.insert("bar", 1e-35 / U_IN_KG);
    map.insert("atm", 101325.0 * 1e-40 / U_IN_KG);

    return map;
}

/// Possible error causes when parsing an unit string.
#[derive(Debug)]
pub enum UnitParsingError {
    /// Error while parsing a power in `x^y` expressions
    PowerParsingError(num::ParseIntError),
    /// Error while parsing the value part of an unit string
    ValueParsingError(num::ParseFloatError),
    /// Parentheses are not balanced in this unit
    UnbalancedParentheses{
        /// The full unit that created this error
        unit:String
    },
    /// Unknown binary operator
    BadBinary{
        /// The operator that created this error
        op: char
    },
    /// This unit was not found
    NotFound{
        /// The full unit that created this error
        unit:String
    },
}

impl From<num::ParseIntError> for UnitParsingError {
    fn from(err: num::ParseIntError) -> UnitParsingError {
        UnitParsingError::PowerParsingError(err)
    }
}

impl From<num::ParseFloatError> for UnitParsingError {
    fn from(err: num::ParseFloatError) -> UnitParsingError {
        UnitParsingError::ValueParsingError(err)
    }
}

impl fmt::Display for UnitParsingError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            UnitParsingError::PowerParsingError(ref err) => err.fmt(f),
            UnitParsingError::ValueParsingError(ref err) => err.fmt(f),
            UnitParsingError::UnbalancedParentheses{ref unit} => write!(f, "Parentheses are not equilibrated in unit: {}.", unit),
            UnitParsingError::BadBinary{ref op} => write!(f, "Bad binary operator {}.", op),
            UnitParsingError::NotFound{ref unit} => write!(f, "Unit '{}' not found.", unit),
        }
    }
}

impl Error for UnitParsingError {
    fn description(&self) -> &str {
        match *self {
            UnitParsingError::PowerParsingError(ref err) => err.description(),
            UnitParsingError::ValueParsingError(ref err) => err.description(),
            UnitParsingError::UnbalancedParentheses{..} => "Parentheses are not equilibrated.",
            UnitParsingError::BadBinary{..} => "Bad binary operator.",
            UnitParsingError::NotFound{..} => "Unit not found.",
        }
    }
}

/// Recursive unit string parsing function. This return the conversion factor
/// for a given unit.
fn conversion(unit: &str) -> Result<f64, UnitParsingError> {
    let unit = unit.trim();
    // First check if we do already have a known unit
    let factors = conversion_factors();
    if let Some(val) = factors.get(unit) {
        return Ok(*val);
    }

    let chars: Vec<char> = unit.chars().collect();
    let strlen = chars.len();

    if strlen == 0 {
        return Ok(1.0);
    }

    // First, check parentheses equilibration and any of '/' or '.' or '*' in the string
    let mut depth = 0;  // parentheses counter
    let mut i = strlen - 1;
    for (j, c) in unit.chars().rev().enumerate() {
        match c {
            '.' | '*' | '/' => {i = j;},
            '(' => {depth += 1;},
            ')' => {depth -= 1;},
            _ => {}
        }
    }
    let i = (strlen - 1) - i; // index of the first operator

    if depth != 0 {
        return Err(UnitParsingError::UnbalancedParentheses{unit: String::from(unit)});
    }

    // Are we enclosed by parentheses?
    if chars[0] == '(' {
        assert!(chars[strlen-1] == ')');
        return conversion(&unit[1..strlen-1])
    }

    // If a product character was found, recurse
    if i > 0 {
        let lhs = try!(conversion(&unit[..i]));
        let rhs = try!(conversion(&unit[i+1..]));

        return match chars[i] {
            '/' => Ok(lhs / rhs),
            '*' | '.' => Ok(lhs * rhs),
            _ => Err(UnitParsingError::BadBinary{op: chars[i]}),
        };
    }

    // Do we have an exponentiation?
    let mut i = strlen - 1;
    while chars[i].is_numeric() || chars[i] == '-' {
        i -= 1;
    }
    if chars[i] == '^' {
        let power = try!(unit[i+1..].parse::<i32>());
        let val = try!(conversion(&unit[0..i]));
        return Ok(val.powi(power));
    }

    return Err(UnitParsingError::NotFound{unit: String::from(unit)});
}

/// Convert the numeric value `val` from the unit `unit` to the internal unit.
///
/// ```
/// let internal = from(10.0, "A");
/// assert!(internal == 10.0);
///
/// let internal = from(1.0, "kJ/mol");
/// assert!(internal == 1000.0);
/// ```
pub fn from(val: f64, unit: &str) -> Result<f64, UnitParsingError> {
    let factor = try!(conversion(unit));
    return Ok(factor * val);
}

/// Parse the string `val` and convert it to the corresponding internal unit
///
/// ```
/// let internal = from_str("10 A");
/// assert!(internal == 10.0);
///
/// let internal = from("1 kJ/mol");
/// assert!(internal == 1000.0);
/// ```
pub fn from_str(val: &str) -> Result<f64, UnitParsingError> {
    let unit = val.split_whitespace().skip(1).collect::<Vec<&str>>().join(" ");
    let factor = try!(conversion(&unit));
    let val: &str = val.split_whitespace().take(1).collect::<Vec<&str>>()[0];
    let val: f64 = try!(val.parse());
    return Ok(factor * val);
}

/// Convert the numeric value `val` (in internal units) to the unit `unit`.
///
/// ```
/// let real = to(10.0, "A");
/// assert!(real == 10.0);
///
/// let real = to(1e-4, "kJ/mol");
/// assert!(internal == 1.0);
/// ```
pub fn to(val: f64, unit: &str) -> Result<f64, UnitParsingError> {
    let factor = try!(conversion(unit));
    return Ok(val / factor);
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn conversions() {
        // Get the factors rights
        assert_eq!(from(10.0, "A").ok(), Some(10.0));
        assert_eq!(from(10.0, "Å").ok(), Some(10.0));
        assert_eq!(from(10.0, "pm").ok(), Some(0.1));
        assert_eq!(from(10.0, "nm").ok(), Some(100.0));

        assert_eq!(from_str("10.0 A").ok(), Some(10.0));
        assert_eq!(from_str("10 A").ok(), Some(10.0));
        assert_eq!(from_str("1e1 A").ok(), Some(10.0));
        assert_eq!(from_str("10").ok(), Some(10.0));

        // Parse stuff
        assert_eq!(from(10.0, "bohr/fs").ok(), Some(10.0*0.52917720859));
        assert_eq!(from(10.0, "kJ/mol/A^2").ok(), Some(0.0010000000007002099));
        assert_eq!(from(10.0, "(Ry / rad^-3   )").ok(), Some(1.3127498789124938));
        assert_eq!(from(10.0, "bar/(m * fs^2)").ok(), Some(6.0221417942167636e-18));

        assert_eq!(from_str("10 bar/(m * fs^2)").ok(), Some(6.0221417942167636e-18));

        // Test the 'to' function too
        assert_eq!(to(25.0, "m").ok(), Some(2.5e-9));
        assert_eq!(to(25.0, "bar").ok(), Some(4.1513469550000005e9));
        assert_eq!(to(25.0, "kJ/mol").ok(), Some(249999.99982494753));
    }

    #[test]
    fn errors() {
        assert!(from(10.0, "(bar/m").is_err());
        assert!(from(10.0, "m^4-8").is_err());
        assert!(from(10.0, "m^z4").is_err());
        assert!(from(10.0, "m/K)").is_err());
        assert!(from(10.0, "HJK").is_err());

        assert!(from_str("10a.0 bar").is_err());
        assert!(from_str("h10").is_err());
    }
}
