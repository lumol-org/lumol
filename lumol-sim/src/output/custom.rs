// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

use std::error;
use std::fmt;
use std::fs::File;
use std::io::{self, BufWriter};
use std::io::prelude::*;
use std::path::{Path, PathBuf};

use caldyn::{Context, Expr};
use caldyn::Error as CaldynError;

use log::error;
use log_once::{warn_once, error_once};

use super::Output;
use lumol_core::{units, System};

/// Possible causes of error when using a custom output
#[derive(Debug)]
pub enum CustomOutputError {
    /// Any IO error
    Io(io::Error),
    /// Error in the mathematical expression
    Expr(CaldynError),
    /// Any other error
    Custom(String),
}

impl From<io::Error> for CustomOutputError {
    fn from(error: io::Error) -> CustomOutputError {
        CustomOutputError::Io(error)
    }
}

impl From<CaldynError> for CustomOutputError {
    fn from(error: CaldynError) -> CustomOutputError {
        CustomOutputError::Expr(error)
    }
}

impl From<String> for CustomOutputError {
    fn from(error: String) -> CustomOutputError {
        CustomOutputError::Custom(error)
    }
}

impl fmt::Display for CustomOutputError {
    fn fmt(&self, fmt: &mut fmt::Formatter<'_>) -> Result<(), fmt::Error> {
        match *self {
            CustomOutputError::Io(ref err) => write!(fmt, "{err}")?,
            CustomOutputError::Expr(ref err) => write!(fmt, "{err}")?,
            CustomOutputError::Custom(ref err) => write!(fmt, "{err}")?,
        }
        Ok(())
    }
}

impl error::Error for CustomOutputError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match *self {
            CustomOutputError::Io(ref err) => Some(err),
            CustomOutputError::Expr(ref err) => Some(err),
            CustomOutputError::Custom(_) => None,
        }
    }
}

/// Helper struct to parse and format custom output strings
struct FormatArgs {
    /// Pairs of "constant string", "format expression"
    args: Vec<(String, Expr)>,
    /// Any remaining tail after the last expression
    tail: String,
}

impl FormatArgs {
    fn new(format: &str) -> Result<FormatArgs, CustomOutputError> {
        let mut args = Vec::new();
        let mut expr = String::new();
        let mut tail = String::new();

        let mut in_expr = false;
        for c in format.chars() {
            match c {
                '{' if !in_expr => {
                    in_expr = true;
                }
                '}' if in_expr => {
                    in_expr = false;
                    let sub_expr = Expr::parse(&expr)?;
                    args.push((tail.clone(), sub_expr));
                    tail.clear();
                    expr.clear();
                }
                '{' if in_expr => {
                    return Err(CustomOutputError::Custom("found { in an expression".into()));
                }
                '}' if !in_expr => {
                    return Err(
                        CustomOutputError::Custom("found } outside of an expression".into()),
                    );
                }
                c => {
                    if in_expr {
                        expr.push(c);
                    } else {
                        tail.push(c);
                    }
                }
            }
        }
        if in_expr {
            return Err(CustomOutputError::Custom("mismatched braces".into()));
        }

        Ok(FormatArgs {
            args: args,
            tail: tail,
        })
    }

    fn format(&self, system: &System) -> Result<String, CustomOutputError> {
        let context = get_output_context(system);
        let mut output = String::new();
        for (string, expr) in &self.args {
            output.push_str(string);
            let value = expr.eval(&context)?;
            output.push_str(&value.to_string());
        }
        output.push_str(&self.tail);
        return Ok(output);
    }
}

fn get_output_context(system: &System) -> Context<'_> {
    let mut context = Context::new();
    context.set_query(move |name| {
        // Get unit conversion factor firsts
        units::CONVERSION_FACTORS.get(name).copied().or_else(|| {
            macro_rules! get_particle_data {
                ($index: ident, $data: ident) => (
                    system.particles()
                          .$data
                          .get($index)
                          .cloned()
                          .unwrap_or_else(|| {
                              warn_once!(
                                  "index out of bound in custom output: \
                                  index is {}, but we only have {} atoms",
                                  $index, system.size()
                              );
                              return num_traits::Zero::zero();
                          })
                );
            }
            if name.contains('[') {
                // vector data
                let (name, index) = parse_index(name);
                match name {
                    // position
                    "x" => Some(get_particle_data!(index, position)[0]),
                    "y" => Some(get_particle_data!(index, position)[1]),
                    "z" => Some(get_particle_data!(index, position)[2]),
                    // velocity
                    "vx" => Some(get_particle_data!(index, velocity)[0]),
                    "vy" => Some(get_particle_data!(index, velocity)[1]),
                    "vz" => Some(get_particle_data!(index, velocity)[2]),
                    // other atomic properties
                    "mass" => Some(get_particle_data!(index, mass)),
                    "charge" => Some(get_particle_data!(index, charge)),
                    _ => None,
                }
            } else {
                // scalar data
                match name {
                    "step" => Some(system.step as f64),
                    "pressure" => Some(system.pressure()),
                    "volume" => Some(system.volume()),
                    "temperature" => Some(system.temperature()),
                    "natoms" => Some(system.size() as f64),
                    "cell.a" => Some(system.cell.a()),
                    "cell.b" => Some(system.cell.b()),
                    "cell.c" => Some(system.cell.c()),
                    "cell.alpha" => Some(system.cell.alpha()),
                    "cell.beta" => Some(system.cell.beta()),
                    "cell.gamma" => Some(system.cell.gamma()),
                    "stress.xx" => Some(system.stress()[0][0]),
                    "stress.yy" => Some(system.stress()[1][1]),
                    "stress.zz" => Some(system.stress()[2][2]),
                    "stress.xy" => Some(system.stress()[0][1]),
                    "stress.xz" => Some(system.stress()[0][2]),
                    "stress.yz" => Some(system.stress()[1][2]),
                    _ => None,
                }
            }
        })
    });

    return context;
}

/// Get the name and index in a string looking like `name[index]`. Everything
/// else is just passed through.
fn parse_index(input: &str) -> (&str, usize) {
    // We can index `input`, because caldyn only works with ASCII data
    let l_brackets = input.match_indices('[').collect::<Vec<_>>();
    let r_brackets = input.match_indices(']').collect::<Vec<_>>();

    if l_brackets.len() != 1 || r_brackets.len() != 1 {
        // More than one bracket
        return (input, 0);
    }

    let start = l_brackets[0].0;
    let end = r_brackets[0].0;
    if start > end {
        // `[` is after `]`
        return (input, 0);
    }

    if let Ok(index) = input[(start + 1)..end].parse() {
        return (&input[..start], index);
    } else {
        // invalid integer value
        return (input, 0);
    }
}

/// The `CustomOutput` writes data into a file from an user-provided template.
///
/// The template string can contain mathematical expressions, using some
/// physical properties of the system. These mathematical expressions must be
/// enclosed in braces (`{}`). Here are some examples:
///
/// - A constant string is reproduced as it is: `some data`;
/// - Anything in braces is replaced by the corresponding values: `{pressure}
///   {volume}`;
/// - Mathematical operators are allowed in braces: `{pressure / volume}`. You
///   can use `+`, `-`, `/`, `*`, `^` for exponentiation and parentheses;
/// - Some properties are arrays of atomic properties `{x[0] + y[20]}`;
/// - Finally, all the properties are given in the internal units. One can
///   specify another unit: `x[0] / nm`.
///
/// Here is a list of all accepted properties:
///
/// - Atomic properties: `x`, `y` and `z` for cartesian coordinates, `vx`, `vy`
///   and `vz` for cartesian components of the velocity , `mass` for the atomic
///   mass, `charge` for the atomic charge.
/// - Physical properties: `pressure`, `volume`, `temperature`, `natoms`, stress
///   tensor components: `stress.xx`, `stress.yy`, `stress.zz`, `stress.xy`,
///   `stress.xz`, `stress.yz`, simulation `step`.
/// - Unit Cell properties: `cell.a`, `cell.b`, `cell.c` are the unit cell
///   vector lengths; `cell.alpha`, `cell.beta` and `cell.gamma` are the unit
///   cell angles.
pub struct CustomOutput {
    file: BufWriter<File>,
    path: PathBuf,
    template: String,
    args: FormatArgs,
}

impl CustomOutput {
    /// Create a new `CustomOutput` writing to the file at `filename` using
    /// the given `template`. The `template` is only partially validated at
    /// this stage.
    pub fn new<P: AsRef<Path>>(
        filename: P,
        template: &str,
    ) -> Result<CustomOutput, CustomOutputError> {
        Ok(CustomOutput {
            file: BufWriter::new(File::create(filename.as_ref())?),
            path: filename.as_ref().to_owned(),
            template: template.into(),
            args: FormatArgs::new(template)?,
        })
    }
}

impl Output for CustomOutput {
    fn setup(&mut self, _: &System) {
        writeln_or_log!(self, "# Custom output");
        writeln_or_log!(self, "# {}", self.template);
    }

    fn write(&mut self, system: &System) {
        if let Ok(formatted) = self.args.format(system) {
            writeln_or_log!(self, "{}", formatted);
        } else {
            error_once!("Could not evaluate custom output {}", self.template);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::tests::{test_output, testing_system};

    fn format(input: &str) -> String {
        FormatArgs::new(input).unwrap().format(&testing_system()).unwrap()
    }

    #[test]
    fn parsing_index() {
        assert_eq!(parse_index("a[6]"), ("a", 6));

        assert_eq!(parse_index("a"), ("a", 0));
        assert_eq!(parse_index("a][6"), ("a][6", 0));
        assert_eq!(parse_index("a[6][2]"), ("a[6][2]", 0));
        assert_eq!(parse_index("a[6]2]"), ("a[6]2]", 0));
        assert_eq!(parse_index("a[6][2"), ("a[6][2", 0));
        assert_eq!(parse_index("a[b]"), ("a[b]", 0));
    }

    #[test]
    fn format_args_parsing() {
        assert!(FormatArgs::new("one {test} two {5 } three!").is_ok());

        assert!(FormatArgs::new("{3 + 4} {").is_err());
        assert!(FormatArgs::new("{3 + 4} }").is_err());
        assert!(FormatArgs::new("{3 + { 4}").is_err());
        assert!(FormatArgs::new("{3 + {} }").is_err());
    }

    #[test]
    fn formating() {
        assert_eq!(format("{3 + 4}"), "7");

        assert_eq!(format("{pressure / bar}"), "10299.991728079816");
        assert_eq!(format("{temperature / K}"), "38083.04389172312");
        assert_eq!(format("{volume / A^3}"), "1000");

        assert_eq!(format("{cell.a / A}"), "10");
        assert_eq!(format("{cell.b / A}"), "10");
        assert_eq!(format("{cell.c / A}"), "10");
        assert_eq!(format("{cell.alpha}"), "90");
        assert_eq!(format("{cell.beta}"), "90");
        assert_eq!(format("{cell.gamma}"), "90");

        assert_eq!(format("{stress.xx / bar}"), "30899.975184239443");
        assert_eq!(format("{stress.yy / bar}"), "0");
        assert_eq!(format("{stress.zz / bar}"), "0");
        assert_eq!(format("{stress.xy / bar}"), "0");
        assert_eq!(format("{stress.xz / bar}"), "0");
        assert_eq!(format("{stress.yz / bar}"), "0");

        assert_eq!(format("{x[1]}"), "1.3");
        assert_eq!(format("{vy[1]}"), "0");
        assert_eq!(format("{vx[0]}"), "0.1");

        assert_eq!(format("{cell.a / bohr}"), "18.897261328856434");
        assert_eq!(format("{cell.a / nm}"), "1");
        assert_eq!(format("{cell.a / m}"), "0.000000001");

        assert_eq!(format("{step}"), "42");
    }

    #[test]
    fn custom() {
        let template = "p {pressure/bar} t {3 * 5} \tff";
        test_output(
            |path| Box::new(CustomOutput::new(path, template).unwrap()),
            "# Custom output
            # p {pressure/bar} t {3 * 5} \tff
            p 10299.991728079816 t 15 \tff
            ",
        );
    }
}
