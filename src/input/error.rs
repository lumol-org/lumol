// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

use std::io;
use std::error;
use std::fmt;
use std::result;

use toml::Parser;
use chemfiles;

use units::ParseError;

/// Custom `Result` type for input files
pub type Result<T> = result::Result<T, Error>;

/// Possible causes of error when reading input files
#[derive(Debug)]
pub enum Error {
    /// Error in the TOML input file
    TOML(String),
    /// IO error
    File(io::Error),
    /// Trajectory creation error
    Trajectory(TrajectoryError),
    /// File content error: missing sections, bad data types
    Config(String),
    /// Unit parsing error
    Unit(ParseError),
}

impl From<io::Error> for Error {
    fn from(err: io::Error) -> Error {Error::File(err)}
}

impl From<TrajectoryError> for Error {
    fn from(err: TrajectoryError) -> Error {Error::Trajectory(err)}
}

impl From<chemfiles::Error> for Error {
    fn from(err: chemfiles::Error) -> Error {
        Error::Trajectory(TrajectoryError::from(err))
    }
}

impl<'a> From<&'a str> for Error {
    fn from(err: &'a str) -> Error {
        Error::Config(String::from(err))
    }
}

impl From<String> for Error {
    fn from(err: String) -> Error {
        Error::Config(err)
    }
}

impl From<ParseError> for Error {
    fn from(err: ParseError) -> Error {Error::Unit(err)}
}

impl fmt::Display for Error {
    fn fmt(&self, fmt: &mut fmt::Formatter) -> result::Result<(), fmt::Error> {
        use std::error::Error as StdError;
        try!(write!(fmt, "{}", self.description()));
        Ok(())
    }
}

impl error::Error for Error {
    fn description(&self) -> &str {
        match *self {
            Error::TOML(ref err) | Error::Config(ref err) => err,
            Error::File(ref err) => err.description(),
            Error::Trajectory(ref err) => err.description(),
            Error::Unit(ref err) => err.description(),
        }
    }

    fn cause(&self) -> Option<&error::Error> {
        match *self {
            Error::TOML(..) | Error::Config(..) => None,
            Error::File(ref err) => Some(err),
            Error::Trajectory(ref err) => Some(err),
            Error::Unit(ref err) => Some(err),
        }
    }
}

pub fn toml_error_to_string(parser: &Parser) -> String {
    let errors = parser.errors.iter().map(|error|{
        let (line, _) = parser.to_linecol(error.lo);
        format!("{} at line {}", error.desc, line + 1)
    }).collect::<Vec<_>>().join("\n    ");

    let plural = if errors.len() == 1 {""} else {"s"};
    return format!("TOML parsing error{}: {}", plural, errors);
}

/******************************************************************************/

/// Possible error when reading and wrinting to trajectories
#[derive(Debug)]
pub struct TrajectoryError(chemfiles::Error);

impl From<chemfiles::Error> for TrajectoryError {
    fn from(err: chemfiles::Error) -> TrajectoryError {
        TrajectoryError(err)
    }
}

impl fmt::Display for TrajectoryError {
    fn fmt(&self, fmt: &mut fmt::Formatter) -> result::Result<(), fmt::Error> {
        self.0.fmt(fmt)
    }
}

impl error::Error for TrajectoryError {
    fn description(&self) -> &str {
        self.0.description()
    }
}
