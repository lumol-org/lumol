// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license
use std::fmt;
use std::io;
use std::path::PathBuf;

use lumol_sim::output::CustomOutputError;
use lumol_core::TrajectoryError;
use lumol_core::units::ParseError;

/// Possible causes of error when reading input files
#[derive(Debug)]
pub enum Error {
    /// Error in the TOML input file
    TOML(Box<dyn std::error::Error>),
    /// IO error, and associated file path
    Io(io::Error, PathBuf),
    /// Error while reading a trajectory file
    Trajectory(TrajectoryError),
    /// File content error: missing sections, bad data types
    Config(String),
    /// Unit parsing error
    Unit(ParseError),
    /// Specific error from the custom outputs
    CustomOutput(CustomOutputError),
}

impl From<(io::Error, PathBuf)> for Error {
    fn from((err, path): (io::Error, PathBuf)) -> Error {
        Error::Io(err, path)
    }
}

impl From<TrajectoryError> for Error {
    fn from(err: TrajectoryError) -> Error {
        Error::Trajectory(err)
    }
}

// impl From<chemfiles::Error> for Error {
//     fn from(err: chemfiles::Error) -> Error {
//         Error::Trajectory(TrajectoryError::from(err))
//     }
// }

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
    fn from(err: ParseError) -> Error {
        Error::Unit(err)
    }
}

impl From<(CustomOutputError, PathBuf)> for Error {
    fn from((err, path): (CustomOutputError, PathBuf)) -> Error {
        match err {
            CustomOutputError::Io(err) => Error::Io(err, path),
            other => Error::CustomOutput(other),
        }
    }
}

impl fmt::Display for Error {
    fn fmt(&self, fmt: &mut fmt::Formatter<'_>) -> Result<(), fmt::Error> {
        match *self {
            Error::Io(ref err, ref path) => {
                match err.kind() {
                    io::ErrorKind::NotFound => {
                        write!(fmt, "can not find '{}'", path.display())
                    }
                    io::ErrorKind::PermissionDenied => {
                        write!(fmt, "permission to access '{}' denied", path.display())
                    }
                    _ => {
                        write!(fmt, "error with '{}': {}", path.display(), err)
                    }
                }
            }
            Error::Trajectory(ref err) => write!(fmt, "{err}"),
            Error::TOML(ref err) => write!(fmt, "{err}"),
            Error::Config(ref err) => write!(fmt, "{err}"),
            Error::Unit(ref err) => write!(fmt, "{err}"),
            Error::CustomOutput(ref err) => write!(fmt, "{err}"),
        }
    }
}

impl std::error::Error for Error {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match *self {
            Error::TOML(..) | Error::Config(..) => None,
            Error::Io(ref err, _) => Some(err),
            Error::Trajectory(ref err) => Some(err),
            Error::Unit(ref err) => Some(err),
            Error::CustomOutput(ref err) => Some(err),
        }
    }
}
