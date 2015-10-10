/* Cymbalum, Molecular Simulation in Rust - Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
 */
//! Reading input files using the YAML format.
//!
//! These files help to define an universe or a simulation, in an human-readable
//! way.

extern crate yaml_rust as yaml;
use self::yaml::ScanError;

use std::io;
use std::result;

use units::UnitParsingError;

#[derive(Debug)]
/// Possible causes of error when reading potential files
pub enum Error{
    /// Error in the YAML input file
    YAMLError(ScanError),
    /// IO error
    FileError(io::Error),
    /// File content error: missing sections, bad data types
    ConfigError{
        /// Error message
        msg: String,
    },
    /// Unit parsing error
    UnitError(UnitParsingError),
}

impl From<ScanError> for Error {
    fn from(err: ScanError) -> Error {Error::YAMLError(err)}
}

impl From<io::Error> for Error {
    fn from(err: io::Error) -> Error {Error::FileError(err)}
}

impl<'a> From<&'a str> for Error {
    fn from(err: &'a str) -> Error {
        Error::ConfigError{msg: String::from(err)}
    }
}

impl From<String> for Error {
    fn from(err: String) -> Error {
        Error::ConfigError{msg: err}
    }
}

impl From<UnitParsingError> for Error {
    fn from(err: UnitParsingError) -> Error {Error::UnitError(err)}
}

/// Custom Result for input files
pub type Result<T> = result::Result<T, Error>;

mod interactions;
pub use self::interactions::read_interactions;
