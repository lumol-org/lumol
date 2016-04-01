// Cymbalum, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license
use toml::{Parser, ParserError, Table};

use std::io::prelude::*;
use std::io;
use std::result;
use std::fs::File;
use std::path::Path;

use system::System;
use units::UnitParsingError;
use potentials::PairPotential;

mod toml;
mod pairs;
mod angles;
#[cfg(test)]
pub mod testing;

use self::pairs::{TwoBody, read_2body};
use self::angles::{read_angles, read_dihedrals};

#[derive(Debug)]
/// Possible causes of error when reading potential files
pub enum Error {
    /// Error in the TOML input file
    TOML(Vec<ParserError>),
    /// IO error
    File(io::Error),
    /// File content error: missing sections, bad data types
    Config{
        /// Error message
        msg: String,
    },
    /// Unit parsing error
    UnitParsing(UnitParsingError),
}

impl From<io::Error> for Error {
    fn from(err: io::Error) -> Error {Error::File(err)}
}

impl From<Vec<ParserError>> for Error {
    fn from(errs: Vec<ParserError>) -> Error {
        Error::TOML(errs)
    }
}

impl<'a> From<&'a str> for Error {
    fn from(err: &'a str) -> Error {
        Error::Config{msg: String::from(err)}
    }
}

impl From<String> for Error {
    fn from(err: String) -> Error {
        Error::Config{msg: err}
    }
}

impl From<UnitParsingError> for Error {
    fn from(err: UnitParsingError) -> Error {Error::UnitParsing(err)}
}

/// Custom `Result` for input files
pub type Result<T> = result::Result<T, Error>;

/// Convert a TOML table to a Rust type. This is the trait to implement in order
/// to use the input files.
pub trait FromToml: Sized {
    /// Do the conversion from `table` to Self.
    fn from_toml(table: &Table) -> Result<Self>;
}

/// Convert a TOML table and a PairPotential to a Rust type. This is intended
/// to be used by potential computation mainly.
pub trait FromTomlWithPairs where Self: Sized {
    /// Do the conversion from `table` and `potential` to Self.
    fn from_toml(table: &Table, potential: Box<PairPotential>) -> Result<Self>;
}

/// Read interactions from the TOML file at `path`, and add them to the
/// `system`. For a full documentation of the input files syntax, see the user
/// manual.
pub fn read_interactions<P: AsRef<Path>>(system: &mut System, path: P) -> Result<()> {
    let mut file = try!(File::open(path));
    let mut buffer = String::new();
    let _ = try!(file.read_to_string(&mut buffer));
    return read_interactions_string(system, &buffer);
}


/// This is the same as `read_interactions`, but directly read a TOML formated
/// string.
pub fn read_interactions_string(system: &mut System, string: &str) -> Result<()> {
    let mut parser = Parser::new(string);
    let config = match parser.parse() {
        Some(config) => config,
        None => {
            return Err(Error::TOML(parser.errors))
        }
    };

    try!(validate(&config));
    if let Some(pairs) = config.get("pairs") {
        let pairs = try!(pairs.as_slice().ok_or(
            Error::from("The 'pairs' section must be an array")
        ));
        try!(read_2body(system, pairs, TwoBody::Pairs));
    }

    if let Some(bonds) = config.get("bonds") {
        let bonds = try!(bonds.as_slice().ok_or(
            Error::from("The 'bonds' section must be an array")
        ));
        try!(read_2body(system, bonds, TwoBody::Bonds));
    }

    if let Some(angles) = config.get("angles") {
        let angles = try!(angles.as_slice().ok_or(
            Error::from("The 'angles' section must be an array")
        ));
        try!(read_angles(system, angles));
    }

    if let Some(dihedrals) = config.get("dihedrals") {
        let dihedrals = try!(dihedrals.as_slice().ok_or(
            Error::from("The 'dihedrals' section must be an array")
        ));
        try!(read_dihedrals(system, dihedrals));
    }

    Ok(())
}

fn validate(config: &Table) -> Result<()> {
    let input = try!(config.get("input").ok_or(
        Error::from("Missing 'input' table")
    ));

    let version = try!(input.lookup("potentials.version").ok_or(
        Error::from("Missing 'potentials.version' key in 'input' table")
    ));

    let version = try!(version.as_integer().ok_or(
        Error::from("'input.potentials.version' must be an integer")
    ));

    if version != 1 {
        return Err(Error::from(
            format!("Only version 1 of input can be read, got {}", version)
        ))
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use system::System;
    use input::read_interactions;
    use input::interactions::testing::bad_interactions;

    #[test]
    fn bad_input() {
        for path in bad_interactions("generic") {
            let mut system = System::new();
            assert!(read_interactions(&mut system, path).is_err());
        }
    }
}
