// Cymbalum, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license
use toml::{Parser, Table, Value};

use std::io::prelude::*;
use std::fs::File;
use std::path::Path;

use system::System;
use potentials::{PairPotential, PairRestriction};

use input::{Error, Result};
use input::validate;
use input::error::toml_error_to_string;

mod toml;
mod pairs;
mod angles;
mod coulomb;

use self::pairs::{TwoBody, read_2body};
use self::angles::{read_angles, read_dihedrals};
use self::coulomb::{read_coulomb, set_charges};

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
// TODO: use restricted privacy for this function
pub fn read_interactions_string(system: &mut System, string: &str) -> Result<()> {
    let mut parser = Parser::new(string);
    let config = match parser.parse() {
        Some(config) => config,
        None => {
            let errors = toml_error_to_string(&parser);
            return Err(Error::TOML(errors));
        }
    };
    try!(validate(&config));
    return read_interactions_toml(system, &config);
}


/// This is the same as `read_interactions`, but directly read a TOML table.
/// The `config` is assumed to be validated by a call to `validate`.
// TODO: use restricted privacy for this function
pub fn read_interactions_toml(system: &mut System, config: &Table) -> Result<()> {
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

    if let Some(coulomb) = config.get("coulomb") {
        let coulomb = try!(coulomb.as_table().ok_or(
            Error::from("The 'coulomb' section must be a table")
        ));
        try!(read_coulomb(system, coulomb));
    }

    if let Some(charges) = config.get("charges") {
        let charges = try!(charges.as_table().ok_or(
            Error::from("The 'charges' section must be a table")
        ));
        try!(set_charges(system, charges));
    }

    Ok(())
}

fn read_restriction(config: &Table) -> Result<Option<PairRestriction>> {
    let restriction = match config.get("restriction") {
        Some(restriction) => restriction,
        None => {return Ok(None)}
    };

    match restriction.clone() {
        Value::String(name) => {
            match &*name {
                "none" => Ok(Some(PairRestriction::None)),
                "intramolecular" | "IntraMolecular" | "intra-molecular"
                    => Ok(Some(PairRestriction::IntraMolecular)),
                "intermolecular" | "InterMolecular" | "inter-molecular"
                    => Ok(Some(PairRestriction::InterMolecular)),
                "exclude12" => Ok(Some(PairRestriction::Exclude12)),
                "exclude13" => Ok(Some(PairRestriction::Exclude13)),
                "exclude14" => Ok(Some(PairRestriction::Exclude14)),
                "scale14" => Err(
                    Error::from("'scale14' restriction must be a table")
                ),
                other => Err(
                    Error::from(format!("Unknown restriction '{}'", other))
                ),
            }
        },
        Value::Table(ref restriction) => {
            if restriction.keys().len() != 1 || restriction.get("scale14").is_none() {
                return Err(Error::from("Restriction table must be 'scale14'"));
            }
            let scale = try!(restriction["scale14"].as_float().ok_or(
                Error::from("'scale14' parameter must be a float")
            ));
            Ok(Some(PairRestriction::Scale14{scaling: scale}))
        }
        _ => Err(Error::from("Restriction must be a table or a string"))
    }
}

#[cfg(test)]
mod tests {
    use system::System;
    use input::read_interactions;
    use input::testing::bad_inputs;

    #[test]
    fn bad_input() {
        for path in bad_inputs("interactions", "generic") {
            let mut system = System::new();
            assert!(read_interactions(&mut system, path).is_err());
        }
    }
}
