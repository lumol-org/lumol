// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license
use toml::{Parser, Table, Value};

use std::io::prelude::*;
use std::fs::File;
use std::path::Path;

use lumol::system::System;
use lumol::potentials::{PairPotential, PairRestriction};

use {Error, Result};
use validate;
use error::toml_error_to_string;

mod toml;
mod pairs;
mod angles;
mod coulomb;

use self::pairs::{read_pairs, read_bonds};
use self::angles::{read_angles, read_dihedrals};
use self::coulomb::{read_coulomb, set_charges};

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
    let config = try!(parser.parse().ok_or(
        Error::TOML(toml_error_to_string(&parser))
    ));

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

        let cutoff = if let Some(global) = config.get("global") {
            let global = try!(global.as_table().ok_or(Error::from(
                "'global' section must be a table"
            )));
            global.get("cutoff")
        } else {
            None
        };

        try!(read_pairs(system, pairs, cutoff));
    }

    if let Some(bonds) = config.get("bonds") {
        let bonds = try!(bonds.as_slice().ok_or(
            Error::from("The 'bonds' section must be an array")
        ));
        try!(read_bonds(system, bonds));
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
    let restriction = config.get("restriction");
    if restriction.is_none() {
        // No restriction found
        return Ok(None);
    };

    match restriction.expect("Unreachable").clone() {
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
    use lumol::system::System;
    use read_interactions;
    use testing::bad_inputs;

    #[test]
    fn bad_input() {
        for path in bad_inputs("interactions", "generic") {
            let mut system = System::new();
            assert!(read_interactions(&mut system, path).is_err());
        }
    }
}
