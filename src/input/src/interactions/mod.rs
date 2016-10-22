// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license
use toml::{Parser, Table, Value};

use std::io::prelude::*;
use std::fs::File;
use std::path::PathBuf;

use lumol::sys::System;
use lumol::energy::PairRestriction;

use {Error, Result};
use validate;
use error::toml_error_to_string;

mod toml;
mod pairs;
mod angles;
mod coulomb;

/// An interaction input file for Lumol.
pub struct InteractionsInput {
    /// The TOML configuration
    config: Table,
}

impl InteractionsInput {
    /// Read interactions from the TOML formatted file at `path`.
    pub fn new<P: Into<PathBuf>>(path: P) -> Result<InteractionsInput> {
        let path = path.into();
        let mut file = try_io!(File::open(&path), path);
        let mut buffer = String::new();
        let _ = try_io!(file.read_to_string(&mut buffer), path);
        return InteractionsInput::from_string(&buffer);
    }

    /// Read the interactions from a TOML formatted string.
    pub fn from_string(string: &str) -> Result<InteractionsInput> {
        let mut parser = Parser::new(string);
        let config = try!(parser.parse().ok_or(
            Error::TOML(toml_error_to_string(&parser))
        ));

        try!(validate(&config));
        return InteractionsInput::from_toml(config);
    }

    /// Read the interactions from a TOML table. This is an internal function,
    /// public because of the code organization.
    // TODO: use restricted privacy here
    pub fn from_toml(config: Table) -> Result<InteractionsInput> {
        Ok(InteractionsInput{
            config: config
        })
    }

    /// Read the interactions from this input into the `system`.
    pub fn read(&self, system: &mut System) -> Result<()> {
        try!(self.read_pairs(system));
        try!(self.read_bonds(system));
        try!(self.read_angles(system));
        try!(self.read_dihedrals(system));
        try!(self.read_coulomb(system));
        try!(self.read_charges(system));
        Ok(())
    }
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
    use lumol::sys::System;
    use InteractionsInput;
    use testing::bad_inputs;

    #[test]
    fn bad_input() {
        let mut system = System::new();
        for path in bad_inputs("interactions", "generic") {
            assert!(
                InteractionsInput::new(path)
                .and_then(|input| input.read(&mut system))
                .is_err()
            );
        }
    }
}
