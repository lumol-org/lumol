// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license
use toml::de::from_str as parse;
use toml::value::{Table, Value};

use std::fs::File;
use std::io::prelude::*;
use std::path::PathBuf;

use lumol_core::energy::PairRestriction;
use lumol_core::System;

use crate::Error;
use crate::validate;

mod potentials;
mod pairs;
mod angles;
mod coulomb;

/// Input file for reading interactions
pub struct InteractionsInput {
    /// The TOML configuration
    config: Table,
}

impl InteractionsInput {
    /// Read interactions from the TOML formatted file at `path`.
    pub fn new<P: Into<PathBuf>>(path: P) -> Result<InteractionsInput, Error> {
        let path = path.into();
        let mut file = try_io!(File::open(&path), path);
        let mut buffer = String::new();
        let _ = try_io!(file.read_to_string(&mut buffer), path);
        return InteractionsInput::from_str(&buffer);
    }

    /// Read the interactions from a TOML formatted string.
    #[allow(clippy::should_implement_trait)]
    pub fn from_str(string: &str) -> Result<InteractionsInput, Error> {
        let config = parse(string).map_err(|err| Error::TOML(Box::new(err)))?;
        validate(&config)?;
        return Ok(InteractionsInput::from_toml(config));
    }

    /// Read the interactions from a TOML table.
    pub(crate) fn from_toml(config: Table) -> InteractionsInput {
        InteractionsInput {
            config: config
        }
    }

    /// Read the interactions from this input into the `system`.
    pub fn read(&self, system: &mut System) -> Result<(), Error> {
        self.read_pairs(system)?;
        self.read_bonds(system)?;
        self.read_angles(system)?;
        self.read_dihedrals(system)?;
        // charges must be read before coulomb
        self.read_charges(system)?;
        self.read_coulomb(system)?;
        Ok(())
    }
}

fn read_restriction(config: &Table) -> Result<Option<PairRestriction>, Error> {
    let restriction = config.get("restriction");
    if restriction.is_none() {
        // No restriction found
        return Ok(None);
    };

    match restriction.expect("Unreachable").clone() {
        Value::String(name) => {
            match &*name {
                "none" => Ok(Some(PairRestriction::None)),
                "intramolecular" | "IntraMolecular" | "intra-molecular" => {
                    Ok(Some(PairRestriction::IntraMolecular))
                }
                "intermolecular" | "InterMolecular" | "inter-molecular" => {
                    Ok(Some(PairRestriction::InterMolecular))
                }
                "exclude12" => Ok(Some(PairRestriction::Exclude12)),
                "exclude13" => Ok(Some(PairRestriction::Exclude13)),
                "exclude14" => Ok(Some(PairRestriction::Exclude14)),
                "scale14" => Err(Error::from("'scale14' restriction must be a table")),
                other => Err(Error::from(format!("Unknown restriction '{other}'"))),
            }
        }
        Value::Table(ref restriction) => {
            if restriction.keys().len() != 1 || restriction.get("scale14").is_none() {
                return Err(Error::from("Restriction table must be 'scale14'"));
            }

            let scale = restriction["scale14"].as_float().ok_or(
                Error::from("'scale14' parameter must be a float")
            )?;

            Ok(Some(PairRestriction::Scale14(scale)))
        }
        _ => Err(Error::from("Restriction must be a table or a string")),
    }
}
