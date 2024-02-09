// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license
use toml::value::{Table, Value};

use lumol_core::System;
use lumol_core::units;

use lumol_core::energy::{BondPotential, PairInteraction, PairPotential};
use lumol_core::energy::{BornMayerHuggins, Buckingham, Gaussian, Morse};
use lumol_core::energy::{Harmonic, LennardJones, NullPotential, Mie};
use lumol_core::energy::TableComputation;

use super::read_restriction;
use crate::{Error, InteractionsInput, FromToml, FromTomlWithData};
use crate::extract;

/// Global settings for the pair interactions
struct GlobalInformation<'a> {
    cutoff: Option<&'a Value>,
    tail: Option<bool>,
}

impl GlobalInformation<'_> {
    fn read(config: &Table) -> Result<GlobalInformation<'_>, Error> {
        match config.get("global") {
            Some(global) => {
                let global = global.as_table().ok_or(
                    Error::from("'global' section must be a table")
                )?;

                let cutoff = global.get("cutoff");
                let tail = global.get("tail_correction")
                    .map(|tail| {
                        tail.as_bool().ok_or(
                            Error::from("the 'tail_correction' section must be a boolean value")
                        )
                    })
                    .map_or(Ok(None), |tail| tail.map(Some))?;

                Ok(GlobalInformation {
                    cutoff: cutoff,
                    tail: tail,
                })
            }
            None => {
                Ok(GlobalInformation {
                    cutoff: None,
                    tail: None,
                })
            }
        }
    }
}

impl InteractionsInput {
    /// Read the "pairs" section from the potential configuration.
    pub(crate) fn read_pairs(&self, system: &mut System) -> Result<(), Error> {
        let Some(pairs) = self.config.get("pairs") else { return Ok(()) };

        let pairs = pairs.as_table().ok_or(
            Error::from("the 'pairs' section must be a table")
        )?;

        let global = GlobalInformation::read(&self.config)?;

        for (key, table) in pairs {
            let atoms = key.split('-').collect::<Vec<_>>();
            if atoms.len() != 2 {
                return Err(Error::from(format!(
                    "expected two atoms for pair potential, got {} ({:?})", atoms.len(), atoms
                )));
            }

            let table = table.as_table().ok_or(
                Error::from(format!(
                    "pair potential associated with {key} must be a table"
                ))
            )?;

            let potential = read_pair_potential(table)?;
            let potential = if let Some(computation) = table.get("computation") {
                let computation = computation.as_table().ok_or(
                    Error::from("'computation' section must be a table")
                )?;
                read_pair_computation(computation, potential)?
            } else {
                potential
            };

            let cutoff = match table.get("cutoff") {
                Some(cutoff) => cutoff,
                None => {
                    global.cutoff.as_ref().ok_or(
                        Error::from("missing 'cutoff' value for pair potential")
                    )?
                }
            };

            let mut interaction = match *cutoff {
                Value::String(ref cutoff) => {
                    let cutoff = units::from_str(cutoff)?;
                    PairInteraction::new(potential, cutoff)
                }
                Value::Table(ref table) => {
                    let shifted = table.get("shifted").ok_or(
                        Error::from("'cutoff' table can only contain 'shifted' key")
                    )?;
                    let cutoff = shifted.as_str().ok_or(
                        Error::from("'cutoff.shifted' value must be a string")
                    )?;
                    let cutoff = units::from_str(cutoff)?;
                    PairInteraction::shifted(potential, cutoff)
                }
                _ => return Err(Error::from("'cutoff' must be a string or a table")),
            };

            let tail = table.get("tail_correction")
                .map(|tail| {
                    tail.as_bool().ok_or(Error::from(
                        "the 'tail_correction' section must be a boolean value"
                    ))
                })
                .map_or(Ok(global.tail), |tail| tail.map(Some))?;

            if let Some(use_tail) = tail {
                if use_tail {
                    interaction.enable_tail_corrections();
                }
            }

            if let Some(restriction) = read_restriction(table)? {
                interaction.set_restriction(restriction);
            }

            system.set_pair_potential((atoms[0], atoms[1]), interaction);
        }
        Ok(())
    }

    /// Read the "bonds" section from the potential configuration.
    pub(crate) fn read_bonds(&self, system: &mut System) -> Result<(), Error> {
        let Some(bonds) = self.config.get("bonds") else { return Ok(()) };

        let bonds = bonds.as_table().ok_or(
            Error::from("the 'pairs' section must be a table")
        )?;

        for (key, table) in bonds {
            let atoms = key.split('-').collect::<Vec<_>>();
            if atoms.len() != 2 {
                return Err(Error::from(format!(
                    "expected two atoms for bond potential, got {} ({:?})", atoms.len(), atoms
                )));
            }

            let table = table.as_table().ok_or(
                Error::from(format!(
                    "bond potential associated with {key} must be a table"
                ))
            )?;

            let potential = read_bond_potential(table)?;
            system.set_bond_potential((atoms[0], atoms[1]), potential);
        }
        Ok(())
    }
}

fn read_pair_potential(table: &Table) -> Result<Box<dyn PairPotential>, Error> {
    match extract::typ(table, "pair potential")? {
        "null" => Ok(Box::new(NullPotential::from_toml(table)?)),
        "harmonic" => Ok(Box::new(Harmonic::from_toml(table)?)),
        "lj" => Ok(Box::new(LennardJones::from_toml(table)?)),
        "buckingham" => Ok(Box::new(Buckingham::from_toml(table)?)),
        "born" => Ok(Box::new(BornMayerHuggins::from_toml(table)?)),
        "morse" => Ok(Box::new(Morse::from_toml(table)?)),
        "gaussian" => Ok(Box::new(Gaussian::from_toml(table)?)),
        "mie" => Ok(Box::new(Mie::from_toml(table)?)),
        other => Err(Error::from(format!("unknown potential type '{other}'"))),
    }
}

fn read_bond_potential(table: &Table) -> Result<Box<dyn BondPotential>, Error> {
    match extract::typ(table, "bond potential")? {
        "null" => Ok(Box::new(NullPotential::from_toml(table)?)),
        "harmonic" => Ok(Box::new(Harmonic::from_toml(table)?)),
        "morse" => Ok(Box::new(Morse::from_toml(table)?)),
        other => Err(Error::from(format!("unknown potential type '{other}'"))),
    }
}

fn read_pair_computation(computation: &Table, potential: Box<dyn PairPotential>) -> Result<Box<dyn PairPotential>, Error> {
    if computation.keys().len() != 1 {
        return Err(Error::from("Missing computation type in computation table"));
    }

    match computation.keys().map(|s| s.as_ref()).next() {
        Some("table") => Ok(Box::new(TableComputation::from_toml(computation, potential)?)),
        Some(other) => Err(Error::from(format!("Unknown computation type '{other}'"))),
        None => unreachable!(),
    }
}
