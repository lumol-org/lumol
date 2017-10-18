// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license
use toml::value::{Value, Table};

use lumol::sys::System;
use lumol::units;

use lumol::energy::{PairPotential, PairInteraction, BondPotential};
use lumol::energy::{Harmonic, LennardJones, NullPotential};
use lumol::energy::{Buckingham, BornMayerHuggins, MorsePotential, Gaussian};
use lumol::energy::TableComputation;

use error::{Error, Result};
use {FromToml, FromTomlWithData};
use extract;
use super::read_restriction;
use super::InteractionsInput;

/// Global settings for the pair interactions
struct GlobalInformation<'a> {
    cutoff: Option<&'a Value>,
    tail: Option<bool>,
}

impl<'a> GlobalInformation<'a> {
    fn read(config: &Table) -> Result<GlobalInformation> {
        match config.get("global") {
            Some(global) => {
                let global = try!(global.as_table()
                    .ok_or(Error::from("'global' section must be a table")));
                let cutoff = global.get("cutoff");
                let tail = if let Some(tail) = global.get("tail_correction") {
                    let tail = try!(tail.as_bool()
                        .ok_or(Error::from("The 'tail_correction' section must be a boolean \
                                            value")));
                    Some(tail)
                } else {
                    None
                };
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
    pub(crate) fn read_pairs(&self, system: &mut System) -> Result<()> {
        let pairs = match self.config.get("pairs") {
            Some(pairs) => pairs,
            None => return Ok(()),
        };

        let pairs = try!(pairs.as_array()
            .ok_or(Error::from("The 'pairs' section must be an array")));

        for pair in pairs {
            let pair = try!(pair.as_table()
                .ok_or(Error::from("Pair potential entry must be a table")));

            let atoms = try!(extract::slice("atoms", pair, "pair potential"));
            if atoms.len() != 2 {
                return Err(Error::from(format!("Wrong size for 'atoms' array in pair \
                                                potential. Should be 2, is {}",
                                               atoms.len())));
            }

            let a = try!(atoms[0]
                .as_str()
                .ok_or(Error::from("The first atom name is not a string in pair potential")));
            let b = try!(atoms[1]
                .as_str()
                .ok_or(Error::from("The second atom name is not a string in pair potential")));

            let potential = try!(read_pair_potential(pair));
            let potential = if let Some(computation) = pair.get("computation") {
                let computation = try!(computation.as_table()
                    .ok_or(Error::from("'computation' section must be a table")));
                try!(read_pair_computation(computation, potential))
            } else {
                potential
            };

            let global = try!(GlobalInformation::read(&self.config));
            let cutoff = match pair.get("cutoff") {
                Some(cutoff) => cutoff,
                None => {
                    try!(global.cutoff
                        .as_ref()
                        .ok_or(Error::from("Missing 'cutoff' value for pair potential")))
                }
            };

            let mut interaction = match *cutoff {
                Value::String(ref cutoff) => {
                    let cutoff = try!(units::from_str(cutoff));
                    PairInteraction::new(potential, cutoff)
                }
                Value::Table(ref table) => {
                    let shifted = try!(table.get("shifted")
                        .ok_or(Error::from("'cutoff' table can only contain 'shifted' key")));
                    let cutoff = try!(shifted.as_str()
                        .ok_or(Error::from("'cutoff.shifted' value must be a string")));
                    let cutoff = try!(units::from_str(cutoff));
                    PairInteraction::shifted(potential, cutoff)
                }
                _ => return Err(Error::from("'cutoff' must be a string or a table")),
            };

            let tail = match pair.get("tail_correction") {
                Some(tail) => {
                    Some(try!(tail.as_bool()
                        .ok_or(Error::from("The 'tail_correction' section must be a boolean \
                                            value"))))
                }
                None => global.tail,
            };

            if let Some(use_tail) = tail {
                if use_tail {
                    interaction.enable_tail_corrections()
                }
            }

            if let Some(restriction) = try!(read_restriction(pair)) {
                interaction.set_restriction(restriction);
            }

            system.add_pair_potential(a, b, interaction);
        }
        Ok(())
    }

    /// Read the "bonds" section from the potential configuration.
    pub(crate) fn read_bonds(&self, system: &mut System) -> Result<()> {
        let bonds = match self.config.get("bonds") {
            Some(bonds) => bonds,
            None => return Ok(()),
        };

        let bonds = try!(bonds.as_array()
            .ok_or(Error::from("The 'bonds' section must be an array")));

        for bond in bonds {
            let bond = try!(bond.as_table()
                .ok_or(Error::from("Bond potential entry must be a table")));

            let atoms = try!(extract::slice("atoms", bond, "bond potential"));
            if atoms.len() != 2 {
                return Err(Error::from(format!("Wrong size for 'atoms' array in bond \
                                                potential. Should be 2, is {}",
                                               atoms.len())));
            }

            let a = try!(atoms[0]
                .as_str()
                .ok_or(Error::from("The first atom name is not a string in pair potential")));
            let b = try!(atoms[1]
                .as_str()
                .ok_or(Error::from("The second atom name is not a string in pair potential")));

            let potential = try!(read_bond_potential(bond));
            system.add_bond_potential(a, b, potential);
        }
        Ok(())
    }
}

fn read_pair_potential(pair: &Table) -> Result<Box<PairPotential>> {
    const KEYWORDS: &'static [&'static str] =
        &["restriction", "computation", "atoms", "cutoff", "tail_correction"];

    let potentials = pair.keys()
        .cloned()
        .filter(|key| !KEYWORDS.contains(&key.as_ref()))
        .collect::<Vec<_>>();

    if potentials.is_empty() {
        return Err(Error::from("Missing potential type in pair potential"));
    }

    if potentials.len() > 1 {
        return Err(Error::from(format!("Got more than one potential type in pair potential: {}",
                                       potentials.join(" and "))));
    }

    let key = &*potentials[0];
    if let Value::Table(ref table) = pair[key] {
        match key {
            "null" => Ok(Box::new(try!(NullPotential::from_toml(table)))),
            "harmonic" => Ok(Box::new(try!(Harmonic::from_toml(table)))),
            "lj" => Ok(Box::new(try!(LennardJones::from_toml(table)))),
            "buckingham" => Ok(Box::new(try!(Buckingham::from_toml(table)))),
            "born" => Ok(Box::new(try!(BornMayerHuggins::from_toml(table)))),
            "morse" => Ok(Box::new(try!(MorsePotential::from_toml(table)))),
            "gaussian" => Ok(Box::new(try!(Gaussian::from_toml(table)))),
            other => Err(Error::from(format!("Unknown potential type '{}'", other))),
        }
    } else {
        Err(Error::from(format!("'{}' potential must be a table", key)))
    }
}

fn read_bond_potential(pair: &Table) -> Result<Box<BondPotential>> {
    let potentials = pair.keys()
        .cloned()
        .filter(|k| k != "atoms")
        .collect::<Vec<_>>();

    if potentials.is_empty() {
        return Err(Error::from("Missing potential type in bond potential"));
    }

    if potentials.len() > 1 {
        return Err(Error::from(format!("Got more than one potential type in bond potential: {}",
                                       potentials.join(" and "))));
    }

    let key = &*potentials[0];
    if let Value::Table(ref table) = pair[key] {
        match key {
            "null" => Ok(Box::new(try!(NullPotential::from_toml(table)))),
            "harmonic" => Ok(Box::new(try!(Harmonic::from_toml(table)))),
            "morse" => Ok(Box::new(try!(MorsePotential::from_toml(table)))),
            other => Err(Error::from(format!("Unknown potential type '{}'", other))),
        }
    } else {
        Err(Error::from(format!("'{}' potential must be a table", key)))
    }
}

/// ***************************************************************************

fn read_pair_computation(computation: &Table,
                         potential: Box<PairPotential>)
                         -> Result<Box<PairPotential>> {
    if computation.keys().len() != 1 {
        return Err(Error::from("Missing computation type in computation table"));
    }

    match computation.keys().map(|s| s.as_ref()).next() {
        Some("table") => Ok(Box::new(try!(TableComputation::from_toml(computation, potential)))),
        Some(other) => Err(Error::from(format!("Unknown computation type '{}'", other))),
        None => unreachable!(),
    }
}
