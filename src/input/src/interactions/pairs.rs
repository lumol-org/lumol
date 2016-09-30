// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license
use toml::{Value, Table};

use lumol::system::System;
use lumol::units;

use lumol::potentials::{PairPotential, PairInteraction, BondPotential};
use lumol::potentials::{Harmonic, LennardJones, NullPotential};
use lumol::potentials::TableComputation;

use error::{Error, Result};
use {FromToml, FromTomlWithData};
use extract;
use super::read_restriction;

/// Read the "pairs" section from the configuration.
pub fn read_pairs(system: &mut System, pairs: &[Value], global_cutoff: Option<&Value>) -> Result<()> {
    for pair in pairs {
        let pair = try!(pair.as_table().ok_or(
            Error::from("Pair potential entry must be a table")
        ));

        let atoms = try!(extract::slice("atoms", pair, "pair potential"));
        if atoms.len() != 2 {
            return Err(Error::from(
                format!("Wrong size for 'atoms' section in pair potential. Should be 2, is {}", atoms.len())
            ));
        }

        let a = try!(atoms[0].as_str().ok_or(Error::from(
            "The first atom name is not a string in pair potential"
        )));
        let b = try!(atoms[1].as_str().ok_or(Error::from(
            "The second atom name is not a string in pair potential"
        )));

        let potential = try!(read_pair_potential(pair));
        let potential = if let Some(computation) = pair.get("computation") {
            let computation = try!(computation.as_table().ok_or(
                Error::from("'computation' section must be a table")
            ));
            try!(read_pair_computation(computation, potential))
        } else {
            potential
        };

        let cutoff = match pair.get("cutoff") {
            Some(cutoff) => cutoff,
            None => try!(global_cutoff.ok_or(Error::from(
                "Missing 'cutoff' value for pair potential"
            )))
        };

        let mut interaction = match *cutoff {
            Value::String(ref cutoff) => {
                let cutoff = try!(units::from_str(cutoff));
                PairInteraction::new(potential, cutoff)
            }
            Value::Table(ref table) => {
                let shifted = try!(table.get("shifted").ok_or(Error::from(
                    "'cutoff' table can only contain 'shifted' key"
                )));
                let cutoff = try!(shifted.as_str().ok_or(Error::from(
                    "'cutoff.shifted' value must be a string"
                )));
                let cutoff = try!(units::from_str(cutoff));
                PairInteraction::shifted(potential, cutoff)
            }
            _ => return Err(Error::from(
                "'cutoff' must be a string or a table in pair potential"
            ))
        };

        if let Some(restriction) = try!(read_restriction(pair)) {
            interaction.set_restriction(restriction);
        }

        system.interactions_mut().add_pair(a, b, interaction);
    }
    Ok(())
}

/// Read the "bonds" section from the configuration.
pub fn read_bonds(system: &mut System, bonds: &[Value]) -> Result<()> {
    for bond in bonds {
        let bond = try!(bond.as_table().ok_or(
            Error::from("Bond potential entry must be a table")
        ));

        let atoms = try!(extract::slice("atoms", bond, "bond potential"));
        if atoms.len() != 2 {
            return Err(Error::from(
                format!("Wrong size for 'atoms' section in bond potential. Should be 2, is {}", atoms.len())
            ));
        }

        let a = try!(atoms[0].as_str().ok_or(Error::from("The first atom name is not a string in pair potential")));
        let b = try!(atoms[1].as_str().ok_or(Error::from("The second atom name is not a string in pair potential")));

        let potential = try!(read_bond_potential(bond));
        system.interactions_mut().add_bond(a, b, potential);
    }
    Ok(())
}


fn read_pair_potential(pair: &Table) -> Result<Box<PairPotential>> {
    const KEYWORDS: &'static[&'static str] = &["restriction", "computation", "atoms", "cutoff"];

    let potentials = pair.keys().cloned()
                    .filter(|key| !KEYWORDS.contains(&key.as_ref()))
                    .collect::<Vec<_>>();

    if potentials.len() != 1 {
        return Err(Error::from(
            format!("Got more than one potential type: {}", potentials.join(" - "))
        ));
    }

    let key = &*potentials[0];
    if let Value::Table(ref table) = pair[key] {
        match key {
            "null" => Ok(Box::new(try!(NullPotential::from_toml(table)))),
            "harmonic" => Ok(Box::new(try!(Harmonic::from_toml(table)))),
            "lj" | "lennardjones" => Ok(Box::new(try!(LennardJones::from_toml(table)))),
            other => Err(
                Error::from(format!("Unknown potential type '{}'", other))
            ),
        }
    } else {
        Err(
            Error::from(format!("potential '{}' must be a table", key))
        )
    }
}

fn read_bond_potential(pair: &Table) -> Result<Box<BondPotential>> {
    let potentials = pair.keys().cloned()
                    .filter(|k| k != "atoms")
                    .collect::<Vec<_>>();

    if potentials.len() != 1 {
        return Err(Error::from(
            format!("Got more than one potential type: {}", potentials.join(" - "))
        ));
    }

    let key = &*potentials[0];
    if let Value::Table(ref table) = pair[key] {
        match key {
            "null" => Ok(Box::new(try!(NullPotential::from_toml(table)))),
            "harmonic" => Ok(Box::new(try!(Harmonic::from_toml(table)))),
            other => Err(
                Error::from(format!("Unknown potential type '{}'", other))
            ),
        }
    } else {
        Err(
            Error::from(format!("potential '{}' must be a table", key))
        )
    }
}

/******************************************************************************/

fn read_pair_computation(computation: &Table, potential: Box<PairPotential>) -> Result<Box<PairPotential>> {
    if computation.keys().len() != 1 {
        return Err(Error::from("Missing computation type in computation table"));
    }

    match computation.keys().map(|s| s.as_ref()).next() {
        Some("table") => Ok(
            Box::new(try!(TableComputation::from_toml(computation, potential)))
        ),
        Some(other) => Err(
            Error::from(format!("Unknown computation type '{}'", other))
        ),
        None => unreachable!()
    }
}

#[cfg(test)]
mod tests {
    use read_interactions;
    use testing::bad_inputs;
    use lumol::system::{Particle, System};
    use std::path::Path;

    #[test]
    fn pairs() {
        let data_root = Path::new(file!()).parent().unwrap().join("data");
        let mut system = System::new();
        system.add_particle(Particle::new("A"));
        system.add_particle(Particle::new("B"));

        read_interactions(&mut system, data_root.join("pairs.toml")).unwrap();

        assert_eq!(system.pair_potentials(0, 1).len(), 11);
    }

    #[test]
    fn bonds() {
        let data_root = Path::new(file!()).parent().unwrap().join("data");
        let mut system = System::new();
        system.add_particle(Particle::new("A"));
        system.add_particle(Particle::new("B"));
        let _ = system.add_bond(0, 1);

        read_interactions(&mut system, data_root.join("bonds.toml")).unwrap();

        assert_eq!(system.bond_potentials(0, 1).len(), 2);
    }

    #[test]
    fn bad_pairs() {
        for path in bad_inputs("interactions", "pairs") {
            let mut system = System::new();
            assert!(read_interactions(&mut system, path).is_err());
        }
    }

    #[test]
    fn bad_bonds() {
        for path in bad_inputs("interactions", "bonds") {
            let mut system = System::new();
            assert!(read_interactions(&mut system, path).is_err());
        }
    }
}
