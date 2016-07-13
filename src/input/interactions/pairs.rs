// Cymbalum, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license
use toml::{Value, Table};

use system::System;
use input::error::{Error, Result};
use input::FromToml;
use super::{FromTomlWithPairs, read_restriction};

use potentials::{PairPotential, PairInteraction, BondPotential};
use potentials::{Harmonic, LennardJones, NullPotential};
use potentials::{TableComputation, CutoffComputation};

/// Read either the "pairs" section from the configuration.
pub fn read_pairs(system: &mut System, pairs: &[Value]) -> Result<()> {
    for pair in pairs {
        let pair = try!(pair.as_table().ok_or(
            Error::from("Pair potential entry must be a table")
        ));

        let atoms = extract_slice!("atoms", pair as "pair potential");
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

        // TODO: get the real cutoff
        let mut interaction = PairInteraction::new(potential, 10000.0);
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

        let atoms = extract_slice!("atoms", bond as "bond potential");
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
    let potentials = pair.keys().cloned()
                    .filter(|k| k != "restriction" && k != "computation" && k != "atoms")
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
        Some("cutoff") => Ok(
            Box::new(try!(CutoffComputation::from_toml(computation, potential)))
        ),
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
    use input::read_interactions;
    use input::testing::bad_inputs;
    use system::{Particle, System};
    use std::path::Path;

    #[test]
    fn pairs() {
        let data_root = Path::new(file!()).parent().unwrap().join("data");
        let mut system = System::new();
        system.add_particle(Particle::new("A"));
        system.add_particle(Particle::new("B"));

        read_interactions(&mut system, data_root.join("pairs.toml")).unwrap();

        assert_eq!(system.pair_potentials(0, 1).len(), 10);
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
