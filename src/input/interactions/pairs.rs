// Cymbalum, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license
use toml::{Value, Table};

use system::System;
use super::{Error, Result, FromToml, FromTomlWithPairs};

use potentials::{PairPotential, PairRestriction};
use potentials::{Harmonic, LennardJones, NullPotential};
use potentials::{TableComputation, CutoffComputation};

pub enum TwoBody {
    Pairs,
    Bonds
}

/// Read either the "pairs" or the "bonds" section from the configuration. The
/// `form` enum set where the interactions will be saved.
pub fn read_2body(system: &mut System, pairs: &[Value], form: TwoBody) -> Result<()> {
    for pair in pairs {
        let pair = try!(pair.as_table().ok_or(
            Error::from("Pair potential entry must be a table")
        ));

        let atoms = try!(pair.get("atoms").ok_or(
            Error::from("Missing 'atoms' section in pair potential")
        ));

        let atoms = try!(atoms.as_slice().ok_or(
            Error::from("'atoms' section must be an array")
        ));

        if atoms.len() != 2 {
            return Err(Error::from(
                format!("Wrong size for 'atoms' section in pair potentials. Should be 2, is {}", atoms.len())
            ));
        }

        let a = try!(atoms[0].as_str().ok_or(Error::from("The first atom name is not a string in pair potential")));
        let b = try!(atoms[1].as_str().ok_or(Error::from("The second atom name is not a string in pair potential")));

        let potential = try!(read_pair_potential(pair));
        let potential = if let Some(computation) = pair.get("computation") {
            let computation = try!(computation.as_table().ok_or(
                Error::from("'computation' section must be a table")
            ));
            try!(read_pair_computation(computation, potential))
        } else {
            potential
        };

        match form {
            TwoBody::Pairs => {
                match try!(read_restriction(pair)) {
                    Some(restriction) => {
                        system.add_pair_interaction_with_restriction(a, b, potential, restriction);
                    },
                    None => system.add_pair_interaction(a, b, potential)
                }
            },
            TwoBody::Bonds => {
                system.add_bond_interaction(a, b, potential);
            }
        }
    }
    Ok(())
}


fn read_pair_potential(pair: &Table) -> Result<Box<PairPotential>> {
    let potentials = pair.keys().cloned()
                    .filter(|k| k != "restriction" && k != "computation" && k != "atoms")
                    .collect::<Vec<_>>();

    if potentials.len() != 1 {
        return Err(Error::from(
            format!("Got more than one potential type: {:?}", potentials)
        ));
    }

    let key = &*potentials[0];
    if let Value::Table(ref table) = pair[key] {
        match key {
            "null" => Ok(Box::new(NullPotential::from_toml(table).unwrap())),
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

/******************************************************************************/
fn read_restriction(config: &Table) -> Result<Option<PairRestriction>> {
    let restriction = match config.get("restriction") {
        Some(restriction) => restriction,
        None => {return Ok(None)}
    };

    match restriction {
        &Value::String(ref name) => {
            match &**name {
                "none" => Ok(Some(PairRestriction::None)),
                "intramolecular" => Ok(Some(PairRestriction::IntraMolecular)),
                "intermolecular" => Ok(Some(PairRestriction::InterMolecular)),
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
        &Value::Table(ref restriction) => {
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
    use input::interactions::testing::bad_interactions;
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
        for path in bad_interactions("pairs") {
            let mut system = System::new();
            assert!(read_interactions(&mut system, path).is_err());
        }
    }

    #[test]
    fn bad_bonds() {
        for path in bad_interactions("bonds") {
            let mut system = System::new();
            assert!(read_interactions(&mut system, path).is_err());
        }
    }
}
