// Cymbalum, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license
use toml::{Table, Value};

use input::error::{Error, Result};
use super::FromToml;
use super::read_restriction;

use system::System;
use potentials::{Wolf, Ewald, CoulombicPotential};

pub fn read_coulomb(system: &mut System, coulomb: &Table) -> Result<()> {
    let solvers = coulomb.keys().cloned()
                         .filter(|key| key != "restriction")
                         .collect::<Vec<_>>();

    if solvers.len() != 1 {
        return Err(Error::from(
            format!("Got more than one coulombic solver: {:?}", solvers)
        ));
    }

    let key = &*solvers[0];
    if let Value::Table(ref table) = coulomb[key] {
        let mut potential: Box<CoulombicPotential> = match key {
            "wolf" => Box::new(try!(Wolf::from_toml(table))),
            "ewald" => Box::new(try!(Ewald::from_toml(table))),
            other => {
                return Err(Error::from(format!("Unknown coulomb solver '{}'", other)))
            },
        };

        match try!(read_restriction(coulomb)) {
            Some(restriction) => {
                potential.set_restriction(restriction);
            }
            None => {}
        }

        system.set_coulomb_interaction(potential);
        Ok(())
    } else {
        Err(
            Error::from(format!("Coulombic solver '{}' must be a table", key))
        )
    }
}

pub fn set_charges(system: &mut System, charges: &Table) -> Result<()> {
    for (name, charge) in charges.iter() {
        let charge = match *charge {
            Value::Integer(val) => val as f64,
            Value::Float(val) => val,
            _ => {
                return Err(Error::from("Charges must be numbers"));
            }
        };

        let mut nchanged = 0;
        for particle in system.iter_mut() {
            if particle.name() == name {
                particle.charge = charge;
                nchanged += 1;
            }
        }

        if nchanged == 0 {
            warn!("No particle with name '{}' was found while setting the charges", name);
        } else {
            info!("Charge set to {} for {} {} particles", charge, nchanged, name);
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use input::read_interactions;
    use input::interactions::testing::bad_interactions;
    use system::{Particle, System};
    use std::path::Path;

    #[test]
    fn ewald() {
        let data_root = Path::new(file!()).parent().unwrap().join("data");
        let mut system = System::new();
        system.add_particle(Particle::new("A"));
        system.add_particle(Particle::new("B"));

        read_interactions(&mut system, data_root.join("ewald.toml")).unwrap();
        assert!(system.coulomb_potential().is_some());

        assert_eq!(system[0].charge, -8.0);
        assert_eq!(system[1].charge, 3.0);
    }

    #[test]
    fn wolf() {
        let data_root = Path::new(file!()).parent().unwrap().join("data");
        let mut system = System::new();
        system.add_particle(Particle::new("A"));
        system.add_particle(Particle::new("B"));

        read_interactions(&mut system, data_root.join("wolf.toml")).unwrap();
        assert!(system.coulomb_potential().is_some());

        assert_eq!(system[0].charge, -2.0);
        assert_eq!(system[1].charge, 2.0);
    }

    #[test]
    fn bad_coulomb() {
        for path in bad_interactions("coulomb") {
            let mut system = System::new();
            assert!(read_interactions(&mut system, path).is_err());
        }
    }
}
