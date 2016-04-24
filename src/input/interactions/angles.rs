// Cymbalum, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license
use toml::{Value, Table};

use system::System;
use input::error::{Error, Result};
use super::FromToml;

use potentials::{Harmonic, CosineHarmonic, Torsion, NullPotential};
use potentials::{AnglePotential, DihedralPotential};

pub fn read_angles(system: &mut System, angles: &[Value]) -> Result<()> {
    for angle in angles {
        let angle = try!(angle.as_table().ok_or(
            Error::from("Angle potential entry must be a table")
        ));

        let atoms = try!(angle.get("atoms").ok_or(
            Error::from("Missing 'atoms' section in angle potential")
        ));

        let atoms = try!(atoms.as_slice().ok_or(
            Error::from("'atoms' section must be an array")
        ));

        if atoms.len() != 3 {
            return Err(Error::from(
                format!("Wrong size for 'atoms' section in angle potentials. Should be 3, is {}", atoms.len())
            ));
        }

        let a = try!(atoms[0].as_str().ok_or(Error::from("The first atom name is not a string in angle potential")));
        let b = try!(atoms[1].as_str().ok_or(Error::from("The second atom name is not a string in angle potential")));
        let c = try!(atoms[2].as_str().ok_or(Error::from("The third atom name is not a string in angle potential")));

        let potential = try!(read_angle_potential(angle));
        system.add_angle_interaction(a, b, c, potential);
    }
    Ok(())
}


fn read_angle_potential(angle: &Table) -> Result<Box<AnglePotential>> {
    let potentials = angle.keys().cloned()
                    .filter(|key| key != "atoms")
                    .collect::<Vec<_>>();

    if potentials.len() != 1 {
        return Err(Error::from(
            format!("Got more than one potential type: {:?}", potentials)
        ));
    }

    let key = &*potentials[0];
    if let Value::Table(ref table) = angle[key] {
        match key {
            "null" => Ok(Box::new(NullPotential::from_toml(table).unwrap())),
            "harmonic" => Ok(Box::new(try!(Harmonic::from_toml(table)))),
            "cosine-harmonic" => Ok(Box::new(try!(CosineHarmonic::from_toml(table)))),
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

pub fn read_dihedrals(system: &mut System, dihedrals: &[Value]) -> Result<()> {
    for dihedral in dihedrals {
        let dihedral = try!(dihedral.as_table().ok_or(
            Error::from("Dihedral angle potential entry must be a table")
        ));

        let atoms = try!(dihedral.get("atoms").ok_or(
            Error::from("Missing 'atoms' section in dihedral angle potential")
        ));

        let atoms = try!(atoms.as_slice().ok_or(
            Error::from("'atoms' section must be an array")
        ));

        if atoms.len() != 4 {
            return Err(Error::from(
                format!("Wrong size for 'atoms' section in dihedral angle potentials. Should be 4, is {}", atoms.len())
            ));
        }

        let a = try!(atoms[0].as_str().ok_or(Error::from("The first atom name is not a string in dihedral angle potential")));
        let b = try!(atoms[1].as_str().ok_or(Error::from("The second atom name is not a string in dihedral angle potential")));
        let c = try!(atoms[2].as_str().ok_or(Error::from("The third atom name is not a string in dihedral angle potential")));
        let d = try!(atoms[3].as_str().ok_or(Error::from("The fourth atom name is not a string in dihedral angle potential")));

        let potential = try!(read_dihedral_potential(dihedral));
        system.add_dihedral_interaction(a, b, c, d, potential);
    }
    Ok(())
}


fn read_dihedral_potential(dihedral: &Table) -> Result<Box<DihedralPotential>> {
    let potentials = dihedral.keys().cloned()
                    .filter(|key| key != "atoms")
                    .collect::<Vec<_>>();

    if potentials.len() != 1 {
        return Err(Error::from(
            format!("Got more than one potential type: {:?}", potentials)
        ));
    }

    let key = &*potentials[0];
    if let Value::Table(ref table) = dihedral[key] {
        match key {
            "null" => Ok(Box::new(NullPotential::from_toml(table).unwrap())),
            "harmonic" => Ok(Box::new(try!(Harmonic::from_toml(table)))),
            "cosine-harmonic" => Ok(Box::new(try!(CosineHarmonic::from_toml(table)))),
            "torsion" => Ok(Box::new(try!(Torsion::from_toml(table)))),
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

#[cfg(test)]
mod tests {
    use input::read_interactions;
    use input::interactions::testing::bad_interactions;
    use system::{Particle, System};
    use std::path::Path;

    #[test]
    fn angles() {
        let data_root = Path::new(file!()).parent().unwrap().join("data");
        let mut system = System::new();
        system.add_particle(Particle::new("A"));
        system.add_particle(Particle::new("B"));
        system.add_particle(Particle::new("C"));
        let _ = system.add_bond(0, 1);
        let _ = system.add_bond(1, 2);

        read_interactions(&mut system, data_root.join("angles.toml")).unwrap();
        assert_eq!(system.angle_potentials(0, 1, 2).len(), 3);
    }

    #[test]
    fn bad_angles() {
        for path in bad_interactions("angles") {
            let mut system = System::new();
            assert!(read_interactions(&mut system, path).is_err());
        }
    }

    #[test]
    fn dihedrals() {
        let data_root = Path::new(file!()).parent().unwrap().join("data");
        let mut system = System::new();
        system.add_particle(Particle::new("A"));
        system.add_particle(Particle::new("B"));
        system.add_particle(Particle::new("C"));
        system.add_particle(Particle::new("D"));
        let _ = system.add_bond(0, 1);
        let _ = system.add_bond(1, 2);
        let _ = system.add_bond(2, 3);

        read_interactions(&mut system, data_root.join("dihedrals.toml")).unwrap();
        assert_eq!(system.dihedral_potentials(0, 1, 2, 3).len(), 4);
    }

    #[test]
    fn bad_dihedrals() {
        for path in bad_interactions("dihedrals") {
            let mut system = System::new();
            assert!(read_interactions(&mut system, path).is_err());
        }
    }
}
