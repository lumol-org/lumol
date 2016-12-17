// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license
use toml::{Value, Table};

use lumol::sys::System;
use lumol::energy::{Harmonic, CosineHarmonic, Torsion, NullPotential};
use lumol::energy::{AnglePotential, DihedralPotential};

use error::{Error, Result};
use FromToml;
use extract;
use super::InteractionsInput;

impl InteractionsInput {
    /// Read the "angles" section from the potential configuration. This is an
    /// internal function, public because of the code organization.
    // TODO: use restricted privacy here
    #[doc(hidden)]
    pub fn read_angles(&self, system: &mut System) -> Result<()> {
        let angles = match self.config.get("angles") {
            Some(angles) => angles,
            None => return Ok(())
        };

        let angles = try!(angles.as_slice().ok_or(
            Error::from("The 'angles' section must be an array")
        ));

        for angle in angles {
            let angle = try!(angle.as_table().ok_or(
                Error::from("Angle potential entry must be a table")
            ));

            let atoms = try!(extract::slice("atoms", angle, "angle potential"));
            if atoms.len() != 3 {
                return Err(Error::from(
                    format!("Wrong size for 'atoms' array in angle potential. Should be 3, is {}", atoms.len())
                ));
            }

            let a = try!(atoms[0].as_str().ok_or(Error::from("The first atom name is not a string in angle potential")));
            let b = try!(atoms[1].as_str().ok_or(Error::from("The second atom name is not a string in angle potential")));
            let c = try!(atoms[2].as_str().ok_or(Error::from("The third atom name is not a string in angle potential")));

            let potential = try!(read_angle_potential(angle));
            system.interactions_mut().add_angle(a, b, c, potential);
        }
        Ok(())
    }

    /// Read the "dihedrals" section from the potential configuration. This is
    /// an internal function, public because of the code organization.
    // TODO: use restricted privacy here
    #[doc(hidden)]
    pub fn read_dihedrals(&self, system: &mut System) -> Result<()> {
        let dihedrals = match self.config.get("dihedrals") {
            Some(dihedrals) => dihedrals,
            None => return Ok(())
        };

        let dihedrals = try!(dihedrals.as_slice().ok_or(
            Error::from("The 'dihedrals' section must be an array")
        ));

        for dihedral in dihedrals {
            let dihedral = try!(dihedral.as_table().ok_or(
                Error::from("dihedral potential entry must be a table")
            ));

            let atoms = try!(extract::slice("atoms", dihedral, "dihedral potential"));
            if atoms.len() != 4 {
                return Err(Error::from(
                    format!("Wrong size for 'atoms' array in dihedral potential. Should be 4, is {}", atoms.len())
                ));
            }

            let a = try!(atoms[0].as_str().ok_or(Error::from("The first atom name is not a string in dihedral potential")));
            let b = try!(atoms[1].as_str().ok_or(Error::from("The second atom name is not a string in dihedral potential")));
            let c = try!(atoms[2].as_str().ok_or(Error::from("The third atom name is not a string in dihedral potential")));
            let d = try!(atoms[3].as_str().ok_or(Error::from("The fourth atom name is not a string in dihedral potential")));

            let potential = try!(read_dihedral_potential(dihedral));
            system.interactions_mut().add_dihedral(a, b, c, d, potential);
        }
        Ok(())
    }
}

fn read_angle_potential(angle: &Table) -> Result<Box<AnglePotential>> {
    let potentials = angle.keys().cloned()
                    .filter(|key| key != "atoms")
                    .collect::<Vec<_>>();

    if potentials.is_empty() {
        return Err(Error::from(
            "Missing potential type in angle potential"
        ));
    }

    if potentials.len() > 1 {
        return Err(Error::from(
            format!("Got more than one potential type in angle potential: {}", potentials.join(" and "))
        ));
    }

    let key = &*potentials[0];
    if let Value::Table(ref table) = angle[key] {
        match key {
            "null" => Ok(Box::new(try!(NullPotential::from_toml(table)))),
            "harmonic" => Ok(Box::new(try!(Harmonic::from_toml(table)))),
            "cosine-harmonic" => Ok(Box::new(try!(CosineHarmonic::from_toml(table)))),
            other => Err(
                Error::from(format!("Unknown potential type '{}'", other))
            ),
        }
    } else {
        Err(
            Error::from(format!("'{}' potential must be a table", key))
        )
    }
}

fn read_dihedral_potential(dihedral: &Table) -> Result<Box<DihedralPotential>> {
    let potentials = dihedral.keys().cloned()
                    .filter(|key| key != "atoms")
                    .collect::<Vec<_>>();

    if potentials.is_empty() {
        return Err(Error::from(
            "Missing potential type in dihedral potential"
        ));
    }

    if potentials.len() > 1 {
        return Err(Error::from(
            format!("Got more than one potential type in dihedral potential: {}", potentials.join(" and "))
        ));
    }

    let key = &*potentials[0];
    if let Value::Table(ref table) = dihedral[key] {
        match key {
            "null" => Ok(Box::new(try!(NullPotential::from_toml(table)))),
            "harmonic" => Ok(Box::new(try!(Harmonic::from_toml(table)))),
            "cosine-harmonic" => Ok(Box::new(try!(CosineHarmonic::from_toml(table)))),
            "torsion" => Ok(Box::new(try!(Torsion::from_toml(table)))),
            other => Err(
                Error::from(format!("Unknown potential type '{}'", other))
            ),
        }
    } else {
        Err(
            Error::from(format!("'{}' potential must be a table", key))
        )
    }
}
