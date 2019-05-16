// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license
use toml::value::{Table, Value};

use lumol_core::energy::{AnglePotential, DihedralPotential};
use lumol_core::energy::{CosineHarmonic, Harmonic, Morse, NullPotential, Torsion};
use lumol_core::System;

use crate::{InteractionsInput, FromToml, Error};
use crate::extract;

impl InteractionsInput {
    /// Read the "angles" section from the potential configuration.
    pub(crate) fn read_angles(&self, system: &mut System) -> Result<(), Error> {
        let angles = match self.config.get("angles") {
            Some(angles) => angles,
            None => return Ok(()),
        };

        let angles = angles.as_array().ok_or(Error::from("The 'angles' section must be an array"))?;

        for angle in angles {
            let angle = angle.as_table().ok_or(Error::from("Angle potential entry must be a table"))?;

            let atoms = extract::slice("atoms", angle, "angle potential")?;
            if atoms.len() != 3 {
                return Err(Error::from(format!(
                    "Wrong size for 'atoms' array in angle potential. Should be 3, is {}",
                    atoms.len()
                )));
            }

            let a = atoms[0].as_str().ok_or(
                Error::from("The first atom name is not a string in angle potential")
            )?;
            let b = atoms[1].as_str().ok_or(
                Error::from("The second atom name is not a string in angle potential")
            )?;
            let c = atoms[2].as_str().ok_or(
                Error::from("The third atom name is not a string in angle potential")
            )?;

            let potential = read_angle_potential(angle)?;
            system.set_angle_potential((a, b, c), potential);
        }
        Ok(())
    }

    /// Read the "dihedrals" section from the potential configuration.
    pub(crate) fn read_dihedrals(&self, system: &mut System) -> Result<(), Error> {
        let dihedrals = match self.config.get("dihedrals") {
            Some(dihedrals) => dihedrals,
            None => return Ok(()),
        };

        let dihedrals = dihedrals.as_array().ok_or(
            Error::from("The 'dihedrals' section must be an array")
        )?;

        for dihedral in dihedrals {
            let dihedral = dihedral.as_table().ok_or(
                Error::from("dihedral potential entry must be a table")
            )?;

            let atoms = extract::slice("atoms", dihedral, "dihedral potential")?;
            if atoms.len() != 4 {
                return Err(Error::from(format!(
                    "Wrong size for 'atoms' array in dihedral potential. Should be 4, is {}",
                    atoms.len()
                )));
            }

            let a = atoms[0].as_str().ok_or(
                Error::from("The first atom name is not a string in dihedral potential")
            )?;
            let b = atoms[1].as_str().ok_or(
                Error::from("The second atom name is not a string in dihedral potential")
            )?;
            let c = atoms[2].as_str().ok_or(
                Error::from("The third atom name is not a string in dihedral potential")
            )?;
            let d = atoms[3].as_str().ok_or(
                Error::from("The fourth atom name is not a string in dihedral potential")
            )?;

            let potential = read_dihedral_potential(dihedral)?;
            system.set_dihedral_potential((a, b, c, d), potential);
        }
        Ok(())
    }
}

fn read_angle_potential(angle: &Table) -> Result<Box<dyn AnglePotential>, Error> {
    let potentials = angle.keys().cloned().filter(|key| key != "atoms").collect::<Vec<_>>();

    if potentials.is_empty() {
        return Err(Error::from("Missing potential type in angle potential"));
    }

    if potentials.len() > 1 {
        return Err(Error::from(format!(
            "Got more than one potential type in angle potential: {}",
            potentials.join(" and ")
        )));
    }

    let key = &*potentials[0];
    if let Value::Table(ref table) = angle[key] {
        match key {
            "null" => Ok(Box::new(NullPotential::from_toml(table)?)),
            "harmonic" => Ok(Box::new(Harmonic::from_toml(table)?)),
            "cosine-harmonic" => Ok(Box::new(CosineHarmonic::from_toml(table)?)),
            "morse" => Ok(Box::new(Morse::from_toml(table)?)),
            other => Err(Error::from(format!("Unknown potential type '{}'", other))),
        }
    } else {
        Err(Error::from(format!("'{}' potential must be a table", key)))
    }
}

fn read_dihedral_potential(dihedral: &Table) -> Result<Box<dyn DihedralPotential>, Error> {
    let potentials = dihedral.keys().cloned().filter(|key| key != "atoms").collect::<Vec<_>>();

    if potentials.is_empty() {
        return Err(Error::from("Missing potential type in dihedral potential"));
    }

    if potentials.len() > 1 {
        return Err(Error::from(format!(
            "Got more than one potential type in dihedral potential: {}",
            potentials.join(" and ")
        )));
    }

    let key = &*potentials[0];
    if let Value::Table(ref table) = dihedral[key] {
        match key {
            "null" => Ok(Box::new(NullPotential::from_toml(table)?)),
            "harmonic" => Ok(Box::new(Harmonic::from_toml(table)?)),
            "cosine-harmonic" => Ok(Box::new(CosineHarmonic::from_toml(table)?)),
            "torsion" => Ok(Box::new(Torsion::from_toml(table)?)),
            "morse" => Ok(Box::new(Morse::from_toml(table)?)),
            other => Err(Error::from(format!("Unknown potential type '{}'", other))),
        }
    } else {
        Err(Error::from(format!("'{}' potential must be a table", key)))
    }
}
