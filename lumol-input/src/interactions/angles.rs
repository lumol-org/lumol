// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license
use toml::value::Table;

use lumol_core::energy::{AnglePotential, DihedralPotential};
use lumol_core::energy::{CosineHarmonic, Harmonic, Morse, NullPotential, Torsion};
use lumol_core::System;

use crate::{InteractionsInput, FromToml, Error};
use crate::extract;

impl InteractionsInput {
    /// Read the "angles" section from the potential configuration.
    pub(crate) fn read_angles(&self, system: &mut System) -> Result<(), Error> {
        let Some(angles) = self.config.get("angles") else { return Ok(()) };

        let angles = angles.as_table().ok_or(
            Error::from("the 'angles' section must be a table")
        )?;

        for (key, table) in angles {
            let atoms = key.split('-').collect::<Vec<_>>();
            if atoms.len() != 3 {
                return Err(Error::from(format!(
                    "expected three atoms for angle potential, got {} ({:?})", atoms.len(), atoms
                )));
            }

            let table = table.as_table().ok_or(
                Error::from(format!(
                    "angle potential associated with {key} must be a table"
                ))
            )?;

            let potential = read_angle_potential(table)?;
            system.set_angle_potential((atoms[0], atoms[1], atoms[2]), potential);
        }
        Ok(())
    }

    /// Read the "dihedrals" section from the potential configuration.
    pub(crate) fn read_dihedrals(&self, system: &mut System) -> Result<(), Error> {
        let Some(dihedrals) = self.config.get("dihedrals") else { return Ok(()) };

        let dihedrals = dihedrals.as_table().ok_or(
            Error::from("the 'dihedrals' section must be a table")
        )?;

        for (key, table) in dihedrals {
            let atoms = key.split('-').collect::<Vec<_>>();
            if atoms.len() != 4 {
                return Err(Error::from(format!(
                    "expected four atoms for dihedral potential, got {} ({:?})", atoms.len(), atoms
                )));
            }

            let table = table.as_table().ok_or(
                Error::from(format!(
                    "dihedral potential associated with {key} must be a table"
                ))
            )?;

            let potential = read_dihedral_potential(table)?;
            system.set_dihedral_potential((atoms[0], atoms[1], atoms[2], atoms[3]), potential);
        }
        Ok(())
    }
}

fn read_angle_potential(table: &Table) -> Result<Box<dyn AnglePotential>, Error> {
    match extract::typ(table, "angle potential")? {
        "null" => Ok(Box::new(NullPotential::from_toml(table)?)),
        "harmonic" => Ok(Box::new(Harmonic::from_toml(table)?)),
        "cosine-harmonic" => Ok(Box::new(CosineHarmonic::from_toml(table)?)),
        "morse" => Ok(Box::new(Morse::from_toml(table)?)),
        other => Err(Error::from(format!("unknown potential type '{other}'"))),
    }
}

fn read_dihedral_potential(table: &Table) -> Result<Box<dyn DihedralPotential>, Error> {
    match extract::typ(table, "dihedral potential")? {
        "null" => Ok(Box::new(NullPotential::from_toml(table)?)),
        "harmonic" => Ok(Box::new(Harmonic::from_toml(table)?)),
        "cosine-harmonic" => Ok(Box::new(CosineHarmonic::from_toml(table)?)),
        "torsion" => Ok(Box::new(Torsion::from_toml(table)?)),
        "morse" => Ok(Box::new(Morse::from_toml(table)?)),
        other => Err(Error::from(format!("unknown potential type '{other}'"))),
    }
}
