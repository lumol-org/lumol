// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license
use toml::Value;

use lumol_core::energy::{CoulombicPotential, Ewald, SharedEwald, Wolf};
use lumol_core::System;

use log::{info, warn};

use super::read_restriction;
use crate::{InteractionsInput, Error, FromToml, FromTomlWithRefData};

impl InteractionsInput {
    /// Read the "coulomb" section from the potential configuration.
    pub(crate) fn read_coulomb(&self, system: &mut System) -> Result<(), Error> {
        let Some(coulomb) = self.config.get("coulomb") else { return Ok(()) };

        let coulomb = coulomb.as_table().ok_or(Error::from("the 'coulomb' section must be a table"))?;

        let solvers = coulomb.keys().filter(|&key| key != "restriction").cloned().collect::<Vec<_>>();

        if solvers.len() != 1 {
            return Err(Error::from(
                format!("got more than one coulombic solver: {}", solvers.join(" and ")),
            ));
        }

        let key = &*solvers[0];
        if let Value::Table(ref table) = coulomb[key] {
            let mut potential: Box<dyn CoulombicPotential> = match key {
                "wolf" => Box::new(Wolf::from_toml(table)?),
                "ewald" => {
                    let ewald = Ewald::from_toml(table, system)?;
                    Box::new(SharedEwald::new(ewald))
                }
                other => return Err(Error::from(format!("unknown coulomb solver '{other}'"))),
            };

            if let Some(restriction) = read_restriction(coulomb)? {
                potential.set_restriction(restriction);
            }

            system.set_coulomb_potential(potential);
            Ok(())
        } else {
            Err(Error::from(format!("coulombic solver '{key}' must be a table")))
        }
    }

    /// Read the "charges" from the potential configuration.
    pub(crate) fn read_charges(&self, system: &mut System) -> Result<(), Error> {
        let Some(charges) = self.config.get("charges") else { return Ok(()) };

        let charges = charges.as_table().ok_or(
            Error::from("the 'charges' section must be a table")
        )?;

        let mut total_charge = 0.0;
        for (name, charge) in charges {
            let charge = match *charge {
                Value::Integer(val) => val as f64,
                Value::Float(val) => val,
                _ => {
                    return Err(Error::from("charges must be numbers"));
                }
            };

            let mut n_changed = 0;
            for particle in system.particles_mut() {
                if particle.name == name {
                    *particle.charge = charge;
                    n_changed += 1;
                    total_charge += charge;
                }
            }

            if n_changed == 0 {
                warn!("No particle with name '{}' was found while setting the charges", name);
            } else {
                info!("Charge set to {:+} for {} {} particles", charge, n_changed, name);
            }
        }

        if total_charge.abs() > 1e-6 {
            warn!("System is not neutral and have a net charge of {:+}", total_charge);
        }
        Ok(())
    }
}
