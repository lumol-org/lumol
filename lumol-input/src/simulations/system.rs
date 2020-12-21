// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license
use toml::value::{Table, Value};

use lumol_core::{System, UnitCell, TrajectoryBuilder};
use lumol_sim::{BoltzmannVelocities, InitVelocities};
use lumol_core::units;

use log::warn;

use crate::{Input, InteractionsInput, Error};
use crate::extract;
use crate::simulations::get_input_path;

impl Input {
    /// Get the the simulated system.
    pub fn read_system(&self) -> Result<System, Error> {
        let config = self.system_table()?;

        let file = extract::str("file", config, "system")?;
        let file = get_input_path(&self.path, file);
        let mut trajectory = TrajectoryBuilder::new().open(file)?;

        let with_cell = self.read_cell()?.map_or(false, |cell| {
            trajectory.set_cell(&cell);
            true
        });

        if config.get("topology").is_some() {
            let topology = extract::str("topology", config, "system")?;
            trajectory.set_topology_file(topology)?;
        }

        let guess_bonds = if let Some(guess_bonds) = config.get("guess_bonds") {
            guess_bonds.as_bool().ok_or(
                Error::from("'guess_bonds' should be a boolean value in system")
            )?
        } else {
            false
        };

        let mut system = if guess_bonds {
            trajectory.read_guess_bonds()?
        } else {
            trajectory.read()?
        };

        self.read_potentials(&mut system)?;
        self.init_velocities(&mut system)?;

        if !with_cell && system.cell.is_infinite() {
            warn!(
                "No unit cell in the system, using an infinite unit cell.\n\
                 You can get rid of this warning by using `cell = []` in the \
                 input file if this is what you want."
            );
        }

        Ok(system)
    }

    fn system_table(&self) -> Result<&Table, Error> {
        let systems = extract::slice("systems", &self.config, "input file")?;

        if systems.is_empty() {
            return Err(Error::from("'systems' array should contain a system"));
        }

        if systems.len() > 1 {
            return Err(Error::from("only one system is supported in input file"));
        }

        let system = systems[0].as_table().ok_or(
            Error::from("'systems' should be an array of tables in input file")
        )?;

        return Ok(system);
    }

    fn read_cell(&self) -> Result<Option<UnitCell>, Error> {
        let config = self.system_table()?;
        if let Some(cell) = config.get("cell") {
            match *cell {
                Value::Array(ref cell) => {
                    if cell.is_empty() {
                        Ok(Some(UnitCell::infinite()))
                    } else if cell.len() == 3 {
                        let a = get_cell_number(&cell[0])?;
                        let b = get_cell_number(&cell[1])?;
                        let c = get_cell_number(&cell[2])?;

                        Ok(Some(UnitCell::ortho(a, b, c)))
                    } else if cell.len() == 6 {
                        let a = get_cell_number(&cell[0])?;
                        let b = get_cell_number(&cell[1])?;
                        let c = get_cell_number(&cell[2])?;
                        let alpha = get_cell_number(&cell[3])?;
                        let beta = get_cell_number(&cell[4])?;
                        let gamma = get_cell_number(&cell[5])?;

                        Ok(Some(UnitCell::triclinic(a, b, c, alpha, beta, gamma)))
                    } else {
                        Err(Error::from("'cell' array must have a size of 3 or 6"))
                    }
                }
                Value::Integer(lenght) => {
                    let lenght = lenght as f64;
                    Ok(Some(UnitCell::cubic(lenght)))
                }
                Value::Float(lenght) => Ok(Some(UnitCell::cubic(lenght))),
                _ => Err(Error::from("'cell' must be a number or an array in system")),
            }
        } else {
            Ok(None)
        }
    }

    fn init_velocities(&self, system: &mut System) -> Result<(), Error> {
        let config = self.system_table()?;

        if let Some(velocities) = config.get("velocities") {
            let velocities = velocities.as_table().ok_or(
                Error::from("'velocities' must be a table in system")
            )?;

            if velocities.get("init").is_some() {
                let temperature = extract::str("init", velocities, "velocities initializer")?;
                let temperature = units::from_str(temperature)?;
                let mut velocities = BoltzmannVelocities::new(temperature);
                velocities.init(system);
            } else {
                warn!("'velocities' key does nothing in this input file");
            }
        }

        Ok(())
    }

    fn read_potentials(&self, system: &mut System) -> Result<(), Error> {
        let config = self.system_table()?;
        if let Some(potentials) = config.get("potentials") {
            if let Some(potentials) = potentials.as_str() {
                let path = get_input_path(&self.path, potentials);
                let input = InteractionsInput::new(path)?;
                input.read(system)?;
            } else if let Some(potentials) = potentials.as_table() {
                let input = InteractionsInput::from_toml(potentials.clone());
                input.read(system)?;
            } else {
                return Err(Error::from("'potentials' must be a string or a table in system"));
            }
        } else {
            warn!("No potentials found in input file");
        }
        Ok(())
    }
}

#[allow(clippy::option_if_let_else)]
fn get_cell_number(value: &Value) -> Result<f64, Error> {
    if let Some(value) = value.as_integer() {
        Ok(value as f64)
    } else if let Some(value) = value.as_float() {
        Ok(value)
    } else {
        Err(Error::from("values must be numbers in 'cell' array"))
    }
}
