// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license
use lumol::sim::Simulation;
use toml::value::Table;

use super::Input;
use error::{Error, Result};
use extract;

impl Input {
    /// Get the the simulation.
    pub fn read_simulation(&self) -> Result<Simulation> {
        let propagator = self.read_propagator()?;
        let mut simulation = Simulation::new(propagator);
        for (output, frequency) in self.read_outputs()? {
            simulation.add_output_with_frequency(output, frequency);
        }

        Ok(simulation)
    }

    /// Get the number of steps in the simulation.
    pub(crate) fn read_nsteps(&self) -> Result<usize> {
        let simulation = self.simulation_table()?;
        let nsteps = simulation.get("nsteps").ok_or(
            Error::from("Missing 'nsteps' key in simulation")
        )?;

        let nsteps = nsteps.as_integer().ok_or(
            Error::from("'nsteps' key must be an integer")
        )?;

        Ok(nsteps as usize)
    }

    /// Get the simulation TOML table.
    pub(crate) fn simulation_table(&self) -> Result<&Table> {
        let simulations = extract::slice("simulations", &self.config, "input file")?;
        if simulations.len() != 1 {
            return Err(Error::from("Only one simulation is supported in the input"));
        }

        let simulation = simulations[0].as_table().ok_or(
            Error::from("Simulations should be tables")
        )?;

        return Ok(simulation);
    }
}
