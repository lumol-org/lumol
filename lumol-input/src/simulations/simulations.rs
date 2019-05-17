// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license
use lumol_sim::Simulation;
use toml::value::Table;

use crate::{Input, Error};
use crate::extract;

impl Input {
    /// Get the the simulation.
    pub fn read_simulation(&self) -> Result<Simulation, Error> {
        let propagator = self.read_propagator()?;
        let mut simulation = Simulation::new(propagator);
        for (output, frequency) in self.read_outputs()? {
            simulation.add_output_with_frequency(output, frequency);
        }

        Ok(simulation)
    }

    /// Get the number of steps in the simulation.
    pub(crate) fn read_nsteps(&self) -> Result<usize, Error> {
        let simulation = self.simulation_table()?;
        let nsteps = simulation.get("nsteps").ok_or(
            Error::from("missing 'nsteps' key in simulation")
        )?;

        let nsteps = nsteps.as_integer().ok_or(
            Error::from("'nsteps' key must be an integer")
        )?;

        Ok(nsteps as usize)
    }

    /// Get the simulation TOML table.
    pub(crate) fn simulation_table(&self) -> Result<&Table, Error> {
        let simulations = extract::slice("simulations", &self.config, "input file")?;
        if simulations.len() != 1 {
            return Err(Error::from("only one simulation is supported in the input"));
        }

        let simulation = simulations[0].as_table().ok_or(
            Error::from("simulations should be tables")
        )?;

        return Ok(simulation);
    }
}
