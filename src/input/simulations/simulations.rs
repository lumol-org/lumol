// Cymbalum, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license
use toml::Table;
use input::error::{Error, Result};
use simulation::Simulation;

use simulation::MolecularDynamics;

pub fn read_simulation(config: &Table) -> Result<Simulation> {
    Ok(Simulation::new(Box::new(MolecularDynamics::new(1.0))))
}

pub fn read_nsteps(config: &Table) -> Result<usize> {
    let simulation = try!(simulation_table(config));
    let nsteps = try!(simulation.get("nsteps").ok_or(
        Error::from("Missing 'nsteps' key in simulation")
    ));

    let nsteps = try!(nsteps.as_integer().ok_or(
        Error::from("'nsteps' key must be an integer")
    ));

    Ok(nsteps as usize)
}

fn simulation_table(config: &Table) -> Result<&Table> {
    let simulations = try!(config.get("simulations").ok_or(
        Error::from("Missing 'simulations' section")
    ));

    let simulations = try!(simulations.as_slice().ok_or(
        Error::from("Missing 'simulations' section is not a table.")
    ));

    if simulations.len() != 1 {
        return Err(Error::from(
            "Only one simulation is supported in the input"
        ));
    }

    let simulation = try!(simulations[0].as_table().ok_or(
        Error::from("Simulations should be tables")
    ));

    return Ok(simulation);
}

#[cfg(test)]
mod tests {
    use input::read_config;
    use input::testing::bad_inputs;
    use std::path::Path;

    #[test]
    fn nsteps() {
        let data_root = Path::new(file!()).parent().unwrap().join("data");
        let config = read_config(data_root.join("simulation.toml")).unwrap();
        assert_eq!(config.nsteps, 1000000);
    }

    #[test]
    fn bad_nsteps() {
        for path in bad_inputs("simulations", "nsteps") {
            assert!(read_config(path).is_err());
        }
    }
}
