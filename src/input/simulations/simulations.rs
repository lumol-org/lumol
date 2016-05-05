// Cymbalum, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license
use toml::Table;
use input::error::{Error, Result};
use simulation::Simulation;

use super::outputs::read_outputs;
use super::propagator::read_propagator;

pub fn read_simulation(config: &Table) -> Result<Simulation> {
    let config = try!(simulation_table(config));
    let propagator = try!(read_propagator(config));

    let mut simulation = Simulation::new(propagator);

    if let Some(outputs) = config.get("outputs") {
        let outputs = try!(outputs.as_slice().ok_or(
            Error::from("'outputs' must be an array")
        ));

        let outputs = try!(read_outputs(outputs));
        for (output, frequency) in outputs {
            simulation.add_output_with_frequency(output, frequency);
        }
    }

    Ok(simulation)
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
    let simulations = extract_slice!("simulations", config as "input file");
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
    use input::testing::{bad_inputs, cleanup};
    use std::path::Path;

    #[test]
    fn nsteps() {
        let path = Path::new(file!()).parent().unwrap()
                                     .join("data")
                                     .join("md.toml");
        let config = read_config(&path).unwrap();
        assert_eq!(config.nsteps, 1000000);
        cleanup(&path);
    }

    #[test]
    fn bad_nsteps() {
        for path in bad_inputs("simulations", "nsteps") {
            assert!(read_config(path).is_err());
        }
    }
}
