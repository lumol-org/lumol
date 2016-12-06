// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license
use toml::Table;
use std::path::PathBuf;

use lumol::out::Output;
use lumol::out::{TrajectoryOutput, CellOutput, EnergyOutput, PropertiesOutput};

use error::{Error, Result};
use FromToml;
use extract;
use super::Input;

impl Input {
    /// Get the the simulation outputs. This is an internal function, public
    /// because of the code organization.
    // TODO: use restricted privacy here
    #[doc(hidden)]
    pub fn read_outputs(&self) -> Result<Vec<(Box<Output>, u64)>> {
        let config = try!(self.simulation_table());
        if let Some(outputs) = config.get("outputs") {
            let outputs = try!(outputs.as_slice().ok_or(
                Error::from("'outputs' must be an array of tables in simulation")
            ));

            let mut result = Vec::new();
            for output in outputs {
                let output = try!(output.as_table().ok_or(
                    Error::from("'outputs' must be an array of tables in simulation")
                ));

                let frequency = match output.get("frequency") {
                    Some(frequency) => {
                        try!(frequency.as_integer().ok_or(
                            Error::from("'frequency' must be an integer in output")
                        )) as u64
                    },
                    None => 1u64
                };

                let output: Box<Output> = match try!(extract::typ(output, "output")) {
                    "Trajectory" | "trajectory" => Box::new(try!(TrajectoryOutput::from_toml(output))),
                    "Energy" | "energy" => Box::new(try!(EnergyOutput::from_toml(output))),
                    "Cell" | "cell" => Box::new(try!(CellOutput::from_toml(output))),
                    "Properties" | "properties" => Box::new(try!(PropertiesOutput::from_toml(output))),
                    other => {
                        return Err(Error::from(
                            format!("Unknown output type '{}'", other)
                        ))
                    }
                };

                result.push((output, frequency));
            }
            Ok(result)
        } else {
            Ok(Vec::new())
        }
    }
}

fn get_file(config: &Table) -> Result<&str> {
    let file = try!(config.get("file").ok_or(
        Error::from("Missing 'file' key in output")
    ));

    file.as_str().ok_or(
        Error::from("'file' must be a string in output")
    )
}

impl FromToml for TrajectoryOutput {
    fn from_toml(config: &Table) -> Result<TrajectoryOutput> {
        let path = try!(get_file(config));
        let output = try!(TrajectoryOutput::new(path));
        Ok(output)
    }
}

impl FromToml for CellOutput {
    fn from_toml(config: &Table) -> Result<CellOutput> {
        let path = try!(get_file(config));
        let output = try_io!(CellOutput::new(path), PathBuf::from(path));
        Ok(output)
    }
}

impl FromToml for EnergyOutput {
    fn from_toml(config: &Table) -> Result<EnergyOutput> {
        let path = try!(get_file(config));
        let output = try_io!(EnergyOutput::new(path), PathBuf::from(path));
        Ok(output)
    }
}

impl FromToml for PropertiesOutput {
    fn from_toml(config: &Table) -> Result<PropertiesOutput> {
        let path = try!(get_file(config));
        let output = try_io!(PropertiesOutput::new(path), PathBuf::from(path));
        Ok(output)
    }
}
