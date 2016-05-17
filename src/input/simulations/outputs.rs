// Cymbalum, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license
use toml::{Value, Table};
use input::error::{Error, Result};
use input::FromToml;

use simulation::Output;
use simulation::{TrajectoryOutput, CellOutput, EnergyOutput, PropertiesOutput};

pub fn read_outputs(outputs: &[Value]) -> Result<Vec<(Box<Output>, u64)>> {
    let mut res = Vec::new();
    for output in outputs {
        let output = try!(output.as_table().ok_or(
            Error::from("All values in outputs should be tables")
        ));

        let frequency = match output.get("frequency") {
            Some(frequency) => {
                try!(frequency.as_integer().ok_or(
                    Error::from("'frequency' must be an integer")
                )) as u64
            },
            None => 1u64
        };

        let output: Box<Output> = match extract_type!(output) {
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

        res.push((output, frequency));
    }
    Ok(res)
}

fn get_file(config: &Table) -> Result<&str> {
    let file = try!(config.get("file").ok_or(
        Error::from("Missing 'file' in output")
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
        let output = try!(CellOutput::new(path));
        Ok(output)
    }
}

impl FromToml for EnergyOutput {
    fn from_toml(config: &Table) -> Result<EnergyOutput> {
        let path = try!(get_file(config));
        let output = try!(EnergyOutput::new(path));
        Ok(output)
    }
}

impl FromToml for PropertiesOutput {
    fn from_toml(config: &Table) -> Result<PropertiesOutput> {
        let path = try!(get_file(config));
        let output = try!(PropertiesOutput::new(path));
        Ok(output)
    }
}

#[cfg(test)]
mod tests {
    use input::read_config;
    use input::testing::{bad_inputs, cleanup};
    use std::path::Path;

    #[test]
    fn outputs() {
        let path = Path::new(file!()).parent().unwrap()
                                     .join("data")
                                     .join("md.toml");
        let config = read_config(&path).unwrap();
        assert_eq!(config.simulation.outputs_len(), 2);
        cleanup(&path);
    }

    #[test]
    fn bad_outputs() {
        for path in bad_inputs("simulations", "outputs") {
            assert!(read_config(path).is_err());
        }
    }
}
