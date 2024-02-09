// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license
use std::path::PathBuf;
use toml::value::Table;

use lumol_sim::output::Output;
use lumol_sim::output::{TrajectoryOutput, PropertiesOutput, EnergyOutput};
use lumol_sim::output::{ForcesOutput, CellOutput, CustomOutput, StressOutput};

use crate::{Input, FromToml, Error};
use crate::extract;

pub type OutputFrequency = (Box<dyn Output>, u64);

impl Input {
    /// Get the the simulation outputs.
    pub(crate) fn read_outputs(&self) -> Result<Vec<OutputFrequency>, Error> {
        let config = self.simulation_table()?;
        if let Some(outputs) = config.get("outputs") {
            let outputs = outputs.as_array().ok_or(
                Error::from("'outputs' must be an array of tables in simulation")
            )?;

            let mut result = Vec::new();
            for output in outputs {
                let output = output.as_table().ok_or(
                    Error::from("'outputs' must be an array of tables in simulation")
                )?;

                let frequency = match output.get("frequency") {
                    Some(frequency) => {
                        frequency.as_integer().ok_or(
                            Error::from("'frequency' must be an integer in output")
                        )? as u64
                    }
                    None => 1,
                };

                let typ = extract::typ(output, "output")?;
                let output: Box<dyn Output> = match &*typ.to_lowercase() {
                    "trajectory" => Box::new(TrajectoryOutput::from_toml(output)?),
                    "properties" => Box::new(PropertiesOutput::from_toml(output)?),
                    "energy" => Box::new(EnergyOutput::from_toml(output)?),
                    "stress" => Box::new(StressOutput::from_toml(output)?),
                    "forces" => Box::new(ForcesOutput::from_toml(output)?),
                    "cell" => Box::new(CellOutput::from_toml(output)?),
                    "custom" => Box::new(CustomOutput::from_toml(output)?),
                    other => return Err(Error::from(format!("unknown output type '{other}'"))),
                };

                result.push((output, frequency));
            }
            Ok(result)
        } else {
            Ok(Vec::new())
        }
    }
}

fn get_file(config: &Table) -> Result<&str, Error> {
    let file = config.get("file").ok_or(
        Error::from("missing 'file' key in output")
    )?;

    file.as_str().ok_or(Error::from("'file' must be a string in output"))
}

impl FromToml for TrajectoryOutput {
    fn from_toml(config: &Table) -> Result<TrajectoryOutput, Error> {
        let path = get_file(config)?;
        let output = TrajectoryOutput::new(path)?;
        Ok(output)
    }
}

impl FromToml for CellOutput {
    fn from_toml(config: &Table) -> Result<CellOutput, Error> {
        let path = get_file(config)?;
        let output = try_io!(CellOutput::new(path), PathBuf::from(path));
        Ok(output)
    }
}

impl FromToml for EnergyOutput {
    fn from_toml(config: &Table) -> Result<EnergyOutput, Error> {
        let path = get_file(config)?;
        let output = try_io!(EnergyOutput::new(path), PathBuf::from(path));
        Ok(output)
    }
}

impl FromToml for PropertiesOutput {
    fn from_toml(config: &Table) -> Result<PropertiesOutput, Error> {
        let path = get_file(config)?;
        let output = try_io!(PropertiesOutput::new(path), PathBuf::from(path));
        Ok(output)
    }
}

impl FromToml for StressOutput {
    fn from_toml(config: &Table) -> Result<StressOutput, Error> {
        let path = get_file(config)?;
        let output = try_io!(StressOutput::new(path), PathBuf::from(path));
        Ok(output)
    }
}

impl FromToml for ForcesOutput {
    fn from_toml(config: &Table) -> Result<ForcesOutput, Error> {
        let path = get_file(config)?;
        let output = try_io!(ForcesOutput::new(path), PathBuf::from(path));
        Ok(output)
    }
}

impl FromToml for CustomOutput {
    fn from_toml(config: &Table) -> Result<CustomOutput, Error> {
        let path = get_file(config)?;
        let template = extract::str("template", config, "custom output")?;
        let output = try_io!(CustomOutput::new(path, template), PathBuf::from(path));
        Ok(output)
    }
}
