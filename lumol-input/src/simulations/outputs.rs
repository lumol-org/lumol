// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license
use std::path::PathBuf;
use toml::value::Table;

use lumol_sim::output::Output;
use lumol_sim::output::{TrajectoryOutput, PropertiesOutput, EnergyOutput};
use lumol_sim::output::{ForcesOutput, CellOutput, CustomOutput, StressOutput};

use super::Input;
use crate::FromToml;
use crate::error::{Error, Result};
use crate::extract;

impl Input {
    /// Get the the simulation outputs.
    pub(crate) fn read_outputs(&self) -> Result<Vec<(Box<dyn Output>, u64)>> {
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
                    other => return Err(Error::from(format!("Unknown output type '{}'", other))),
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
    let file = config.get("file").ok_or(
        Error::from("Missing 'file' key in output")
    )?;

    file.as_str().ok_or(Error::from("'file' must be a string in output"))
}

impl FromToml for TrajectoryOutput {
    fn from_toml(config: &Table) -> Result<TrajectoryOutput> {
        let path = get_file(config)?;
        let output = TrajectoryOutput::new(path)?;
        Ok(output)
    }
}

impl FromToml for CellOutput {
    fn from_toml(config: &Table) -> Result<CellOutput> {
        let path = get_file(config)?;
        let output = try_io!(CellOutput::new(path), PathBuf::from(path));
        Ok(output)
    }
}

impl FromToml for EnergyOutput {
    fn from_toml(config: &Table) -> Result<EnergyOutput> {
        let path = get_file(config)?;
        let output = try_io!(EnergyOutput::new(path), PathBuf::from(path));
        Ok(output)
    }
}

impl FromToml for PropertiesOutput {
    fn from_toml(config: &Table) -> Result<PropertiesOutput> {
        let path = get_file(config)?;
        let output = try_io!(PropertiesOutput::new(path), PathBuf::from(path));
        Ok(output)
    }
}

impl FromToml for StressOutput {
    fn from_toml(config: &Table) -> Result<StressOutput> {
        let path = get_file(config)?;
        let output = try_io!(StressOutput::new(path), PathBuf::from(path));
        Ok(output)
    }
}

impl FromToml for ForcesOutput {
    fn from_toml(config: &Table) -> Result<ForcesOutput> {
        let path = get_file(config)?;
        let output = try_io!(ForcesOutput::new(path), PathBuf::from(path));
        Ok(output)
    }
}

impl FromToml for CustomOutput {
    fn from_toml(config: &Table) -> Result<CustomOutput> {
        let path = get_file(config)?;
        let template = extract::str("template", config, "custom output")?;
        let output = try_io!(CustomOutput::new(path, template), PathBuf::from(path));
        Ok(output)
    }
}
