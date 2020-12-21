// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license
use toml::de::from_str as parse_toml;
use toml::value::Table;

use std::fs::File;
use std::io::prelude::*;
use std::path::{Path, PathBuf};

use lumol_sim::Simulation;
use lumol_core::System;

use crate::Error;
use crate::validate;

mod logging;
mod system;
mod outputs;
mod propagator;
#[allow(clippy::module_inception)]
mod simulations;
mod min;
mod md;
mod mc;

pub use self::logging::setup_default_logger;

/// A configuration about how to run a single simulation. This contains the
/// system to simulate, the simulation itself and the number of steps to run
/// the simulation.
pub struct Config {
    /// The simulated system
    pub system: System,
    /// The simulation object
    pub simulation: Simulation,
    /// The simulation duration
    pub nsteps: usize,
}

/// An input file for Lumol.
pub struct Input {
    /// The input file path
    path: PathBuf,
    /// The TOML configuration
    config: Table,
}

impl Input {
    /// Read the file at `Path` and create a new `Input` from it.
    pub fn new<P: Into<PathBuf>>(path: P) -> Result<Input, Error> {
        let path = path.into();
        let mut file = try_io!(File::open(&path), path);
        let mut buffer = String::new();
        let _ = try_io!(file.read_to_string(&mut buffer), path);
        return Input::from_str(path, &buffer);
    }

    /// Read the `Input` from a TOML formatted string.
    pub fn from_str(path: PathBuf, string: &str) -> Result<Input, Error> {
        let config = parse_toml(string).map_err(|err| { Error::TOML(Box::new(err)) })?;
        validate(&config)?;
        Ok(Input {
            path: path,
            config: config,
        })
    }

    /// Read input file and get the corresponding `Config`
    pub fn read(&self) -> Result<Config, Error> {
        self.setup_logging()?;
        let system = self.read_system()?;
        let simulation = self.read_simulation()?;
        let nsteps = self.read_nsteps()?;

        Ok(Config {
            system: system,
            simulation: simulation,
            nsteps: nsteps,
        })
    }
}

fn get_input_path<P1: AsRef<Path>, P2: AsRef<Path>>(root: P1, path: P2) -> PathBuf {
    let path = PathBuf::from(path.as_ref());
    if path.is_absolute() {
        path
    } else {
        let parent = root.as_ref().parent().expect("Could not get parent path");
        parent.join(path)
    }
}
