// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 Lumol's contributors — BSD license
use toml::de::from_str as parse;
use toml::value::Table;

use std::io::prelude::*;
use std::path::{Path, PathBuf};
use std::fs::File;

use validate;
use error::{Error, Result};

use lumol::sim::Simulation;
use lumol::sys::System;

mod logging;
mod system;
mod outputs;
mod propagator;
mod simulations;
mod min;
mod md;
mod mc;

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
    pub fn new<P: Into<PathBuf>>(path: P) -> Result<Input> {
        let path = path.into();
        let mut file = try_io!(File::open(&path), path);
        let mut buffer = String::new();
        let _ = try_io!(file.read_to_string(&mut buffer), path);
        return Input::from_str(path, &buffer);
    }

    /// Read the `Input` from a TOML formatted string.
    // TODO: use restricted privacy here
    #[doc(hidden)]
    pub fn from_str(path: PathBuf, string: &str) -> Result<Input> {
        let config = try!(parse(string).map_err(|err| {
            Error::TOML(Box::new(err))
        }));
        try!(validate(&config));
        Ok(Input{path: path, config: config.clone()})
    }

    /// Read input file and get the corresponding `Config`
    pub fn read(&self) -> Result<Config> {
        try!(self.setup_logging());
        let system = try!(self.read_system());
        let simulation = try!(self.read_simulation());
        let nsteps = try!(self.read_nsteps());

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
