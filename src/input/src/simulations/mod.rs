// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license
use toml::{Parser, Table};

use std::io::prelude::*;
use std::path::{Path, PathBuf};
use std::fs::File;

use validate;
use error::{Error, Result};
use error::toml_error_to_string;

use lumol::sim::Simulation;
use lumol::sys::System;

mod system;
mod outputs;
mod propagator;
mod simulations;
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
    pub fn from_str(path: PathBuf, buffer: &str) -> Result<Input> {
        let mut parser = Parser::new(buffer);
        let config = try!(parser.parse().ok_or(
            Error::TOML(toml_error_to_string(&parser))
        ));

        try!(validate(&config));
        Ok(Input{path: path, config: config})
    }

    /// Read input file and get the corresponding `Config`
    pub fn read(&self) -> Result<Config> {
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
