// Cymbalum, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

//! Input capacities for cymbalum
//!
//! This module provide input files reader for two types of files:
//!
//! * Configuration files uses the TOML format to define a `System` or a
//!   `Simulation` in an human-readable way, without writing any code. This
//!   type of input is further divided in potentials input files and whole
//!   simulation input files;
//! * Structures files defines the positions and the names of particles in a
//!   `System` or in a specific `Molecule`.
//!
//! # Reading configuration files
//!
//! The `read_interactions` function read interactions into a `System`, and the
//! `read_config` function reads a whole simulation (`Simulation` and `System`
//! objects).
//!
//! # Reading and writing trajectories
//!
//! The `Trajectory` type allow to read and write trajectories on the disk.
//! One can also read a single molecule from a file using the `read_molecule`
//! function.

mod error;
pub use self::error::{Error, TrajectoryError};

#[macro_use]
mod macros;

mod interactions;
pub use self::interactions::read_interactions;
pub use self::interactions::FromTomlWithPairs;

mod simulations;
pub use self::simulations::read_config;

mod chemfiles;
pub use self::chemfiles::Trajectory;
pub use self::chemfiles::{guess_bonds, read_molecule};

#[cfg(test)]
pub mod testing;

#[cfg(test)]
pub use self::interactions::read_interactions_string;

use toml::Table;

/// Convert a TOML table to a Rust type.
pub trait FromToml: Sized {
    /// Do the conversion from `table` to Self.
    fn from_toml(table: &Table) -> Result<Self, Error>;
}

/// Convert a TOML table and some additional data to a Rust type.
pub trait FromTomlWithData: Sized {
    /// The type of the additional data needed.
    type Data;
    /// Do the conversion from `table` and `data` to Self.
    fn from_toml(table: &Table, data: Self::Data) -> Result<Self, Error>;
}
