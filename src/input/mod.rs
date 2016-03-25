// Cymbalum, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

//! Input and output capacities.
//!
//! # Input files
//!
//! This module provide readers for two types of file:
//! * Configuration files uses the YAML format to define the `System` or the
//!   `Simulation` in an human-readable way, and without writing any code;
//! * Structures files defines the positions and the names of particles in the
//!   `System` or in a specific `Molecule`.
//!
//! # Reading YAML files
//!
//! The `read_interactions` function read interactions into a `System`.
//!
//! # Reading and writing trajectories
//!
//! The `Trajectory` type allow to read and write trajectories on the disk.
//! One can also read a single molecule from a file using the `read_molecule`
//! function.

mod interactions;
pub use self::interactions::read_interactions;
pub use self::interactions::Error as InteractionError;

mod chemfiles;
pub use self::chemfiles::{Trajectory, TrajectoryResult};
pub use self::chemfiles::{guess_bonds, read_molecule};

#[cfg(test)]
pub use self::interactions::read_interactions_string;
