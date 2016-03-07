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

mod adaptators;

/******************************************************************************/
use std::path::Path;
use system::{System, Molecule, Particle};
use chemfiles;

use self::adaptators::{ToChemfiles, ToCymbalum};

/// A Trajectory is a file containing one or more successives simulation steps
pub struct Trajectory(chemfiles::Trajectory);
pub use chemfiles::Error as TrajectoryError;
/// Result type for all Trajectory operations
pub type TrajectoryResult<T> = Result<T, TrajectoryError>;

impl Trajectory {
    /// Open an existing file at `path` for reading.
    pub fn open<P: AsRef<Path>>(path: P) -> TrajectoryResult<Trajectory> {
        let traj = try!(chemfiles::Trajectory::open(path));
        return Ok(Trajectory(traj));
    }

    /// Create a new file at `path` for writing, and overwrite any existing file.
    pub fn create<P: AsRef<Path>>(path: P) -> TrajectoryResult<Trajectory> {
        let traj = try!(chemfiles::Trajectory::create(path));
        return Ok(Trajectory(traj));
    }

    /// Read the next step of the trajectory
    pub fn read(&mut self) -> TrajectoryResult<System> {
        let mut frame = try!(chemfiles::Frame::new(0));
        try!(self.0.read(&mut frame));
        return frame.to_cymbalum();
    }

    /// Read the next step of the trajectory, and guess the bonds of the
    /// resulting System.
    pub fn read_guess_bonds(&mut self) -> TrajectoryResult<System> {
        let mut frame = try!(chemfiles::Frame::new(0));
        try!(self.0.read(&mut frame));
        try!(frame.guess_topology(true));
        return frame.to_cymbalum();
    }

    /// Write the system to the trajectory.
    pub fn write(&mut self, system: &System) -> TrajectoryResult<()> {
        let frame = try!(system.to_chemfiles());
        return self.0.write(&frame);
    }
}

/// Read a the first molecule from the file at `path`. If no bond information
/// exists in the file, bonds are guessed.
pub fn read_molecule<P: AsRef<Path>>(path: P) -> TrajectoryResult<(Molecule, Vec<Particle>)> {
    let path = path.as_ref();

    let mut trajectory = try!(chemfiles::Trajectory::open(path));
    let mut frame = try!(chemfiles::Frame::new(0));
    try!(trajectory.read(&mut frame));

    // Only guess the topology when we have no bond information
    let topology = try!(frame.topology());
    if try!(topology.bonds_count()) == 0 {
        try!(frame.guess_topology(true));
    }
    let system = try!(frame.to_cymbalum());

    assert!(system.size() != 0, "No molecule in the file at {:?}", path);
    let molecule = system.molecule(0).clone();
    let mut particles = Vec::new();
    for i in &molecule {
        particles.push(system[i].clone());
    }
    return Ok((molecule, particles));
}


#[cfg(test)]
pub mod testing {
    //! Some utilities for internal unit tests
    use system::System;
    use super::TrajectoryResult;
    use super::adaptators::{ToChemfiles, ToCymbalum};

    pub use super::interactions::read_interactions_string;

    /// Guess bonds in a system
    pub fn guess_bonds(system: System) -> TrajectoryResult<System> {
        let mut frame = try!(system.to_chemfiles());
        try!(frame.guess_topology(true));
        return frame.to_cymbalum();
    }
}
