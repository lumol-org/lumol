// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

//! The `system` module provide a way to store data about a simulated system.
//! These systems are represented by an `System` instance, made of a list of
//! `Particle`, an enclosing `UnitCell` and some interactions.

pub use sys2::{PeriodicTable, ElementData};
pub use sys2::{Particle, ParticleKind};
pub use sys2::Composition;
pub use sys2::{UnitCell, CellShape};
pub use sys2::{Bond, Angle, Dihedral};
pub use sys2::Connectivity;
pub use sys2::{CONNECT_12, CONNECT_13, CONNECT_14, CONNECT_FAR};
pub use sys2::Molecule;
pub use sys2::molecule_type;

mod interactions;

mod energy;
pub use self::energy::EnergyEvaluator;

mod cache;
pub use self::cache::EnergyCache;

mod systems;
pub use self::systems::System;
pub use self::systems::Permutations;

mod chfl;
pub use self::chfl::{Trajectory, TrajectoryError};
pub use self::chfl::{guess_bonds, read_molecule};
pub use self::chfl::ToChemfiles;

pub mod veloc;
pub mod compute;
