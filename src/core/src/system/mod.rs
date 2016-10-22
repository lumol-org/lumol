// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

//! The `system` module provide a way to store data about a simulated system.
//! These systems are represented by an `System` instance, made of a list of
//! `Particle`, an enclosing `UnitCell` and some interactions.

mod periodic;
pub use self::periodic::{PeriodicTable, ElementData};

mod particles;
pub use self::particles::{Particle, ParticleKind};

mod cells;
pub use self::cells::{UnitCell, CellShape};

pub mod velocities;
pub use self::velocities::InitVelocities;
pub use self::velocities::{BoltzmannVelocities, UniformVelocities};

mod molecules;
pub use self::molecules::{Bond, Angle, Dihedral};
pub use self::molecules::Connectivity;
pub use self::molecules::{CONNECT_12, CONNECT_13, CONNECT_14, CONNECT_FAR};
pub use self::molecules::Molecule;
pub use self::molecules::molecule_type;

mod interactions;

mod energy;
pub use self::energy::EnergyEvaluator;

mod cache;
pub use self::cache::EnergyCache;

mod systems;
pub use self::systems::System;
pub use self::systems::Permutations;

pub mod compute;
pub use self::compute::Compute;

mod chfl;
pub use self::chfl::{Trajectory, TrajectoryError};
pub use self::chfl::{guess_bonds, read_molecule};
pub use self::chfl::ToChemfiles;
