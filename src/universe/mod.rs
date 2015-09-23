/*
 * Cymbalum, Molecular Simulation in Rust
 * Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
*/

//! The `universe` module provide a way to store data about a simulated system.
//! These systems are represented by an `Universe` instance, made of a list of
//! `Particle`, an enclosing `UnitCell` and some interactions.

mod periodic;
pub use self::periodic::{PeriodicTable, ElementData};
pub use self::periodic::PERIODIC_TABLE;

mod particles;
pub use self::particles::Particle;

mod cells;
pub use self::cells::UnitCell;

pub mod velocities;
pub use self::velocities::InitVelocities;
pub use self::velocities::{BoltzmanVelocities, UniformVelocities};

mod topologies;
pub use self::topologies::{Bond, Angle, Dihedral};
pub use self::topologies::Topology;

mod interactions;
mod universes;
pub use self::universes::Universe;

pub mod chemharp;
