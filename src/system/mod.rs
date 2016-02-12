/* Cymbalum, Molecular Simulation in Rust - Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
 */

//! The `system` module provide a way to store data about a simulated system.
//! These systems are represented by an `System` instance, made of a list of
//! `Particle`, an enclosing `UnitCell` and some interactions.

mod periodic;
pub use self::periodic::{PeriodicTable, ElementData};

mod particles;
pub use self::particles::Particle;

mod cells;
pub use self::cells::{UnitCell, CellType};

pub mod velocities;
pub use self::velocities::InitVelocities;
pub use self::velocities::{BoltzmanVelocities, UniformVelocities};

mod molecules;
pub use self::molecules::{Bond, Angle, Dihedral};
pub use self::molecules::Connectivity;
pub use self::molecules::{CONNECT_12, CONNECT_13, CONNECT_14, CONNECT_FAR};
pub use self::molecules::Molecule;
pub use self::molecules::moltype;

mod interactions;
pub use self::interactions::PairInteraction;

mod energy;
pub use self::energy::EnergyEvaluator;

mod cache;
pub use self::cache::EnergyCache;

mod systems;
pub use self::systems::System;

pub mod chemfiles;
