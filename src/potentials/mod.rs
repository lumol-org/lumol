/*
 * Cymbalum, Molecular Simulation in Rust
 * Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
*/

//! Energy and forces computation in simulations.
//!
//! In order to compute an energy or a force, two things are needed: a function
//! mapping distance/angles to energy/force, and a way to compute theses
//! functions. The functions implements the `PotentialFunction` trait, and the
//! `PotentialComputation` trait provides ways to effectively compute the
//! energy: with a cutoff, from a table, *etc.*

pub mod functions;
pub use self::functions::{PotentialFunction, PairPotential};
pub use self::functions::{LennardJones, Harmonic};

pub mod computations;
pub use self::computations::PotentialComputation;
pub use self::computations::{DirectComputation, TableComputation, CutoffComputation};
