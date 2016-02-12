/* Cymbalum, Molecular Simulation in Rust - Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
 */
//! A propagator is responsible for updating the system during a simulation
use system::System;

/// The propagator trait is the main algorithm of a simulation, i.e. the one
/// which update the system. The main function here is `propagate`, which
/// should propagate the simulation for one step.
pub trait Propagator {
    /// Setup code, preparing all the meta-informations needed about the
    /// simulation
    fn setup(&mut self, _: &System) {}

    /// Propagate the system for one simulation step.
    fn propagate(&mut self, system: &mut System);

    /// Finish the simulation, and maybe output some informations about it
    fn finish(&mut self, _: &System) {}
}
