/*
 * Cymbalum, Molecular Simulation in Rust
 * Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
*/

//! While running a simulation, we often want to have control over some
//! simulation parameters: the temperature, the pressure, etc. This is the goal
//! of the control algorithms, all implmenting of the `Control` trait.

use ::universe::Universe;

/// Trait for controling some parameters in an universe during a simulation.
pub trait Control {
    /// Function called once at the beggining of the simulation, which allow
    /// for some setup of the control algorithm if needed.
    fn setup(&mut self, _: &Universe) {}

    /// Do your job, control algorithm!
    fn control(&mut self, universe: &mut Universe);

    /// Function called once at the end of the simulation.
    fn finish(&mut self, _: &Universe) {}
}
