/*
 * Cymbalum, Molecular Simulation in Rust
 * Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
*/

use ::universe::Universe;

use super::Propagator;

/// The Simulation struct holds all the needed algorithms for running the
/// simulation. It should be use together with an Universe to perform the
/// simulation.
pub struct Simulation {
    propagator: Box<Propagator>,
}

impl Simulation {
    /// Create a new Simulation from a Propagator.
    pub fn new<P>(propagator: P) -> Simulation where P: Propagator + 'static {
        Simulation {
            propagator: Box::new(propagator),
        }
    }

    /// Run the simulation on Universe for `nsteps` steps.
    pub fn run(&mut self, universe: &mut Universe, nsteps: usize) {
        self.propagator.setup(universe);
        for _ in 0..nsteps {
            self.propagator.propagate(universe);
        }
        self.propagator.finish(universe);
    }
}
