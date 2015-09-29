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
use super::Output;

/// The Simulation struct holds all the needed algorithms for running the
/// simulation. It should be use together with an Universe to perform the
/// simulation.
pub struct Simulation {
    propagator: Box<Propagator>,
    outputs: Vec<Box<Output>>
}

impl Simulation {
    /// Create a new Simulation from a Propagator.
    pub fn new<P>(propagator: P) -> Simulation where P: Propagator + 'static {
        Simulation {
            propagator: Box::new(propagator),
            outputs: Vec::new(),
        }
    }

    /// Run the simulation on Universe for `nsteps` steps.
    pub fn run(&mut self, universe: &mut Universe, nsteps: usize) {
        self.setup(universe);
        for _ in 0..nsteps {
            self.propagator.propagate(universe);
            universe.increment_step();
            for output in &mut self.outputs {
                output.write(universe);
            }
        }
        self.finish(universe);
    }

    /// Add a new `Output` algorithm in the outputs list
    pub fn add_output<O>(&mut self, out: O) where O: Output + 'static {
        self.outputs.push(Box::new(out));
    }

    fn setup(&mut self, universe: &mut Universe) {
        self.propagator.setup(universe);
        for output in &mut self.outputs {
            output.setup(universe);
        }
    }

    fn finish(&mut self, universe: &mut Universe) {
        self.propagator.finish(universe);
        for output in &mut self.outputs {
            output.finish(universe);
        }
    }
}
