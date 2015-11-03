/* Cymbalum, Molecular Simulation in Rust - Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
 */

//! Metropolis Monte-Carlo propagator implementation
extern crate rand;
use self::rand::SeedableRng;

use constants::K_BOLTZMANN;
use universe::Universe;
use simulation::Propagator;

use super::MCMove;

/// Metropolis Monte-Carlo propagator
pub struct MonteCarlo {
    /// Boltzmann factor: beta = 1/(kB * T)
    beta: f64,
    /// List of possible Monte-Carlo moves
    moves: Vec<Box<MCMove>>,
    /// Cummulative frequencies of the Monte-Carlo moves
    frequencies: Vec<f64>,
    /// Random number generator for the simulation. All random state will be
    /// taken from this.
    rng: Box<rand::Rng>,
}

impl MonteCarlo {
    /// Create a new Monte-Carlo propagator at temperature `T`.
    pub fn new(T: f64) -> MonteCarlo {
        let mut rng = Box::new(rand::XorShiftRng::new_unseeded());
        rng.reseed([2015u32, 42u32, 3u32, 12u32]);
        return MonteCarlo::from_rng(T, rng);
    }

    /// Create a Monte-Carlo propagator at temperature `T`, using the `rng`
    /// random number generator.
    pub fn from_rng(T: f64, rng: Box<rand::Rng>) -> MonteCarlo {
        assert!(T >= 0.0, "Monte-Carlo temperature must be positive");
        MonteCarlo {
            beta: 1.0 / (K_BOLTZMANN * T),
            moves: Vec::new(),
            frequencies: Vec::new(),
            rng: rng
        }
    }

    /// Add a the `mcmove` Monte-Carlo move to this propagator, with frequency
    /// `frequency`.
    pub fn add(&mut self, mcmove: Box<MCMove>, frequency: f64) -> &MonteCarlo {
        self.moves.push(mcmove);
        self.frequencies.push(frequency);
        return self;
    }

    /// Get the temperature of the simulation
    pub fn temperature(&self) -> f64 {
        1.0 / (self.beta * K_BOLTZMANN)
    }

    /// Set the temperature of the simulation
    pub fn set_temperature(&mut self, T: f64) {
        self.beta = 1.0 / (T * K_BOLTZMANN);
    }

    fn normalize_frequencies(&mut self) {
        assert!(self.frequencies.len() == self.moves.len());
        let mut frequency_sum = 0.0;
        for frequency in &self.frequencies {
            frequency_sum += *frequency;
        }

        for frequency in &mut self.frequencies {
            *frequency /= frequency_sum;
        }

        let last = self.frequencies.len() - 1;
        self.frequencies[last] = 1.0;
    }
}

impl Propagator for MonteCarlo {
    fn setup(&mut self, _: &Universe) {
        self.normalize_frequencies();
    }

    fn propagate(&mut self, universe: &mut Universe) {
        let mcmove = {
            let probability = self.rng.next_f64();
            // Get the index of the first move with frequency >= move_f.
            let (i, _) = self.frequencies.iter()
                                         .enumerate()
                                         .find(|&(_, f)| probability <= *f)
                                         .expect("Could not find a move in MonteCarlo moves list");
            &mut self.moves[i]
        };
        trace!("Selected move is '{}'", mcmove.describe());

        mcmove.prepare(universe, &mut self.rng);

        let cost = mcmove.cost(universe, self.beta);
        trace!("    --> Move cost is {}", cost);

        let accepted = if cost <= 0.0 {
            true
        } else {
            let r = self.rng.next_f64();
            if r < f64::exp(-cost) {
                true
            } else {
                false
            }
        };

        if accepted {
            trace!("    --> Move was accepted");
            mcmove.apply(universe);
        } else {
            trace!("    --> Move was rejected");
            mcmove.restore(universe);
        }
    }
}
