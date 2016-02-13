// Cymbalum, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

//! Metropolis Monte-Carlo propagator implementation
extern crate rand;
use self::rand::SeedableRng;

use constants::K_BOLTZMANN;
use system::{System, EnergyCache};
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
    /// Cache for faster energy computation
    cache: EnergyCache,
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
            rng: rng,
            cache: EnergyCache::new(),
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
    fn setup(&mut self, system: &System) {
        self.normalize_frequencies();
        self.cache.init(system);
    }

    fn propagate(&mut self, system: &mut System) {
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

        if !mcmove.prepare(system, &mut self.rng) {
            trace!("    --> Can not perform the move");
            return;
        }

        let cost = mcmove.cost(system, self.beta, &mut self.cache);
        trace!("    --> Move cost is {}", cost);

        let accepted = cost <= 0.0 || self.rng.next_f64() < f64::exp(-cost);

        if accepted {
            trace!("    --> Move was accepted");
            mcmove.apply(system);
            self.cache.update(system);
        } else {
            trace!("    --> Move was rejected");
            mcmove.restore(system);
        }
    }
}
