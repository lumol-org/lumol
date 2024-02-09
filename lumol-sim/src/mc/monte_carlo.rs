// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

//! Metropolis Monte Carlo propagator implementation
use std::ops::{Deref, DerefMut};

use rand::{self, Rng, SeedableRng};

use log::{warn, info, trace};

use lumol_core::consts::K_BOLTZMANN;
use lumol_core::{DegreesOfFreedom, EnergyCache, System};

use crate::propagator::{Propagator, TemperatureStrategy};
use super::{MCDegreeOfFreedom, MCMove};

/// This struct keeps a move and some statistics on the move (number of times
/// it was called, how often it was accepted, ...)
///
/// These statistics can be used to compute a scaling factor to increase the
/// efficiency of a move.
struct Move {
    /// The actual move implementation
    mcmove: Box<dyn MCMove>,
    /// Count the total number of times the move was called.
    total_attempted: u64,
    /// Count the total number of times the move was accepted.
    total_accepted: u64,
    /// Count the number of times the move was accepted since the last update.
    accepted: u64,
    /// Count the number of times the move was called since the last update.
    attempted: u64,
    /// The target fraction of accepted over attempted moves.
    target_acceptance: Option<f64>,
}

impl Move {
    /// Create a move, initializing all counts to zero and setting the
    /// `target_acceptance`.
    pub fn new(mcmove: Box<dyn MCMove>, target_acceptance: Option<f64>) -> Move {
        let mut counter = Move {
            mcmove: mcmove,
            total_attempted: 0,
            total_accepted: 0,
            accepted: 0,
            attempted: 0,
            target_acceptance: None,
        };
        counter.set_acceptance(target_acceptance);
        counter
    }

    /// Set the target acceptance for the move.
    pub fn set_acceptance(&mut self, target_acceptance: Option<f64>) {
        // Check if `target_acceptance` has a valid value.
        if let Some(acceptance) = target_acceptance {
            assert!(
                0.0 < acceptance && acceptance < 1.0,
                "target acceptance ratio must be between 0 and 1, got {}",
                acceptance
            );
        }
        self.target_acceptance = target_acceptance;
    }

    /// Increase counters for attempt.
    pub fn reject(&mut self) {
        self.total_attempted += 1;
        self.attempted += 1;
    }

    /// Increase counters to track the number of times the move was accepted.
    pub fn accept(&mut self) {
        self.total_attempted += 1;
        self.attempted += 1;
        self.accepted += 1;
        self.total_accepted += 1;
    }

    /// Reset counters for attempts and acceptance since the last update.
    pub fn update(&mut self) {
        self.accepted = 0;
        self.attempted = 0;
    }

    /// Return fraction of total number of accepted over total number of attempted moves.
    pub fn acceptance(&self) -> f64 {
        if self.total_attempted == 0 {
            0.0
        } else {
            self.total_accepted as f64 / self.total_attempted as f64
        }
    }

    /// Compute a scaling factor according to the desired acceptance.
    pub fn scaling_factor(&self) -> Option<f64> {
        if let Some(ta) = self.target_acceptance {
            if self.attempted == 0 {
                return None;
            };
            let quotient = self.accepted as f64 / self.attempted as f64 / ta;
            // Limit the change to 20%
            if quotient > 1.2 {
                Some(1.2)
            } else if quotient < 0.8 {
                Some(0.8)
            } else {
                Some(quotient)
            }
        } else {
            None
        }
    }
}

impl Deref for Move {
    type Target = Box<dyn MCMove>;
    fn deref(&self) -> &Self::Target {
        &self.mcmove
    }
}

impl DerefMut for Move {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.mcmove
    }
}

/// Metropolis Monte Carlo propagator
pub struct MonteCarlo {
    /// Random number generator for the simulation. All random state will be
    /// taken from this.
    rng: Box<dyn rand::RngCore>,
    /// Boltzmann factor: beta = 1/(kB * T)
    beta: f64,
    /// List of possible Monte Carlo moves
    moves: Vec<Move>,
    /// Cummulative frequencies of the Monte Carlo moves
    frequencies: Vec<f64>,
    /// Specifies the number of moves after which an update of a move's
    /// amplitude is performed.
    update_frequency: u64,
    /// Cache for faster energy computation
    cache: EnergyCache,
}

/// Builder for `MonteCarlo` struct
pub struct MonteCarloBuilder {
    rng: Box<dyn rand::RngCore>,
    beta: f64,
    moves: Vec<Move>,
    frequencies: Vec<f64>,
}

impl MonteCarloBuilder {
    /// Create a new Monte Carlo propagator at temperature `T`.
    pub fn new(temperature: f64) -> MonteCarloBuilder {
        let rng = Box::new(rand_xorshift::XorShiftRng::from_seed([
            0xeb, 0xa8, 0xe4, 0x29, 0xca, 0x60, 0x44, 0xb0,
            0xd3, 0x77, 0xc6, 0xa0, 0x21, 0x71, 0x37, 0xf7,
        ]));
        return MonteCarloBuilder::from_rng(temperature, rng);
    }

    /// Create a Monte Carlo propagator at temperature `T`, using the `rng`
    /// random number generator.
    pub fn from_rng(temperature: f64, rng: Box<dyn rand::RngCore>) -> MonteCarloBuilder {
        assert!(temperature > 0.0, "Monte Carlo temperature must be positive, got {}", temperature);
        MonteCarloBuilder {
            beta: 1.0 / (K_BOLTZMANN * temperature),
            moves: Vec::new(),
            frequencies: Vec::new(),
            rng: rng,
        }
    }

    /// Add the `mcmove` Monte Carlo move to the propagator. `frequency`
    /// describes how frequently the move will be called, and
    /// `target_acceptance` is the desired acceptance ratio of the move.
    pub fn add(
        &mut self,
        mcmove: Box<dyn MCMove>,
        frequency: f64,
        target_acceptance: impl Into<Option<f64>>,
    ) {
        self.moves.push(Move::new(mcmove, target_acceptance.into()));
        self.frequencies.push(frequency);
    }

    fn normalize_frequencies(&mut self) {
        assert_eq!(self.frequencies.len(), self.moves.len());
        if self.frequencies.is_empty() {
            warn!("No move in the Monte Carlo simulation, did you forget to specify them?");
            return;
        }

        // Normalize the frequencies
        let sum = self.frequencies.iter().fold(0.0, |sum, &f| sum + f);
        for frequency in &mut self.frequencies {
            *frequency /= sum;
        }
        // Make the frequencies vector contain cumulative frequencies
        for i in 1..self.frequencies.len() {
            self.frequencies[i] += self.frequencies[i - 1];
        }
        let last = self.frequencies.len() - 1;
        self.frequencies[last] = 1.0;
    }

    /// Normalize the frequencies for all moves and get the corresponding
    /// `MonteCarlo` propagator
    pub fn finish(mut self) -> MonteCarlo {
        self.normalize_frequencies();
        MonteCarlo {
            rng: self.rng,
            beta: self.beta,
            moves: self.moves,
            frequencies: self.frequencies,
            update_frequency: 0,
            cache: EnergyCache::new(),
        }
    }
}

impl MonteCarlo {
    /// Set the number of times a move has to be called before its amplitude
    /// is updated. This value is applied to all moves.
    pub fn set_amplitude_update_frequency(&mut self, frequency: u64) {
        self.update_frequency = frequency;
    }

    /// Get the temperature of the simulation
    pub fn temperature(&self) -> f64 {
        1.0 / (self.beta * K_BOLTZMANN)
    }

    /// Set the temperature of the simulation
    pub fn set_temperature(&mut self, temperature: f64) {
        self.beta = 1.0 / (temperature * K_BOLTZMANN);
    }
}

impl Propagator for MonteCarlo {
    fn temperature_strategy(&self) -> TemperatureStrategy {
        TemperatureStrategy::External(self.temperature())
    }

    fn degrees_of_freedom(&self, system: &System) -> DegreesOfFreedom {
        if self.moves.is_empty() {
            return DegreesOfFreedom::Particles;
        }

        let mut mc_dof = self.moves[0].mcmove.degrees_of_freedom();
        for other in &self.moves[1..] {
            mc_dof = mc_dof.combine(other.mcmove.degrees_of_freedom());
        }

        match mc_dof {
            MCDegreeOfFreedom::Particles => DegreesOfFreedom::Particles,
            MCDegreeOfFreedom::AllMolecules => DegreesOfFreedom::Molecules,
            MCDegreeOfFreedom::Molecules(hashes) => {
                // This is only a warning and not an error, because the
                // composition of the system could change during the simulation
                //
                // They might also be some moves working with molecules that
                // does not exist (yet) in the system.
                let composition = system.composition();
                for (hash, _) in composition.all_molecules() {
                    if !hashes.contains(&hash) {
                        warn!(
                            "the molecules with hash {:?} are not simulated by \
                             this set of Monte Carlo moves",
                            hash
                        );
                    }
                }
                DegreesOfFreedom::Molecules
            }
        }
    }

    fn setup(&mut self, system: &System) {
        self.cache.init(system);
        for mc_move in &mut self.moves {
            mc_move.mcmove.setup(system);
        }
    }

    fn propagate(&mut self, system: &mut System) {
        let current_move = {
            let probability: f64 = self.rng.gen();
            // Get the index of the first move with frequency >= probability.
            let (i, _) = self.frequencies.iter()
                             .enumerate()
                             .find(|&(_, f)| probability <= *f)
                             .expect("Could not find a move in MonteCarlo moves list");
            &mut self.moves[i]
        };
        trace!("Selected move is '{}'", current_move.describe());

        if !current_move.prepare(system, &mut self.rng) {
            trace!("    --> Can not perform the move");
            return;
        }

        // compute cost
        let cost = current_move.cost(system, self.beta, &mut self.cache);
        trace!("    --> Move cost is {}", cost);

        // apply metropolis criterion
        let accepted = cost <= 0.0 || self.rng.gen::<f64>() < f64::exp(-cost);

        if accepted {
            trace!("    --> Move was accepted");
            current_move.apply(system);
            self.cache.update(system);
            current_move.accept();
        } else {
            trace!("    --> Move was rejected");
            current_move.restore(system);
            current_move.reject();
        }

        // Do the adjustments for the selected move as needed
        if current_move.attempted == self.update_frequency {
            let factor = current_move.scaling_factor();
            current_move.update_amplitude(factor);
            current_move.update();
        }
    }

    /// Print some informations about moves to screen
    fn finish(&mut self, _: &System) {
        info!("Monte Carlo simulation summary");
        for mc_move in &self.moves {
            info!(
                "    {}: {} attempts -- {:2.1} % accepted",
                mc_move.describe(),
                mc_move.total_attempted,
                mc_move.acceptance() * 100.0
            );
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::RngCore;
    use crate::propagator::Propagator;
    use crate::mc::{MCDegreeOfFreedom, MCMove};
    use lumol_core::{EnergyCache, System};

    struct DummyMove;
    impl MCMove for DummyMove {
        fn describe(&self) -> &str {
            "dummy"
        }
        fn degrees_of_freedom(&self) -> MCDegreeOfFreedom {
            MCDegreeOfFreedom::Particles
        }
        fn setup(&mut self, _: &System) {}
        fn prepare(&mut self, _: &mut System, _: &mut dyn RngCore) -> bool {
            true
        }
        fn cost(&self, _: &System, _: f64, _: &mut EnergyCache) -> f64 {
            0.0
        }
        fn apply(&mut self, _: &mut System) {}
        fn restore(&mut self, _: &mut System) {}
        fn update_amplitude(&mut self, _: Option<f64>) {}
    }

    #[test]
    fn frequencies() {
        let mut builder = MonteCarloBuilder::new(100.0);
        builder.add(Box::new(DummyMove), 13.0, None);
        builder.add(Box::new(DummyMove), 2.0, None);
        builder.add(Box::new(DummyMove), 5.0, None);

        let mut mc = builder.finish();
        mc.setup(&System::new());
        let mut last_frequency = 0.0;
        for &f in &mc.frequencies {
            assert!(f > last_frequency);
            last_frequency = f;
        }
        assert_eq!(mc.frequencies.last(), Some(&1.0));

        assert_eq!(mc.frequencies.len(), 3);
        assert_eq!(mc.frequencies[0], 0.65);
        assert_eq!(mc.frequencies[1], 0.75);
        assert_eq!(mc.frequencies[2], 1.0);
    }

    #[test]
    #[should_panic(expected = "Monte Carlo temperature must be positive, got -1")]
    fn negative_temperature() {
        let _ = MonteCarloBuilder::new(-1.0);
    }

    #[test]
    #[should_panic(expected = "Monte Carlo temperature must be positive, got 0")]
    fn zero_temperature() {
        let _ = MonteCarloBuilder::new(0.0);
    }

    #[test]
    #[should_panic(expected = "target acceptance ratio must be between 0 and 1, got 1.1")]
    fn too_large_acceptance() {
        MonteCarloBuilder::new(100.0).add(Box::new(DummyMove), 1.0, 1.1);
    }

    #[test]
    #[should_panic(expected = "target acceptance ratio must be between 0 and 1, got -0.1")]
    fn negative_acceptance() {
        MonteCarloBuilder::new(100.0).add(Box::new(DummyMove), 1.0, -0.1);
    }

    #[test]
    fn valid_acceptance() {
        let mut builder = MonteCarloBuilder::new(100.0);
        builder.add(Box::new(DummyMove), 1.0, 0.5);

        let mut mc = builder.finish();
        assert_eq!(mc.moves[0].target_acceptance, Some(0.5));
        mc.moves[0].set_acceptance(None);
        assert_eq!(mc.moves[0].target_acceptance, None);
    }

    #[test]
    fn scaling_factor() {
        let mut counter = Move::new(Box::new(DummyMove), Some(0.5));
        assert_eq!(counter.scaling_factor(), None);
        counter.attempted = 100;
        counter.accepted = 100;
        assert_eq!(counter.scaling_factor(), Some(1.2));
        counter.accepted = 0;
        assert_eq!(counter.scaling_factor(), Some(0.8));
        counter.accepted = 55;
        assert_eq!(counter.scaling_factor(), Some(1.1));
    }
}
