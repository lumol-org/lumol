// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

//! Metropolis Monte-Carlo propagator implementation
use rand::{self, SeedableRng};

use consts::K_BOLTZMANN;
use sys::{System, EnergyCache};
use sim::{Propagator, TemperatureStrategy};

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
    /// Flag checking if the moves frequencies has been converted to
    /// cummulative frequencies or not yet.
    initialized: bool,
}

impl MonteCarlo {
    /// Create a new Monte-Carlo propagator at temperature `T`.
    pub fn new(temperature: f64) -> MonteCarlo {
        let mut rng = Box::new(rand::XorShiftRng::new_unseeded());
        rng.reseed([2015u32, 42u32, 3u32, 12u32]);
        return MonteCarlo::from_rng(temperature, rng);
    }

    /// Create a Monte-Carlo propagator at temperature `T`, using the `rng`
    /// random number generator.
    pub fn from_rng(temperature: f64, rng: Box<rand::Rng>) -> MonteCarlo {
        assert!(temperature >= 0.0, "Monte-Carlo temperature must be positive");
        MonteCarlo {
            beta: 1.0 / (K_BOLTZMANN * temperature),
            moves: Vec::new(),
            frequencies: Vec::new(),
            rng: rng,
            cache: EnergyCache::new(),
            initialized: false,
        }
    }

    /// Add a the `mcmove` Monte-Carlo move to this propagator, with frequency
    /// `frequency`. All calls to this function should happen before any
    /// simulation run.
    ///
    /// # Panics
    ///
    /// If called after a simulation run.
    pub fn add(&mut self, mcmove: Box<MCMove>, frequency: f64) {
        if self.initialized {
            fatal_error!(
                "Monte-Carlo simulation has already been initialized, \
                we can not add new moves."
            );
        }
        self.moves.push(mcmove);
        self.frequencies.push(frequency);
    }

    /// Get the temperature of the simulation
    pub fn temperature(&self) -> f64 {
        1.0 / (self.beta * K_BOLTZMANN)
    }

    /// Set the temperature of the simulation
    pub fn set_temperature(&mut self, temperature: f64) {
        self.beta = 1.0 / (temperature * K_BOLTZMANN);
    }

    fn normalize_frequencies(&mut self) {
        assert!(self.frequencies.len() == self.moves.len());
        if self.frequencies.is_empty() {
            warn!(
                "No move in the Monte-Carlo simulation, \
                did you forget to specify them?"
            );
            return;
        }

        if self.initialized {
            error!("This Monte-Carlo simulation has already been initialized.");
            return;
        }

        self.initialized = true;
        // Normalize the frequencies
        let sum = self.frequencies.iter().fold(0.0, |sum, &f| sum + f);
        for frequency in &mut self.frequencies {
            *frequency /= sum;
        }
        // Make the frequencies vector contain cummulative frequencies
        for i in 1..self.frequencies.len() {
            self.frequencies[i] += self.frequencies[i - 1];
        }
        let last = self.frequencies.len() - 1;
        self.frequencies[last] = 1.0;
    }
}

impl Propagator for MonteCarlo {
    fn temperature_strategy(&self) -> TemperatureStrategy {
        TemperatureStrategy::External(self.temperature())
    }

    fn setup(&mut self, system: &System) {
        self.normalize_frequencies();
        self.cache.init(system);
    }

    fn propagate(&mut self, system: &mut System) {
        let mcmove = {
            let probability = self.rng.next_f64();
            // Get the index of the first move with frequency >= probability.
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

#[cfg(test)]
mod tests {
    use sim::mc::{MonteCarlo, MCMove};
    use sim::Propagator;
    use sys::{System, EnergyCache};
    use rand::Rng;

    struct DummyMove;
    impl MCMove for DummyMove {
        fn describe(&self) -> &str {"dummy"}
        fn prepare(&mut self, _: &mut System, _: &mut Box<Rng>) -> bool {true}
        fn cost(&self, _: &System, _: f64, _: &mut EnergyCache) -> f64 {0.0}
        fn apply(&mut self, _: &mut System) {}
        fn restore(&mut self, _: &mut System) {}
    }

    #[test]
    fn frequencies() {
        let mut mc = MonteCarlo::new(100.0);
        mc.add(Box::new(DummyMove), 13.0);
        mc.add(Box::new(DummyMove), 2.0);
        mc.add(Box::new(DummyMove), 5.0);

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
    #[should_panic]
    fn add_after_init() {
        let mut mc = MonteCarlo::new(100.0);
        mc.add(Box::new(DummyMove), 1.0);
        mc.setup(&System::new());
        mc.add(Box::new(DummyMove), 1.0);
    }
}
