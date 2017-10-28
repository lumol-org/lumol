// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

//! Metropolis Monte Carlo propagator implementation
use rand::{self, SeedableRng};

use consts::K_BOLTZMANN;
use sys::{System, EnergyCache};
use sim::{Propagator, TemperatureStrategy};

use super::MCMove;

/// Metropolis Monte Carlo propagator
pub struct MonteCarlo {
    /// Boltzmann factor: beta = 1/(kB * T)
    beta: f64,
    /// List of possible Monte Carlo moves
    moves: Vec<(Box<MCMove>, MoveCounter)>,
    /// Cummulative frequencies of the Monte Carlo moves
    frequencies: Vec<f64>,
    /// Specifies the number of moves after which an update of a move's
    /// amplitude is performed.
    update_frequency: u64,
    /// Random number generator for the simulation. All random state will be
    /// taken from this.
    rng: Box<rand::Rng>,
    /// Cache for faster energy computation
    cache: EnergyCache,
    /// Flag checking if the moves frequencies has been converted to
    /// cumulative frequencies or not yet.
    initialized: bool,
}

impl MonteCarlo {
    /// Create a new Monte Carlo propagator at temperature `T`.
    pub fn new(temperature: f64) -> MonteCarlo {
        let mut rng = Box::new(rand::XorShiftRng::new_unseeded());
        rng.reseed([2015u32, 42u32, 3u32, 12u32]);
        return MonteCarlo::from_rng(temperature, rng);
    }

    /// Create a Monte Carlo propagator at temperature `T`, using the `rng`
    /// random number generator.
    pub fn from_rng(temperature: f64, rng: Box<rand::Rng>) -> MonteCarlo {
        assert!(temperature >= 0.0,
                "Monte Carlo temperature must be positive");
        MonteCarlo {
            beta: 1.0 / (K_BOLTZMANN * temperature),
            moves: Vec::new(),
            frequencies: Vec::new(),
            update_frequency: 0,
            rng: rng,
            cache: EnergyCache::new(),
            initialized: false,
        }
    }

    /// Add the `mcmove` Monte Carlo move to this propagator, with frequency
    /// `frequency`. All calls to this function should happen before any
    /// simulation run.
    ///
    /// # Panics
    ///
    /// If called after a simulation run.
    pub fn add(&mut self, mcmove: Box<MCMove>, frequency: f64) {
        if self.initialized {
            fatal_error!("Monte Carlo simulation has already been initialized, we can not add \
                          new moves.");
        }
        self.moves.push((mcmove, MoveCounter::new(None)));
        self.frequencies.push(frequency);
    }

    /// Add the `mcmove` Monte Carlo move to the propagator.
    /// `frequency` describes how frequent a move is called, `target_acceptance`
    /// is the desired acceptance ratio of the move.
    ///
    /// # Panics
    ///
    /// If called after a simulation run.
    /// If `target_acceptance` is either negative or larger than one.
    pub fn add_move_with_acceptance(&mut self,
                                    mcmove: Box<MCMove>,
                                    frequency: f64,
                                    target_acceptance: f64) {
        if self.initialized {
            fatal_error!("Monte Carlo simulation has already been initialized, we can not add \
                          new moves.");
        }
        self.moves.push((mcmove, MoveCounter::new(Some(target_acceptance))));
        self.frequencies.push(frequency);
    }

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

    fn normalize_frequencies(&mut self) {
        assert_eq!(self.frequencies.len(), self.moves.len());
        if self.frequencies.is_empty() {
            warn!("No move in the Monte Carlo simulation, did you forget to specify them?");
            return;
        }

        if self.initialized {
            error!("This Monte Carlo simulation has already been initialized.");
            return;
        }

        self.initialized = true;
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
}

impl Propagator for MonteCarlo {
    fn temperature_strategy(&self) -> TemperatureStrategy {
        TemperatureStrategy::External(self.temperature())
    }

    fn setup(&mut self, system: &System) {
        self.normalize_frequencies();
        self.cache.init(system);
        for mc_move in &mut self.moves {
            mc_move.0.setup(system)
        }
    }

    fn propagate(&mut self, system: &mut System) {
        let mcmove = {
            let probability = self.rng.next_f64();
            // Get the index of the first move with frequency >= probability.
            let (i, _) = self.frequencies
                .iter()
                .enumerate()
                .find(|&(_, f)| probability <= *f)
                .expect("Could not find a move in MonteCarlo moves list");
            &mut self.moves[i]
        };
        trace!("Selected move is '{}'", mcmove.0.describe());

        if !mcmove.0.prepare(system, &mut self.rng) {
            trace!("    --> Can not perform the move");
            return;
        }

        // compute cost
        let cost = mcmove.0.cost(system, self.beta, &mut self.cache);
        trace!("    --> Move cost is {}", cost);

        // apply metropolis criterion
        let accepted = cost <= 0.0 || self.rng.next_f64() < f64::exp(-cost);

        if accepted {
            trace!("    --> Move was accepted");
            mcmove.0.apply(system);
            self.cache.update(system);
            mcmove.1.accept();
        } else {
            trace!("    --> Move was rejected");
            mcmove.0.restore(system);
            mcmove.1.reject();
        }

        // Do the adjustments for the selected move as needed
        if mcmove.1.nattempted == self.update_frequency {
            mcmove.0.update_amplitude(mcmove.1.compute_scaling_factor());
            mcmove.1.update();
        }
    }

    /// Print some informations about moves to screen
    fn finish(&mut self, _: &System) {
        info!("Monte Carlo simulation summary");
        for mc_move in &self.moves {
            info!("Statistics for move: {}", mc_move.0.describe());
            mc_move.1.print_info();
        }
    }
}

/// This struct keeps track of the number of times a move was called
/// and how often it was accepted.
///
/// These statistics can be used to compute a scaling factor to
/// increase the efficiency of a move.
///
/// # Example
///
/// ```
/// use lumol_core::sim::mc::MoveCounter;
///
/// // Create a new counter with a target acceptance of 0.5 (=50 %)
/// let mut counter = MoveCounter::new(Some(0.5));
/// // Move was attempted and accepted.
/// counter.accept();
/// // Increase scaling factor since 100 % of moves where accepted.
/// assert_eq!(counter.compute_scaling_factor(), Some(1.2));
/// // Update the counters to track acceptance for move with new
/// // scaling factor.
/// counter.update();
/// // Directly calling for a new scaling factor returns
/// // `None` since there were no attempts since last update.
/// assert_eq!(counter.compute_scaling_factor(), None);
/// // The overall acceptance is still 1.0 (=100 %).
/// assert_eq!(counter.current_acceptance(), 1.0);
///
/// // The scaling factor is computed based on the acceptance
/// // since the last update was performed, not the overall acceptance
/// counter.reject();
/// // The scaling factor will be smaller: the acceptance
/// // since the last update is zero.
/// assert_eq!(counter.compute_scaling_factor(), Some(0.8));
/// assert_eq!(counter.current_acceptance(), 0.5);
/// ```
pub struct MoveCounter {
    /// Count the total number of times the move was called.
    pub nattempted_total: u64,
    /// Count the total number of times the move was accepted.
    pub naccepted_total: u64,
    /// Count the number of times the move was accepted since the last update.
    pub naccepted: u64,
    /// Count the number of times the move was called since the last update.
    pub nattempted: u64,
    /// The target fraction of accepted over attempted moves.
    target_acceptance: Option<f64>,
}

impl MoveCounter {
    /// Create a new counter for the move, initializing all counts to zero and
    /// setting the `target_acceptance`.
    pub fn new(target_acceptance: Option<f64>) -> MoveCounter {
        let mut counter = MoveCounter {
            nattempted_total: 0,
            naccepted_total: 0,
            naccepted: 0,
            nattempted: 0,
            target_acceptance: None,
        };
        counter.set_acceptance(target_acceptance);
        counter
    }

    /// Set the target acceptance for the move counter.
    pub fn set_acceptance(&mut self, target_acceptance: Option<f64>) {
        // Check if `target_acceptance` has a valid value.
        if let Some(acceptance) = target_acceptance {
            assert!(0.0 < acceptance && acceptance < 1.0,
                    "The target acceptance ratio has to be a value between 0 and 1")
        }
        self.target_acceptance = target_acceptance;
    }

    /// Increase counters for attempt.
    #[inline]
    pub fn reject(&mut self) {
        self.nattempted_total += 1;
        self.nattempted += 1;
    }

    /// Increase counters to track the number of times the move was accepted.
    #[inline]
    pub fn accept(&mut self) {
        self.nattempted_total += 1;
        self.nattempted += 1;
        self.naccepted += 1;
        self.naccepted_total += 1;
    }

    /// Reset counters for attempts and acceptance since the last update.
    #[inline]
    pub fn update(&mut self) {
        self.naccepted = 0;
        self.nattempted = 0;
    }

    /// Print the total number of attempts and the current acceptance to `info` log.
    pub fn print_info(&self) {
        info!("  Attempts  : {}", self.nattempted_total);
        info!("  Acceptance: {:.4} %", self.current_acceptance() * 100.0);
    }

    /// Return fraction of total number of accepted over total number of attempted moves.
    pub fn current_acceptance(&self) -> f64 {
        if self.nattempted_total != 0 {
            self.naccepted_total as f64 / self.nattempted_total as f64
        } else {
            0.0
        }
    }

    /// Compute a scaling factor according to the desired acceptance.
    pub fn compute_scaling_factor(&self) -> Option<f64> {
        // Check if there exists an target_acceptance
        if let Some(ta) = self.target_acceptance {
            // Capture division by zero
            if self.nattempted == 0 {
                return None;
            };
            let quotient = self.naccepted as f64 / self.nattempted as f64 / ta;
            // Limit the change
            match quotient {
                _ if quotient > 1.2 => Some(1.2),
                _ if quotient < 0.8 => Some(0.8),
                _ => Some(quotient),
            }
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use sim::mc::{MonteCarlo, MCMove, MoveCounter};
    use sim::Propagator;
    use sys::{System, EnergyCache};
    use rand::Rng;

    struct DummyMove;
    impl MCMove for DummyMove {
        fn describe(&self) -> &str {
            "dummy"
        }
        fn setup(&mut self, _: &System) {}
        fn prepare(&mut self, _: &mut System, _: &mut Box<Rng>) -> bool {
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

    #[test]
    #[should_panic]
    fn too_large_acceptance() {
        let mut mc = MonteCarlo::new(100.0);
        mc.add_move_with_acceptance(Box::new(DummyMove), 1.0, 0.5);
        mc.moves[0].1.set_acceptance(Some(1.1));
    }

    #[test]
    #[should_panic]
    fn negative_acceptance() {
        let mut mc = MonteCarlo::new(100.0);
        mc.add_move_with_acceptance(Box::new(DummyMove), 1.0, -0.1);
    }

    #[test]
    fn valid_acceptance() {
        let mut mc = MonteCarlo::new(100.0);
        mc.add_move_with_acceptance(Box::new(DummyMove), 1.0, 0.5);
        assert_eq!(mc.moves[0].1.target_acceptance, Some(0.5));
        mc.moves[0].1.set_acceptance(None);
        assert_eq!(mc.moves[0].1.target_acceptance, None);
    }

    #[test]
    fn scaling_factor() {
        let mut counter = MoveCounter::new(Some(0.5));
        assert_eq!(counter.compute_scaling_factor(), None);
        counter.nattempted = 100;
        counter.naccepted = 100;
        assert_eq!(counter.compute_scaling_factor(), Some(1.2));
        counter.naccepted = 0;
        assert_eq!(counter.compute_scaling_factor(), Some(0.8));
        counter.naccepted = 55;
        assert_eq!(counter.compute_scaling_factor(), Some(1.1));
    }
}
