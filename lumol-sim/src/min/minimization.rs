// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

//! Energy minimization algorithms
use std::f64;

use log::info;

use lumol_core::{System, DegreesOfFreedom};

use crate::propagator::{Propagator, TemperatureStrategy};

/// Tolerance criteria used for energy minimization
pub struct Tolerance {
    /// Potential energy of the system
    pub energy: f64,
    /// Maximal squared norm of the force acting on an atom
    pub force2: f64,
}

/// The `Minimizer` trait define minimization interface.
///
/// A minimizer is an algorithm responsible for finding new configurations of
/// lower energy.
pub trait Minimizer {
    /// Setup the minimizer. This function is called once at the begining of
    /// every simulation run.
    fn setup(&mut self, _: &System) {}

    /// Find a new configuration of lower energy, and return the corresponding
    /// values for energy and forces.
    fn minimize(&mut self, system: &mut System) -> Tolerance;

    /// Get the number of degrees of freedom simulated by this minimizer
    ///
    /// This function is called once at thr beginning of the simulation
    fn degrees_of_freedom(&self, system: &System) -> DegreesOfFreedom;
}

/// Minimization propagator for simulations.
///
/// The minimization stops when the energy difference between the previous and
/// the current step is lower than the energy criterion, or when the maximal
/// squared norm of the atomic force is lower than the force criterion.
pub struct Minimization {
    minimizer: Box<dyn Minimizer>,
    is_converged: bool,
    last_energy: f64,
    tolerance: Tolerance,
}

impl Minimization {
    /// Create a new `Minimization` using the given `minimizer` and specific
    /// energy and force `tolerance`.
    pub fn new(minimizer: Box<dyn Minimizer>, tolerance: Tolerance) -> Minimization {
        Minimization {
            minimizer: minimizer,
            is_converged: false,
            last_energy: 0.0,
            tolerance: tolerance,
        }
    }

    /// Check if the minimization has converged.
    pub fn converged(&self) -> bool {
        self.is_converged
    }
}

impl Propagator for Minimization {
    fn temperature_strategy(&self) -> TemperatureStrategy {
        TemperatureStrategy::None
    }

    fn degrees_of_freedom(&self, system: &System) -> DegreesOfFreedom {
        self.minimizer.degrees_of_freedom(system)
    }

    fn setup(&mut self, system: &System) {
        self.is_converged = false;
        self.last_energy = system.potential_energy();
        self.minimizer.setup(system);
    }

    fn propagate(&mut self, system: &mut System) {
        if self.is_converged {
            return;
        }

        let result = self.minimizer.minimize(system);

        if result.force2 < self.tolerance.force2 {
            self.is_converged = true;
            info!("Minimization converged on force tolerance");
        }

        if (self.last_energy - result.energy).abs() < self.tolerance.energy {
            self.is_converged = true;
            info!("Minimization converged on energy tolerance");
        }

        self.last_energy = result.energy;
    }
}
