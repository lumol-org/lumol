// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

//! Energy minimization algorithms
use utils;
use sys::System;
use sim::{Propagator, TemperatureStrategy};

use std::f64;

/// Criteria used for energy minimization
pub struct EnergyCriteria {
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
    /// energy criteria.
    fn minimize(&mut self, system: &mut System) -> EnergyCriteria;
}

/// Minimization propagator for simulations.
///
/// The minimization stops when the energy difference between the previous and
/// the current step is lower than the energy criterion, or when the maximal
/// squared norm of the atomic force is lower than the force criterion.
pub struct Minimization {
    minimizer: Box<Minimizer>,
    is_converged: bool,
    last_energy: f64,
    criteria: EnergyCriteria
}

impl Minimization {
    /// Create a new `Minimization` using the given `minimizer`.
    pub fn new(minimizer: Box<Minimizer>) -> Minimization {
        let critera = EnergyCriteria {
            energy: utils::unit_from(1e-5, "kJ/mol"),
            force2: utils::unit_from(1e-5, "kJ^2/mol^2/A^2"),
        };
        return Minimization::with_criteria(minimizer, critera);
    }

    /// Create a new `Minimization` using the given `minimizer` and specific
    /// energy and force `criteria`.
    pub fn with_criteria(minimizer: Box<Minimizer>, criteria: EnergyCriteria) -> Minimization {
        Minimization {
            minimizer: minimizer,
            is_converged: false,
            last_energy: 0.0,
            criteria: criteria
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

        if result.force2 < self.criteria.force2 {
            self.is_converged = true;
            info!("Minimization converged on force tolerance");
        }

        if (self.last_energy - result.energy).abs() < self.criteria.energy {
            self.is_converged = true;
            info!("Minimization converged on energy tolerance");
        }

        self.last_energy = result.energy;
    }
}
