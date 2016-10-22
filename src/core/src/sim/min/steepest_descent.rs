// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

use utils;
use sys::System;
use sim::{Propagator, TemperatureStrategy};

use std::f64;

/// Steepest gradient descent for energy minimization.
///
/// This algorithm is very rough, and will not converge in all the situations.
/// However it is easy to use, and simple enough to be implemented quickly.
pub struct SteepestDescent {
    /// Dumping factor
    gamma: f64,
    /// Force norm convergence criterion
    force_crit: f64,
    /// Energy convergence criterion
    energy_crit: f64,
    /// Has the last minimization converged?
    is_converged: bool
}

impl Default for SteepestDescent {
    fn default() -> SteepestDescent {
        SteepestDescent::new()
    }
}

impl SteepestDescent {
    /// Create a `SteepestDescent` with sensible default values for energy and
    /// force convergence criteria.
    ///
    /// The default for force criterion is `1e-5 kJ/mol/Å^2`, and `1e-5
    /// kJ/mol/Å^2` for the energy criterion.
    pub fn new() -> SteepestDescent {
        let delta_f = utils::unit_from(1e-5, "kJ/mol/A");
        let delta_e = utils::unit_from(1e-5, "kJ/mol");
        SteepestDescent::with_criteria(delta_f, delta_e)
    }

    /// Create a new `SteepestDescent` with the given `force` and `energy`
    /// convergence criterion.
    pub fn with_criteria(force: f64, energy: f64) -> SteepestDescent {
        let gamma = utils::unit_from(0.1, "fs^2/u");
        SteepestDescent{
            gamma: gamma,
            energy_crit: energy,
            force_crit: force*force,
            is_converged: false
        }
    }

    /// Has the minimization converged so far ?
    pub fn converged(&self) -> bool {
        self.is_converged
    }
}

impl Propagator for SteepestDescent {
    fn temperature_strategy(&self) -> TemperatureStrategy {
        TemperatureStrategy::None
    }

    fn setup(&mut self, _: &System) {
        self.is_converged = false;
    }

    fn propagate(&mut self, system: &mut System) {
        let forces = system.forces();
        // Find the largest non-NaN in forces, or NaN otherwise
        let max_force_norm2 = forces.iter().map(|&f| f.norm2()).fold(f64::NAN, f64::max);
        if max_force_norm2 < self.force_crit {
            self.is_converged = true;
            return;
        }

        let energy = system.potential_energy();
        if energy < self.energy_crit {
            self.is_converged = true;
            return;
        }

        // Store the current coordinates
        let mut positions = Vec::with_capacity(system.size());
        for particle in system.iter() {
            positions.push(particle.position);
        }
        let positions = positions;

        let mut gamma_changed = false;
        // Update coordinates, reducing gamma until we find a configuration of
        // lower energy
        loop {
            for (i, p) in system.iter_mut().enumerate() {
                p.position = positions[i] + self.gamma * forces[i];
            }
            let new_energy = system.potential_energy();
            if new_energy <= energy {
                break;
            }
            self.gamma /= 2.0;
            gamma_changed = true;
        }

        // If we had a successful iteration without needing to reduce gamma,
        // we can increase it slightly
        if !gamma_changed {
            self.gamma *= 1.1;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use sys::*;
    use types::*;
    use energy::*;
    use sim::Propagator;

    #[test]
    fn minization() {
        let mut system = System::from_cell(UnitCell::cubic(20.0));;
        system.add_particle(Particle::new("Cl"));
        system[0].position = Vector3D::zero();
        system.add_particle(Particle::new("Cl"));
        system[1].position = Vector3D::new(0.0, 0.0, 2.0);

        system.interactions_mut().add_pair("Cl", "Cl",
            PairInteraction::new(Box::new(Harmonic{x0: 2.3, k: 0.1}), 10.0)
        );

        let mut minization = SteepestDescent::new();
        for _ in 0..100 {
            minization.propagate(&mut system);
        }
        assert!(minization.converged());
        assert_approx_eq!(system.distance(0, 1), 2.3, 1e-4);
    }
}
