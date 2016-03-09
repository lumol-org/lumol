// Cymbalum, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

use units;
use System;
use simulation::Propagator;

/// Steepest gradient descent for energy minization.
///
/// This algorithm is very rough, and will not converge in all the situations.
/// However it is easy to use, and simple enough to be implemented quickly.
pub struct SteepestDescent {
    /// Dumping factor
    gamma: f64,
    /// Force norm convergence criterium
    force_crit: f64,
    /// Energy convergence criterium
    energy_crit: f64,
    /// Has the last minization converged?
    is_converged: bool
}

impl SteepestDescent {
    /// Create a GradientDescent with sensible default values for energy and
    /// force convergence criteria.
    ///
    /// The default for force criterium is `1e-5 kJ/mol/Å^2`, and `1e-5
    /// kJ/mol/Å^2` for the energy criterium.
    pub fn new() -> SteepestDescent {
        let delta_f = units::from(1e-5, "kJ/mol/A").unwrap();
        let delta_e = units::from(1e-5, "kJ/mol").unwrap();
        SteepestDescent::with_criteria(delta_f, delta_e)
    }

    /// Create a new `GradientDescent` with the force convergence criterium of
    /// `force`, and the energy convergence criterium of `energy`.
    pub fn with_criteria(force: f64, energy: f64) -> SteepestDescent {
        let gamma = units::from(0.1, "fs^2/u").unwrap();
        SteepestDescent{
            gamma: gamma,
            energy_crit: energy,
            force_crit: force*force,
            is_converged: false
        }
    }

    /// Has the minization converged so far ?
    pub fn converged(&self) -> bool {
        self.is_converged
    }
}

impl Propagator for SteepestDescent {
    fn setup(&mut self, _: &System) {
        self.is_converged = false;
    }

    fn propagate(&mut self, system: &mut System) {
        let forces = system.forces();
        // Get the maximal value in the vector.
        // How can you know that fold(cte, cte) will do it ?
        let max_force_norm2 = forces.iter().map(|&f| f.norm2()).fold(-1./0. /* -inf */, f64::max);
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
    use system::*;
    use types::*;
    use potentials::*;
    use simulation::Propagator;

    #[test]
    fn minization() {
        let mut system = System::from_cell(UnitCell::cubic(20.0));;
        system.add_particle(Particle::new("Cl"));
        system[0].position = Vector3D::new(0.0, 0.0, 0.0);
        system.add_particle(Particle::new("Cl"));
        system[1].position = Vector3D::new(0.0, 0.0, 2.0);
        system.add_pair_interaction("Cl", "Cl", Box::new(Harmonic{x0: 2.3, k: 0.1}));

        let mut minization = SteepestDescent::new();
        for _ in 0..100 {
            minization.propagate(&mut system);
        }
        assert!(minization.converged());
        assert_approx_eq!(system.distance(0, 1), 2.3, 1e-4);
    }
}
