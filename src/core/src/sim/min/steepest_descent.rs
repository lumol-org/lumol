// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

use sys::System;
use sys::zip_particle::*;
use utils;

use super::{Minimizer, Tolerance};

use std::f64;

/// Steepest descent minimization algorithm.
///
/// This method propagates the system along the gradient of energy to find a
/// minimum. Although easy to use, it will not converge in all situations.
pub struct SteepestDescent {
    /// Damping factor
    gamma: f64,
}

impl SteepestDescent {
    /// Create a new `SteepestDescent` minimizer
    pub fn new() -> SteepestDescent {
        SteepestDescent {
            gamma: utils::unit_from(0.1, "fs^2/u"),
        }
    }
}

impl Minimizer for SteepestDescent {
    fn minimize(&mut self, system: &mut System) -> Tolerance {
        // Store the current coordinates
        let prevpos = system.particles().position.to_vec();

        let mut gamma_changed = false;
        let forces = system.forces();
        let initial_energy = system.potential_energy();
        let mut energy;
        // Update coordinates, reducing gamma until we find a configuration of
        // lower energy
        loop {
            for (position, prevpos, force) in
                system.particles_mut().zip_mut((&mut Position, &prevpos, &forces))
            {
                *position = prevpos + self.gamma * force;
            }

            energy = system.potential_energy();
            if energy <= initial_energy {
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

        return Tolerance {
            energy: energy,
            force2: forces.iter().map(|&f| f.norm2()).fold(f64::NAN, f64::max),
        };
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use energy::{Harmonic, PairInteraction};
    use sim::{Minimization, Propagator};
    use sys::System;
    use utils::system_from_xyz;

    fn testing_system() -> System {
        let mut system = system_from_xyz(
            "2
            cell: 20.0
            Cl 0.0 0.0 0.0
            Cl 0.0 0.0 2.0
            ",
        );

        let pair = PairInteraction::new(Box::new(Harmonic { x0: 2.3, k: 0.1 }), 10.0);
        system.add_pair_potential(("Cl", "Cl"), pair);
        return system;
    }

    #[test]
    fn minization() {
        let mut system = testing_system();

        let mut minization = Minimization::new(Box::new(SteepestDescent::new()));
        for _ in 0..100 {
            minization.propagate(&mut system);
        }
        assert!(minization.converged());
        assert_relative_eq!(system.distance(0, 1), 2.3, epsilon = 1e-3);
    }
}
