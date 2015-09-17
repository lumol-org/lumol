/*
 * Cymbalum, Molecular Simulation in Rust
 * Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
*/
//! Gradient-descent based energy minization.

use ::units;
use ::Universe;
use super::Propagator;
use super::{Compute, PotentialEnergy, Forces};

/// Gradient-descent based energy minization.
///
/// This algorithm is very rough, and will not converge in all the situations.
/// However it is easy to use, and simple enough to be implemented quickly.
pub struct GradientDescent {
    /// Dumping factor
    gamma: f64,
    /// Max force norm convergence criterium
    df2: f64,
    /// Energy convergence criterium
    dE: f64,
    /// Has the last minization converged?
    is_converged: bool
}

impl GradientDescent {
    /// Create a GradientDescent with sensible default values for energy and
    /// force convergence criteria. The force criterium is `1e-5 kJ/mol/Å^2`,
    /// and the energy criterium is `1e-5 kJ/mol/Å^2`.
    pub fn new() -> GradientDescent {
        let delta_f = units::from(1e-5, "kJ/mol/A").unwrap();
        let delta_E = units::from(1e-5, "kJ/mol").unwrap();
        GradientDescent::with_criteria(delta_f, delta_E)
    }

    /// Create a new `GradientDescent` with the force convergence criterium of
    /// `force`, and the energy convergence criterium of `energy`.
    pub fn with_criteria(force: f64, energy: f64) -> GradientDescent {
        let gamma = units::from(0.1, "fs^2/u").unwrap();
        GradientDescent{
            gamma: gamma,
            dE: energy,
            df2: force*force,
            is_converged: false
        }
    }

    /// Has the minization converged so far ?
    pub fn converged(&self) -> bool {
        self.is_converged
    }
}

impl Propagator for GradientDescent {
    fn setup(&mut self, _: &Universe) {
        self.is_converged = false;
    }

    fn propagate(&mut self, universe: &mut Universe) {
        let forces = Forces.compute(&universe);
        // Get the maximal value in the vector.
        // How can you know that fold(cte, cte) will do it ?
        let max_force_norm2 = forces.iter().map(|&f| f.norm2()).fold(-1./0. /* -inf */, f64::max);
        if max_force_norm2 < self.df2 {
            self.is_converged = true;
            return;
        }

        let energy = PotentialEnergy.compute(&universe);
        if energy < self.dE {
            self.is_converged = true;
            return;
        }

        // Store the current coordinates
        let mut positions = Vec::with_capacity(universe.size());
        for p in universe.iter() {
            positions.push(p.position().clone());
        }
        let positions = positions;

        let old_gamma = self.gamma;
        // Update coordinates, reducing gamma until we find a configuration of
        // lower energy
        loop {
            for (i, p) in universe.iter_mut().enumerate() {
                p.set_position(positions[i] + self.gamma * forces[i]);
            }
            let new_energy = PotentialEnergy.compute(&universe);
            if new_energy <= energy {
                break;
            }
            self.gamma /= 2.0;
        }

        // If we had a successful iteration without needing to reduce gamma,
        // we increase it slightly
        if self.gamma == old_gamma {
            self.gamma *= 1.1;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ::universe::*;
    use ::types::*;
    use ::potentials::*;
    use super::super::Propagator;

    #[test]
    fn minization() {
        let mut universe = Universe::from_cell(UnitCell::cubic(20.0));;
        universe.add_particle(Particle::new("Cl"));
        universe[0].set_position(Vector3D::new(0.0, 0.0, 0.0));
        universe.add_particle(Particle::new("Cl"));
        universe[1].set_position(Vector3D::new(0.0, 0.0, 2.0));
        universe.add_pair_interaction("Cl", "Cl", Harmonic{x0: 2.3, k: 0.1});

        let mut minization = GradientDescent::new();
        for _ in 0..100 {
            minization.propagate(&mut universe);
        }
        assert!(minization.converged());
        assert_approx_eq!(universe.distance(0, 1), 2.3, 1e-4);
    }
}
