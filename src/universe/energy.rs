/* Cymbalum, Molecular Simulation in Rust - Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
 */
//! Primitives for energy computation.
//!
//! This module provides simple function to compute separated components of the
//! potential energy of an `Universe`.
use universe::Universe;

/// An helper struct to evaluate energy components of an universe.
pub struct EnergyEvaluator<'a> {
    universe: &'a Universe,
}

impl<'a> EnergyEvaluator<'a> {
    /// Create a new energy evaluator acting on the given `universe`.
    pub fn new(universe: &'a Universe) -> EnergyEvaluator<'a> {
        EnergyEvaluator {
            universe: universe,
        }
    }
}

impl<'a> EnergyEvaluator<'a> {
    /// Compute the energy associated with the pair of particles `i, j` at
    /// distance `r`
    #[inline]
    pub fn pair(&self, r: f64, i: usize, j: usize) -> f64 {
        let mut energy = 0.0;
        for &(ref potential, ref restriction) in self.universe.pair_potentials(i, j) {
            if !restriction.is_excluded_pair(self.universe, i, j) {
                let s = restriction.scaling(self.universe, i, j);
                energy += s * potential.energy(r);
            }
        }
        return energy;
    }

    /// Compute the energy of all the pairs in the system
    pub fn pairs(&self) -> f64 {
        let mut energy = 0.0;
        for i in 0..self.universe.size() {
            for j in (i+1)..self.universe.size() {
                let r = self.universe.wraped_vector(i, j).norm();
                energy += self.pair(r, i, j);
            }
        }
        return energy;
    }

    /// Compute the energy associated with the bonded particles `i, j` at
    /// distance `r`
    #[inline]
    pub fn bond(&self, r: f64, i: usize, j: usize) -> f64 {
        let mut energy = 0.0;
        for potential in self.universe.bond_potentials(i, j) {
            energy += potential.energy(r);
        }
        return energy;
    }

    /// Compute the energy of all the bonds in the system
    pub fn bonds(&self) -> f64 {
        let mut energy = 0.0;
        for molecule in self.universe.molecules() {
            for bond in molecule.bonds() {
                let (i, j) = (bond.i(), bond.j());
                let r = self.universe.wraped_vector(i, j).norm();
                energy += self.bond(r, i, j);
            }
        }
        return energy;
    }

    /// Compute the energy associated with the angle `i, j, k` at angle `theta`
    #[inline]
    pub fn angle(&self, theta: f64, i: usize, j: usize, k: usize) -> f64 {
        let mut energy = 0.0;
        for potential in self.universe.angle_potentials(i, j, k) {
            energy += potential.energy(theta);
        }
        return energy;
    }

    /// Compute the energy of all the angles in the system
    pub fn angles(&self) -> f64 {
        let mut energy = 0.0;
        for molecule in self.universe.molecules() {
            for angle in molecule.angles() {
                let (i, j, k) = (angle.i(), angle.j(), angle.k());
                let theta = self.universe.angle(i, j, k);
                energy += self.angle(theta, i, j, k);
            }
        }
        return energy;
    }

    /// Compute the energy associated with the dihedral angle `i, j, k, m` at
    /// angle `phi`
    #[inline]
    pub fn dihedral(&self, phi: f64, i: usize, j: usize, k: usize, m: usize) -> f64 {
        let mut energy = 0.0;
        for potential in self.universe.dihedral_potentials(i, j, k, m) {
            energy += potential.energy(phi);
        }
        return energy;
    }

    /// Compute the energy of all the dihedral angles in the system
    pub fn dihedrals(&self) -> f64 {
        let mut energy = 0.0;
        for molecule in self.universe.molecules() {
            for dihedral in molecule.dihedrals() {
                let (i, j, k, m) = (dihedral.i(), dihedral.j(), dihedral.k(), dihedral.m());
                let phi = self.universe.dihedral(i, j, k, m);
                energy += self.dihedral(phi, i, j, k, m);
            }
        }
        return energy;
    }

    /// Compute the energy of the molecule interaction, *i.e.* the sum `bonds +
    /// angles + dihedrals`.
    #[inline]
    pub fn molecules(&self) -> f64 {
        return self.bonds() + self.angles() + self.dihedrals();
    }

    /// Compute the energy of the electrostatic interactions
    #[inline]
    pub fn coulomb(&self) -> f64 {
        if let Some(ref coulomb) = *self.universe.coulomb_potential() {
            coulomb.energy(self.universe)
        } else {
            0.0
        }
    }

    /// Compute the energy of the global potentials
    #[inline]
    pub fn global(&self) -> f64 {
        let mut energy = 0.0;
        for global in self.universe.global_potentials() {
            energy += global.energy(self.universe);
        }
        return energy;
    }
}

#[cfg(test)]
mod tests {
    // This code is tested in the src/simulation/compute.rs file, as it is used
    // by the PotentialEnergy compute.
}
