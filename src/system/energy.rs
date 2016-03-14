// Cymbalum, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

//! Primitives for energy computation.
//!
//! This module provides simple function to compute separated components of the
//! potential energy of an `System`.
use system::System;

/// An helper struct to evaluate energy components of a system.
pub struct EnergyEvaluator<'a> {
    system: &'a System,
}

impl<'a> EnergyEvaluator<'a> {
    /// Create a new energy evaluator acting on the given `system`.
    pub fn new(system: &'a System) -> EnergyEvaluator<'a> {
        EnergyEvaluator {
            system: system,
        }
    }
}

impl<'a> EnergyEvaluator<'a> {
    /// Compute the energy associated with the pair of particles `i, j` at
    /// distance `r`
    #[inline]
    pub fn pair(&self, r: f64, i: usize, j: usize) -> f64 {
        let mut energy = 0.0;
        for &(ref potential, ref restriction) in self.system.pair_potentials(i, j) {
            let info = restriction.informations(self.system, i, j);
            if !info.excluded {
                energy += info.scaling * potential.energy(r);
            }
        }
        return energy;
    }

    /// Compute the energy of all the pairs in the system
    pub fn pairs(&self) -> f64 {
        let mut energy = 0.0;
        for i in 0..self.system.size() {
            for j in (i+1)..self.system.size() {
                let r = self.system.wraped_vector(i, j).norm();
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
        for potential in self.system.bond_potentials(i, j) {
            energy += potential.energy(r);
        }
        return energy;
    }

    /// Compute the energy of all the bonds in the system
    pub fn bonds(&self) -> f64 {
        let mut energy = 0.0;
        for molecule in self.system.molecules() {
            for bond in molecule.bonds() {
                let (i, j) = (bond.i(), bond.j());
                let r = self.system.wraped_vector(i, j).norm();
                energy += self.bond(r, i, j);
            }
        }
        return energy;
    }

    /// Compute the energy associated with the angle `i, j, k` at angle `theta`
    #[inline]
    pub fn angle(&self, theta: f64, i: usize, j: usize, k: usize) -> f64 {
        let mut energy = 0.0;
        for potential in self.system.angle_potentials(i, j, k) {
            energy += potential.energy(theta);
        }
        return energy;
    }

    /// Compute the energy of all the angles in the system
    pub fn angles(&self) -> f64 {
        let mut energy = 0.0;
        for molecule in self.system.molecules() {
            for angle in molecule.angles() {
                let (i, j, k) = (angle.i(), angle.j(), angle.k());
                let theta = self.system.angle(i, j, k);
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
        for potential in self.system.dihedral_potentials(i, j, k, m) {
            energy += potential.energy(phi);
        }
        return energy;
    }

    /// Compute the energy of all the dihedral angles in the system
    pub fn dihedrals(&self) -> f64 {
        let mut energy = 0.0;
        for molecule in self.system.molecules() {
            for dihedral in molecule.dihedrals() {
                let (i, j, k, m) = (dihedral.i(), dihedral.j(), dihedral.k(), dihedral.m());
                let phi = self.system.dihedral(i, j, k, m);
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
        if let Some(coulomb) = self.system.coulomb_potential() {
            coulomb.borrow_mut().energy(self.system)
        } else {
            0.0
        }
    }

    /// Compute the energy of the global potentials
    #[inline]
    pub fn global(&self) -> f64 {
        let mut energy = 0.0;
        for global in self.system.global_potentials() {
            energy += global.borrow_mut().energy(self.system);
        }
        return energy;
    }
}

#[cfg(test)]
mod tests {
    // This code is tested in the src/simulation/compute.rs file, as it is used
    // by the PotentialEnergy compute.
}
