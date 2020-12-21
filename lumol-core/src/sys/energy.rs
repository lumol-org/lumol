// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

//! Primitives for energy computation.
//!
//! This module provides simple function to compute separated components of the
//! potential energy of an `System`.

use std::f64::consts::PI;

use rayon::prelude::*;

use crate::BondPath;
use crate::System;

/// An helper struct to evaluate energy components of a system.
pub struct EnergyEvaluator<'a> {
    system: &'a System,
}

impl<'a> EnergyEvaluator<'a> {
    /// Create a new `EnergyEvaluator` acting on the given `system`.
    pub fn new(system: &'a System) -> EnergyEvaluator<'a> {
        EnergyEvaluator { system: system }
    }
}

impl<'a> EnergyEvaluator<'a> {
    /// Compute the energy associated with the pair of particles `i, j` at
    /// distance `r`
    #[inline]
    pub fn pair(&self, path: BondPath, r: f64, i: usize, j: usize) -> f64 {
        match self.system.pair_potential(i, j) {
            Some(potential) => {
                let info = potential.restriction().information(path);
                if !info.excluded {
                    info.scaling * potential.energy(r)
                } else {
                    0.0
                }
            }
            None => 0.0
        }
    }

    /// Compute the energy of all the pairs in the system
    pub fn pairs(&self) -> f64 {
        let energies = (0..self.system.size()).into_par_iter().map(|i| {
            let mut local_energy = 0.0;

            for j in (i + 1)..self.system.size() {
                let r = self.system.nearest_image(i, j).norm();
                let path = self.system.bond_path(i, j);
                local_energy += self.pair(path, r, i, j);
            }
            local_energy
        });
        return energies.sum();
    }

    /// Compute the energy due to long range corrections for the pairs
    #[inline]
    pub fn pairs_tail(&self) -> f64 {
        if self.system.cell.is_infinite() {
            return 0.0;
        }
        let mut energy = 0.0;
        let volume = self.system.volume();
        let composition = self.system.composition();
        for (i, ni) in composition.all_particles() {
            for (j, nj) in composition.all_particles() {
                let two_pi_density = 2.0 * PI * (ni as f64) * (nj as f64) / volume;
                if let Some(potential) = self.system.interactions().pair((i, j)) {
                    energy += two_pi_density * potential.tail_energy();
                }
            }
        }
        return energy;
    }

    /// Compute the energy associated with the bonded particles `i, j` at
    /// distance `r`
    #[inline]
    pub fn bond(&self, r: f64, i: usize, j: usize) -> f64 {
        self.system.bond_potential(i, j)
                   .map_or(0.0, |potential| potential.energy(r))
    }

    /// Compute the energy of all the bonds in the system
    pub fn bonds(&self) -> f64 {
        let mut energy = 0.0;
        for molecule in self.system.molecules() {
            for bond in molecule.bonds() {
                let (i, j) = (bond.i(), bond.j());
                let r = self.system.nearest_image(i, j).norm();
                energy += self.bond(r, i, j);
            }
        }
        return energy;
    }

    /// Compute the energy associated with the angle `i, j, k` at angle `theta`
    #[inline]
    pub fn angle(&self, theta: f64, i: usize, j: usize, k: usize) -> f64 {
        self.system.angle_potential(i, j, k)
                   .map_or(0.0, |potential| potential.energy(theta))
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
        self.system.dihedral_potential(i, j, k, m)
                   .map_or(0.0, |potential| potential.energy(phi))
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

    /// Compute the energy of the electrostatic interactions
    #[inline]
    pub fn coulomb(&self) -> f64 {
        self.system.coulomb_potential()
            .map_or(0.0, |coulomb| coulomb.energy(self.system))
    }

    /// Compute the energy of the global potentials
    #[inline]
    pub fn global(&self) -> f64 {
        let mut energy = 0.0;
        for global in self.system.global_potentials() {
            energy += global.energy(self.system);
        }
        return energy;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Harmonic, LennardJones, NullPotential, PairInteraction};
    use crate::{System, UnitCell};
    use crate::utils::system_from_xyz;
    use crate::units;

    use approx::assert_ulps_eq;

    fn testing_system() -> System {
        let mut system = system_from_xyz(
            "4
            cell: 10.0
            F 0.0 0.0 0.0
            F 1.0 0.0 0.0
            F 1.0 1.0 0.0
            F 2.0 1.0 0.0
            ",
        );
        assert!(system.add_bond(0, 1).is_empty());
        assert!(system.add_bond(1, 2).is_empty());
        assert!(system.add_bond(2, 3).is_empty());

        let mut pair = PairInteraction::new(
            Box::new(LennardJones {
                epsilon: units::from(100.0, "kJ/mol/A^2").unwrap(),
                sigma: units::from(0.8, "A").unwrap(),
            }),
            5.0,
        );
        pair.enable_tail_corrections();

        system.set_pair_potential(("F", "F"), pair);

        system.set_bond_potential(
            ("F", "F"),
            Box::new(Harmonic {
                k: units::from(100.0, "kJ/mol/A^2").unwrap(),
                x0: units::from(2.0, "A").unwrap(),
            }),
        );

        system.set_angle_potential(
            ("F", "F", "F"),
            Box::new(Harmonic {
                k: units::from(100.0, "kJ/mol/deg^2").unwrap(),
                x0: units::from(88.0, "deg").unwrap(),
            }),
        );

        system.set_dihedral_potential(
            ("F", "F", "F", "F"),
            Box::new(Harmonic {
                k: units::from(100.0, "kJ/mol/deg^2").unwrap(),
                x0: units::from(185.0, "deg").unwrap(),
            }),
        );

        // unused interaction to check that we do handle this right
        system.set_pair_potential(("H", "O"), PairInteraction::new(Box::new(NullPotential), 0.0));

        return system;
    }

    #[test]
    #[allow(clippy::unreadable_literal)]
    fn pairs() {
        let system = testing_system();
        let evaluator = EnergyEvaluator::new(&system);
        assert_ulps_eq!(evaluator.pairs(), units::from(-258.3019360389957, "kJ/mol").unwrap());
        assert_ulps_eq!(evaluator.pairs_tail(), -0.0000028110338032153973);
    }

    #[test]
    fn pairs_tail_infinite_cell() {
        let mut system = testing_system();
        system.cell = UnitCell::infinite();

        let evaluator = EnergyEvaluator::new(&system);
        assert_eq!(evaluator.pairs_tail(), 0.0);
    }

    #[test]
    fn bonds() {
        let system = testing_system();
        let evaluator = EnergyEvaluator::new(&system);
        assert_ulps_eq!(evaluator.bonds(), units::from(150.0, "kJ/mol").unwrap());
    }

    #[test]
    fn angles() {
        let system = testing_system();
        let evaluator = EnergyEvaluator::new(&system);
        assert_ulps_eq!(evaluator.angles(), units::from(400.0, "kJ/mol").unwrap());
    }

    #[test]
    fn dihedrals() {
        let system = testing_system();
        let evaluator = EnergyEvaluator::new(&system);
        assert_ulps_eq!(evaluator.dihedrals(), units::from(1250.0, "kJ/mol").unwrap(), max_ulps = 15);
    }
}
