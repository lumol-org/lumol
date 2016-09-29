// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

//! Primitives for energy computation.
//!
//! This module provides simple function to compute separated components of the
//! potential energy of an `System`.
use system::System;
use potentials::Potential;

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
        for potential in self.system.pair_potentials(i, j) {
            let info = potential.restriction().informations(self.system, i, j);
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
                let r = self.system.nearest_image(i, j).norm();
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
                let r = self.system.nearest_image(i, j).norm();
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

    /// Compute the energy of the electrostatic interactions
    #[inline]
    pub fn coulomb(&self) -> f64 {
        if let Some(coulomb) = self.system.interactions().coulomb() {
            coulomb.borrow_mut().energy(self.system)
        } else {
            0.0
        }
    }

    /// Compute the energy of the global potentials
    #[inline]
    pub fn global(&self) -> f64 {
        let mut energy = 0.0;
        for global in self.system.interactions().globals() {
            energy += global.borrow_mut().energy(self.system);
        }
        return energy;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use system::{System, Particle, UnitCell};
    use potentials::{Harmonic, LennardJones, PairInteraction};
    use types::Vector3D;
    use utils::unit_from;

    fn testing_system() -> System {
        let mut system = System::from_cell(UnitCell::cubic(10.0));;
        system.add_particle(Particle::new("F"));
        system[0].position = Vector3D::new(0.0, 0.0, 0.0);
        system.add_particle(Particle::new("F"));
        system[1].position = Vector3D::new(1.0, 0.0, 0.0);
        system.add_particle(Particle::new("F"));
        system[2].position = Vector3D::new(1.0, 1.0, 0.0);
        system.add_particle(Particle::new("F"));
        system[3].position = Vector3D::new(2.0, 1.0, 0.0);

        assert!(system.add_bond(0, 1).is_empty());
        assert!(system.add_bond(1, 2).is_empty());
        assert!(system.add_bond(2, 3).is_empty());

        let lj = Box::new(LennardJones {
            epsilon: unit_from(100.0, "kJ/mol/A^2"),
            sigma: unit_from(0.8, "A")
        });
        system.interactions_mut().add_pair("F", "F", PairInteraction::new(lj, 5.0));

        system.interactions_mut().add_bond("F", "F",
            Box::new(Harmonic{
                k: unit_from(100.0, "kJ/mol/A^2"),
                x0: unit_from(2.0, "A")
        }));

        system.interactions_mut().add_angle("F", "F", "F",
            Box::new(Harmonic{
                k: unit_from(100.0, "kJ/mol/deg^2"),
                x0: unit_from(88.0, "deg")
        }));

        system.interactions_mut().add_dihedral("F", "F", "F", "F",
            Box::new(Harmonic{
                k: unit_from(100.0, "kJ/mol/deg^2"),
                x0: unit_from(185.0, "deg")
        }));

        return system;
    }

    #[test]
    fn pairs() {
        let system = testing_system();
        let evaluator = EnergyEvaluator::new(&system);
        assert_approx_eq!(evaluator.pairs(), unit_from(-258.3019360389957, "kJ/mol"));
    }

    #[test]
    fn bonds() {
        let system = testing_system();
        let evaluator = EnergyEvaluator::new(&system);
        assert_approx_eq!(evaluator.bonds(), unit_from(150.0, "kJ/mol"));
    }

    #[test]
    fn angles() {
        let system = testing_system();
        let evaluator = EnergyEvaluator::new(&system);
        assert_approx_eq!(evaluator.angles(), unit_from(400.0, "kJ/mol"));
    }

    #[test]
    fn dihedrals() {
        let system = testing_system();
        let evaluator = EnergyEvaluator::new(&system);
        assert_approx_eq!(evaluator.dihedrals(), unit_from(1250.0, "kJ/mol"));
    }
}
