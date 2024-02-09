// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license
use std::f64::consts::{PI, FRAC_2_SQRT_PI};

use rayon::prelude::*;

use crate::math::erfc;
use crate::consts::FOUR_PI_EPSILON_0;
use crate::PairRestriction;
use crate::utils::ThreadLocalVec;
use crate::Configuration;
use crate::{Matrix3, Vector3D};

use super::{CoulombicPotential, GlobalCache, GlobalPotential};

/// Wolf summation for coulombic interactions.
///
/// This is a fast, direct, pairwise summation for coulombic potential
/// [Wolf1999].
///
/// # Examples
///
/// ```
/// # use lumol_core::sys::{Particle, Molecule, UnitCell, System};
/// # use lumol_core::energy::Wolf;
/// # use lumol_core::types::Vector3D;
///
/// // A relatively large cutoff is needed for Wolf summation
/// let wolf = Wolf::new(12.0);
///
/// // Setup a system containing a NaCl pair
/// let mut system = System::with_cell(UnitCell::cubic(30.0));
///
/// let mut na = Particle::new("Na");
/// na.charge = 1.0;
/// na.position = Vector3D::new(0.0, 0.0, 0.0);
///
/// let mut cl = Particle::new("Cl");
/// cl.charge = -1.0;
/// cl.position = Vector3D::new(2.0, 0.0, 0.0);
///
/// system.add_molecule(Molecule::new(na));
/// system.add_molecule(Molecule::new(cl));
///
/// // Use Wolf summation for electrostatic interactions
/// system.set_coulomb_potential(Box::new(wolf));
///
/// assert_eq!(system.potential_energy(), -0.0729290269539354);
/// ```
///
/// [Wolf1999]: Wolf, D. et al. J. Chem. Phys. 110, 8254 (1999).
#[derive(Clone)]
pub struct Wolf {
    /// Splitting parameter
    alpha: f64,
    /// Cutoff radius in real space
    cutoff: f64,
    /// Energy constant for wolf computation
    energy_constant: f64,
    /// Force constant for wolf computation
    force_constant: f64,
    /// Restriction scheme
    restriction: PairRestriction,
}

impl Wolf {
    /// Create a new Wolf summation, using a real-space cutoff of `cutoff`.
    pub fn new(cutoff: f64) -> Wolf {
        assert!(cutoff > 0.0, "Got a negative cutoff in Wolf summation");
        let alpha = PI / cutoff;

        let alpha_cutoff = alpha * cutoff;
        let alpha_cutoff_2 = alpha_cutoff * alpha_cutoff;

        let energy_constant = erfc(alpha_cutoff) / cutoff;
        let force_constant = erfc(alpha_cutoff) / (cutoff * cutoff) + FRAC_2_SQRT_PI * alpha * f64::exp(-alpha_cutoff_2) / cutoff;
        Wolf {
            alpha: alpha,
            cutoff: cutoff,
            energy_constant: energy_constant,
            force_constant: force_constant,
            restriction: PairRestriction::None,
        }
    }

    /// Compute the energy for the pair of particles with charge `qi` and `qj`,
    /// at the distance of `rij`. The `scaling` parameter comes from the
    /// restriction associated with this potential.
    #[inline]
    fn energy_pair(&self, qiqj: f64, rij: f64) -> f64 {
        if rij > self.cutoff {
            0.0
        } else {
            qiqj * (erfc(self.alpha * rij) / rij - self.energy_constant) / FOUR_PI_EPSILON_0
        }
    }

    /// Compute the energy for self interaction of a particle with charge `qi`
    #[inline]
    fn energy_self(&self, qi: f64) -> f64 {
        qi * qi * 0.5 * (self.energy_constant + self.alpha * FRAC_2_SQRT_PI) / FOUR_PI_EPSILON_0
    }

    /// Compute the force over the distance for the pair of particles with
    /// charge `qi` and `qj`, at the distance `rij`.
    #[inline]
    fn force_pair(&self, qiqj: f64, rij: f64) -> f64 {
        if rij > self.cutoff {
            0.0
        } else {
            let rij2 = rij * rij;
            let alpha_rij = self.alpha * rij;
            let exp_alpha_rij = f64::exp(-alpha_rij * alpha_rij);
            let factor = erfc(alpha_rij) / rij2 + self.alpha * FRAC_2_SQRT_PI * exp_alpha_rij / rij;
            return qiqj * (factor - self.force_constant) / (rij * FOUR_PI_EPSILON_0);
        }
    }
}

impl GlobalCache for Wolf {
    fn move_molecule_cost(
        &self,
        configuration: &Configuration,
        molecule_id: usize,
        new_positions: &[Vector3D],
    ) -> f64 {
        let mut old_energy = 0.0;
        let mut new_energy = 0.0;

        let charges = configuration.particles().charge;
        let positions = configuration.particles().position;

        // Iterate over all interactions between a particle in the moved
        // molecule and a particle in another molecule
        let molecule = configuration.molecule(molecule_id);
        for (i, part_i) in molecule.indexes().enumerate() {
            let qi = charges[part_i];
            if qi == 0.0 {
                continue;
            }

            for (_, other_molecule) in configuration.molecules().enumerate().filter(|(id, _)| molecule_id != *id) {
                for part_j in other_molecule.indexes() {
                    let qj = charges[part_j];
                    if qj == 0.0 {
                        continue;
                    }

                    let path = configuration.bond_path(part_i, part_j);
                    let info = self.restriction.information(path);
                    if info.excluded {
                        continue;
                    }

                    let old_r = configuration.distance(part_i, part_j);
                    let new_r = configuration.cell.distance(&new_positions[i], &positions[part_j]);

                    old_energy += info.scaling * self.energy_pair(qi * qj, old_r);
                    new_energy += info.scaling * self.energy_pair(qi * qj, new_r);
                }
            }
        }

        return new_energy - old_energy;
    }

    fn update(&self) {
        // Nothing to do
    }
}

impl GlobalPotential for Wolf {
    fn cutoff(&self) -> Option<f64> {
        Some(self.cutoff)
    }

    fn energy(&self, configuration: &Configuration) -> f64 {
        let natoms = configuration.size();
        let charges = configuration.particles().charge;

        let energies = (0..natoms).into_par_iter().map(|i| {
            let mut energy = 0.0;
            let qi = charges[i];
            if qi == 0.0 {
                return 0.0;
            }

            for j in i + 1..natoms {
                let qj = charges[j];
                if qj == 0.0 {
                    continue;
                }

                let path = configuration.bond_path(i, j);
                let info = self.restriction.information(path);
                if info.excluded {
                    continue;
                }

                let rij = configuration.distance(i, j);
                energy += info.scaling * self.energy_pair(qi * qj, rij);
            }

            return energy - self.energy_self(qi);
        });
        return energies.sum();
    }

    fn forces(&self, configuration: &Configuration, forces: &mut [Vector3D]) {
        assert_eq!(forces.len(), configuration.size());

        let natoms = configuration.size();
        let charges = configuration.particles().charge;
        // To avoid race conditions, each thread needs its own local forces Vec
        let thread_local_forces = ThreadLocalVec::with_size(natoms);

        (0..natoms).into_par_iter().for_each(|i| {
            // Get the thread local forces Vec
            let mut forces = thread_local_forces.borrow_mut();

            let mut force_i = Vector3D::zero();
            let qi = charges[i];
            if qi == 0.0 {
                return;
            }
            for j in i + 1..natoms {
                let qj = charges[j];
                if qj == 0.0 {
                    continue;
                }

                let path = configuration.bond_path(i, j);
                let info = self.restriction.information(path);
                if info.excluded {
                    continue;
                }

                let rij = configuration.nearest_image(i, j);
                let force = info.scaling * self.force_pair(qi * qj, rij.norm()) * rij;
                force_i += force;
                forces[j] -= force;
            }
            forces[i] += force_i;
        });

        // At this point all the forces are computed, but the results are
        // scattered across all thread local Vecs, here we gather them.
        thread_local_forces.sum_into(forces);
    }

    fn atomic_virial(&self, configuration: &Configuration) -> Matrix3 {
        let natoms = configuration.size();
        let charges = configuration.particles().charge;

        let virials = (0..natoms).into_par_iter().map(|i| {
            let qi = charges[i];
            if qi == 0.0 {
                return Matrix3::zero();
            }
            let mut local_virial = Matrix3::zero();

            for j in i + 1..natoms {
                let qj = charges[j];
                if qj == 0.0 {
                    continue;
                }

                let path = configuration.bond_path(i, j);
                let info = self.restriction.information(path);
                if info.excluded {
                    continue;
                }

                let rij = configuration.nearest_image(i, j);
                let force = info.scaling * self.force_pair(qi * qj, rij.norm()) * rij;
                local_virial += force.tensorial(&rij);
            }

            local_virial
        });

        return virials.sum();
    }

    fn molecular_virial(&self, configuration: &Configuration) -> Matrix3 {
        let charges = configuration.particles().charge;
        let virials = configuration.molecules().enumerate().par_bridge().map(|(i, molecule_i)| {
            let mut local_virial = Matrix3::zero();
            let ri = molecule_i.center_of_mass();

            for molecule_j in configuration.molecules().skip(i + 1) {
                let rj = molecule_j.center_of_mass();
                let mut r_ij = ri - rj;
                configuration.cell.vector_image(&mut r_ij);

                for part_a in molecule_i.indexes() {
                    let q_a = charges[part_a];
                    if q_a == 0.0 {
                        continue;
                    }

                    for part_b in molecule_j.indexes() {
                        let q_b = charges[part_b];
                        if q_b == 0.0 {
                            continue;
                        }

                        let path = configuration.bond_path(part_a, part_b);
                        let info = self.restriction.information(path);
                        if info.excluded {
                            continue;
                        }

                        let r_ab = configuration.nearest_image(part_a, part_b);
                        let force = info.scaling * self.force_pair(q_a * q_b, r_ab.norm()) * r_ab;
                        let w_ab = force.tensorial(&r_ab);
                        local_virial += w_ab * (r_ab * r_ij) / r_ab.norm2();
                     }
                 }
             }
             return local_virial;
         });
         return virials.sum();
     }
}

impl CoulombicPotential for Wolf {
    fn set_restriction(&mut self, restriction: PairRestriction) {
        self.restriction = restriction;
    }
}

#[cfg(test)]
mod tests {
    pub use super::*;
    use crate::{System, Matrix3};
    use crate::GlobalPotential;
    use crate::utils::system_from_xyz;

    use approx::{assert_ulps_eq, assert_relative_eq};

    pub fn testing_system() -> System {
        let mut system = system_from_xyz(
            "2
            cell: 20.0
            Cl 0.0 0.0 0.0
            Na 1.5 0.0 0.0
            ",
        );
        system.particles_mut().charge[0] = -1.0;
        system.particles_mut().charge[1] = 1.0;
        return system;
    }

    #[test]
    #[allow(clippy::unreadable_literal)]
    fn energy() {
        const E_BRUTE_FORCE: f64 = -0.09262397663346732;

        let system = testing_system();
        let wolf = Wolf::new(8.0);

        let e = wolf.energy(&system);
        // Wolf is not very good for heterogeneous systems
        assert_ulps_eq!(e, E_BRUTE_FORCE, epsilon = 1e-2);
    }

    #[test]
    fn forces() {
        let mut system = testing_system();
        let wolf = Wolf::new(8.0);

        let mut forces = vec![Vector3D::zero(); system.size()];
        wolf.forces(&system, &mut forces);
        let norm = (forces[0] + forces[1]).norm();
        // Total force should be null
        assert_ulps_eq!(norm, 0.0);

        // Finite difference computation of the force
        let e = wolf.energy(&system);
        let eps = 1e-9;
        system.particles_mut().position[0][0] += eps;

        let e1 = wolf.energy(&system);
        let mut forces = vec![Vector3D::zero(); system.size()];
        wolf.forces(&system, &mut forces);
        assert_relative_eq!((e - e1) / eps, forces[0][0], epsilon = 1e-6);
    }

    #[test]
    fn atomic_virial() {
        let system = testing_system();
        let wolf = Wolf::new(8.0);

        let mut forces = vec![Vector3D::zero(); system.size()];
        wolf.forces(&system, &mut forces);
        let force = forces[0][0];
        let expected = Matrix3::new([[-force * 1.5, 0.0, 0.0], [0.0; 3], [0.0; 3]]);

        assert_eq!(wolf.atomic_virial(&system), expected);
    }

    #[test]
    fn atomic_virial_finite_differences() {
        fn scale(system: &mut System, i: usize, j: usize, eps: f64) {
            let mut scaling = Matrix3::one();
            scaling[i][j] += eps;
            let old_cell = system.cell;
            let new_cell = system.cell.scale(scaling);

            for position in system.particles_mut().position {
                *position = new_cell.cartesian(&old_cell.fractional(position));
            }
            system.cell = new_cell;
        }

        let eps = 1e-9;
        let mut system = testing_system();
        let wolf = Wolf::new(8.0);

        let virial = wolf.atomic_virial(&system);

        let mut finite_diff = Matrix3::zero();
        let energy_0 = wolf.energy(&system);

        scale(&mut system, 0, 0, eps);
        let energy_1 = wolf.energy(&system);
        finite_diff[0][0] = (energy_0 - energy_1) / eps;

        scale(&mut system, 1, 0, eps);
        let energy_2 = wolf.energy(&system);
        finite_diff[1][0] = (energy_1 - energy_2) / eps;

        scale(&mut system, 2, 0, eps);
        let energy_3 = wolf.energy(&system);
        finite_diff[2][0] = (energy_2 - energy_3) / eps;

        assert_relative_eq!(virial, finite_diff, epsilon = 1e-5);
    }

    mod cache {
        use super::*;
        use crate::{CoulombicPotential, GlobalCache, GlobalPotential, PairRestriction};
        use crate::System;
        use crate::Vector3D;
        use crate::utils::system_from_xyz;

        pub fn testing_system() -> System {
            let mut system = system_from_xyz(
                "6
                cell: 20.0
                O  0.0  0.0  0.0
                H -0.7 -0.7  0.3
                H  0.3 -0.3 -0.8
                O  2.0  2.0  0.0
                H  1.3  1.3  0.3
                H  2.3  1.7 -0.8
                ",
            );
            assert!(system.add_bond(0, 1).is_empty());
            assert!(system.add_bond(0, 2).is_empty());
            assert!(system.add_bond(3, 4).is_empty());
            assert!(system.add_bond(3, 5).is_empty());
            assert!(system.molecules().count() == 2);

            for particle in system.particles_mut() {
                if particle.name == "O" {
                    *particle.charge = -0.8476;
                } else if particle.name == "H" {
                    *particle.charge = 0.4238;
                }
            }
            return system;
        }

        #[test]
        #[allow(clippy::unreadable_literal)]
        fn move_rigid_molecule() {
            let mut system = testing_system();
            let mut wolf = Wolf::new(8.0);
            wolf.set_restriction(PairRestriction::InterMolecular);

            let check = wolf.clone();

            let old_energy = check.energy(&system);

            let new_positions = &[
                Vector3D::new(4.0, 0.0, -2.0),
                Vector3D::new(3.010010191494968, 0.19045656166589708, -2.1166435218719863),
                Vector3D::new(4.0761078062722484, -0.8995901989882638, -2.0703212322750546),
            ];
            let cost = wolf.move_molecule_cost(&system, 0, new_positions);

            system.particles_mut().position[0] = new_positions[0];
            system.particles_mut().position[1] = new_positions[1];
            system.particles_mut().position[2] = new_positions[2];
            let new_energy = check.energy(&system);
            assert_ulps_eq!(cost, new_energy - old_energy);
        }
    }
}
