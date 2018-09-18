// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license
use std::f64::consts::PI;

use consts::ELCC;
use energy::{PairRestriction, RestrictionInfo};
use math::*;
use parallel::ThreadLocalStore;
use parallel::prelude::*;
use sys::Configuration;
use types::{Matrix3, Vector3D, Zero};

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
/// assert_eq!(system.potential_energy(), -0.07292902695393541);
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
    energy_cst: f64,
    /// Force constant for wolf computation
    force_cst: f64,
    /// Restriction scheme
    restriction: PairRestriction,
}

impl Wolf {
    /// Create a new Wolf summation, using a real-space cutoff of `cutoff`.
    pub fn new(cutoff: f64) -> Wolf {
        assert!(cutoff > 0.0, "Got a negative cutoff in Wolf summation");
        let alpha = PI / cutoff;
        let e_cst = erfc(alpha * cutoff) / cutoff;
        let f_cst = erfc(alpha * cutoff) / (cutoff * cutoff)
            + 2.0 * alpha / sqrt(PI) * exp(-alpha * alpha * cutoff * cutoff) / cutoff;
        Wolf {
            alpha: alpha,
            cutoff: cutoff,
            energy_cst: e_cst,
            force_cst: f_cst,
            restriction: PairRestriction::None,
        }
    }

    /// Compute the energy for the pair of particles with charge `qi` and `qj`,
    /// at the distance of `rij`. The `scaling` parameter comes from the
    /// restriction associated with this potential.
    #[inline]
    fn energy_pair(&self, info: RestrictionInfo, qiqj: f64, rij: f64) -> f64 {
        if rij > self.cutoff || info.excluded {
            return 0.0;
        }
        info.scaling * qiqj * (erfc(self.alpha * rij) / rij - self.energy_cst) / ELCC
    }

    /// Compute the energy for self interaction of a particle with charge `qi`
    #[inline]
    fn energy_self(&self, qi: f64) -> f64 {
        qi * qi * (self.energy_cst / 2.0 + self.alpha / sqrt(PI)) / ELCC
    }

    /// Compute the force for self the pair of particles with charge `qi` and
    /// `qj`, at the distance of `rij`. The `scaling` parameter comes from the
    /// restriction associated with this potential.
    #[inline]
    fn force_pair(&self, info: RestrictionInfo, qiqj: f64, rij: Vector3D) -> Vector3D {
        let d = rij.norm();
        if d > self.cutoff || info.excluded {
            return Vector3D::zero();
        }
        let factor = erfc(self.alpha * d) / (d * d)
            + 2.0 * self.alpha / sqrt(PI) * exp(-self.alpha * self.alpha * d * d) / d;
        return info.scaling * qiqj * (factor - self.force_cst) * rij.normalized() / ELCC;
    }
}

impl GlobalCache for Wolf {
    fn move_molecule_cost(
        &self,
        configuration: &Configuration,
        molecule_id: usize,
        new_positions: &[Vector3D],
    ) -> f64 {
        let mut old_energynergy = 0.0;
        let mut new_energynergy = 0.0;

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

                    let old_r = configuration.distance(part_i, part_j);
                    let new_r = configuration.cell.distance(&new_positions[i], &positions[part_j]);

                    let path = configuration.bond_path(part_i, part_j);
                    let info = self.restriction.information(path);

                    old_energynergy += self.energy_pair(info, qi * qj, old_r);
                    new_energynergy += self.energy_pair(info, qi * qj, new_r);
                }
            }
        }

        return new_energynergy - old_energynergy;
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

        let energies = (0..natoms).par_map(|i| {
            let mut local_energy = 0.0;
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

                let rij = configuration.distance(i, j);
                local_energy += self.energy_pair(info, qi * qj, rij);
            }

            local_energy - self.energy_self(qi)
        });
        return energies.sum();
    }

    fn forces(&self, configuration: &Configuration, forces: &mut [Vector3D]) {
        assert_eq!(forces.len(), configuration.size());

        // To avoid race conditions, each thread needs its
        // own local forces Vec
        let natoms = configuration.size();
        let charges = configuration.particles().charge;
        let thread_forces_store = ThreadLocalStore::new(|| vec![Vector3D::zero(); natoms]);

        (0..natoms).into_par_iter().for_each(|i| {
            // Get the thread local forces Vec
            let mut thread_forces = thread_forces_store.borrow_mut();

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

                let rij = configuration.nearest_image(i, j);
                let force = self.force_pair(info, qi * qj, rij);
                thread_forces[i] += force;
                thread_forces[j] -= force;
            }
        });

        // At this point all the forces are computed, but the
        // results are scattered across all thread local Vecs,
        // here we gather them.
        thread_forces_store.sum_local_values(forces);
    }

    fn atomic_virial(&self, configuration: &Configuration) -> Matrix3 {
        let natoms = configuration.size();
        let charges = configuration.particles().charge;

        let virials = (0..natoms).par_map(|i| {
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

                let rij = configuration.nearest_image(i, j);
                let force = self.force_pair(info, qi * qj, rij);
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

                        let r_ab = configuration.nearest_image(part_a, part_b);
                        let force = self.force_pair(info, q_a * q_b, r_ab);
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
    use energy::GlobalPotential;
    use sys::System;
    use types::{Matrix3, One};
    use utils::system_from_xyz;

    const E_BRUTE_FORCE: f64 = -0.09262397663346732;

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
    fn energy() {
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
            let old_cell = system.cell.clone();
            let new_cell = system.cell.scale(scaling);

            for position in system.particles_mut().position {
                *position = new_cell.cartesian(&old_cell.fractional(&position));
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
        use energy::{CoulombicPotential, GlobalCache, GlobalPotential, PairRestriction};
        use sys::System;
        use types::Vector3D;

        pub fn testing_system() -> System {
            use utils::system_from_xyz;
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
        fn move_rigid_molecule() {
            let mut system = testing_system();
            let mut wolf = Wolf::new(8.0);
            wolf.set_restriction(PairRestriction::InterMolecular);

            let check = wolf.clone();

            let old_energy = check.energy(&system);

            let new_positions = &[
                Vector3D::new(4.0, 0.0, -2.0),
                Vector3D::new(3.0100101914949688, 0.19045656166589708, -2.11664352187198654),
                Vector3D::new(4.07610780627224825, -0.8995901989882638, -2.07032123227505443),
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
