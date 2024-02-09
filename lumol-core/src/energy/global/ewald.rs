// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license
use std::ops::{Index, IndexMut, Deref, Range};
use std::sync::{RwLock, RwLockReadGuard, RwLockWriteGuard};
use std::f64::consts::{PI, FRAC_2_SQRT_PI};
use std::f64;

use rayon::prelude::*;

use soa_derive::StructOfArray;

use log::{warn, info};
use log_once::warn_once;

use crate::math::{erf, erfc};
use crate::{Configuration, UnitCell, CellShape};
use crate::{Matrix3, Vector3D, Array3, Complex};
use crate::consts::FOUR_PI_EPSILON_0;
use crate::{PairRestriction, RestrictionInfo};
use crate::utils::ThreadLocalVec;

use super::{GlobalPotential, CoulombicPotential, GlobalCache};

/// 3D array with negative indexing on the first dimension, for use in Ewald
/// phase factors.
///
/// # Examples
///
///  ```ignore
///  let array = Ewald3DArray::zeros((-6..5, 8, 2));
///
///  // Negative numbers are allowed for indexing into the array, as long as
///  // they fit in the range
///  array[(-4, 2, 1)] = Complex::polar(1.2, 0.0);
///  array[(3, 0, 0)] = Complex::polar(2.2, 0.0);
///  ```
#[derive(Clone, Debug)]
struct Ewald3DArray {
    data: Array3<Complex>,
    offset: isize,
}

impl Ewald3DArray {
    /// Create an array of the given size, and initialize it with zeros
    pub fn zeros((range, j, k): (Range<isize>, usize, usize)) -> Ewald3DArray {
        let offset = -range.start;
        let i = (range.end + offset) as usize;
        Ewald3DArray {
            data: Array3::zeros((i, j, k)),
            offset: offset
        }
    }

    /// Resize and fill the array with zeros if the new size does not match the
    /// old one.
    pub fn resize_if_different(&mut self, (range, j, k): (Range<isize>, usize, usize)) {
        self.offset = -range.start;
        let i = (range.end + self.offset) as usize;
        self.data.resize_if_different((i, j, k));
    }
}

impl Index<(isize, usize, usize)> for Ewald3DArray {
    type Output = Complex;
    fn index(&self, (i, j, k): (isize, usize, usize)) -> &Complex {
        let i = (i + self.offset) as usize;
        &self.data[(i, j, k)]
    }
}

impl IndexMut<(isize, usize, usize)> for Ewald3DArray {
    fn index_mut(&mut self, (i, j, k): (isize, usize, usize)) -> &mut Complex {
        let i = (i + self.offset) as usize;
        &mut self.data[(i, j, k)]
    }
}

/// Various parameters used by Ewald calculations.
///
/// They are grouped in a struct for easier passing as function arguments.
#[derive(Clone, Debug)]
pub struct EwaldParameters {
    /// Splitting parameter between k-space and real space
    pub alpha: f64,
    /// Cutoff radius in real space
    pub rc: f64,
    /// Number of points to use in k-space
    pub kmax: isize,
    /// Spherical cutoff in k-space
    pub kmax2: f64,
}

/// Various pre-factors used by Ewald computation
///
/// All of these factors only depend on the unit cell, and can be re-used if the
/// unit cell do not change.
///
/// Computing the factors account for the `\vec k = 0` and `k2 > kmax2` cases,
/// so iterating over the values in these vectors will give all the needed
/// k-points, and only them.
#[derive(Clone, Debug, StructOfArray)]
#[soa_derive(Clone)]
struct EwaldFactor {
    /// Index of the vector in k-space
    index: (isize, isize, isize),
    /// Energetic pre-factor: `4 π / V exp(- k² / (4 α²)) / k²`
    energy: f64,
    /// Electric field/force pre-factor: `8 π / V exp(- k² / (4 α²)) / k² \vec k / k`
    field: Vector3D,
    /// Virial pre-factor: `𝟙 - 2 (1 / k² + 1 / (4 α²)) \vec k ⊗ \vec k / k²`
    virial: Matrix3,
}

/// Compute the factors for the given `cell` and Ewald `parameters`
fn compute_ewald_factors(factors: &mut EwaldFactorVec, cell: &UnitCell, parameters: &EwaldParameters) {

    let kmax = parameters.kmax;
    let kmax3d = 4 * kmax * kmax * kmax + 6 * kmax * kmax + 3 * kmax;

    factors.clear();
    factors.reserve(kmax3d as usize);

    match cell.shape() {
        CellShape::Infinite => panic!("Ewald is not defined with infinite unit cell"),
        // TODO: there is a faster algorithm for orthorhombic cell
        CellShape::Orthorhombic | CellShape::Triclinic => compute_ewald_factors_triclinic(factors, cell, parameters),
    }
}

fn compute_ewald_factors_triclinic(factors: &mut EwaldFactorVec, cell: &UnitCell, parameters: &EwaldParameters) {
    let alpha_sq_inv_fourth = 0.25 / (parameters.alpha * parameters.alpha);
    let four_pi_v = 4.0 * PI / cell.volume();

    let factor_from_k_vector = |k_vector: Vector3D, k2: f64, index: (isize, isize, isize)| {
        let energy = four_pi_v * f64::exp(- k2 * alpha_sq_inv_fourth) / k2;
        let field = 2.0 * energy * k_vector;
        let virial_factor = -2.0 * (1.0 / k2 + alpha_sq_inv_fourth);
        let virial = Matrix3::one() + virial_factor * k_vector.tensorial(&k_vector);
        let virial = energy * virial;

        return EwaldFactor { index, energy, field, virial };
    };

    let kmax = parameters.kmax;
    // k-vectors with a positive `ikx`
    for ikx in 1..kmax {
        for iky in -kmax..kmax {
            for ikz in -kmax..kmax {
                let k_vector = cell.k_vector([ikx as f64, iky as f64, ikz as f64]);
                let k2 = k_vector.norm2();
                if k2 > parameters.kmax2 {
                    continue;
                }

                factors.push(factor_from_k_vector(k_vector, k2, (ikx, iky, ikz)));
            }
        }
    }

    // k-vectors with `ikx = 0`
    for iky in 1..kmax {
        for ikz in -kmax..kmax {
            let k_vector = cell.k_vector([0.0, iky as f64, ikz as f64]);
            let k2 = k_vector.norm2();
            if k2 > parameters.kmax2 {
                continue;
            }

            factors.push(factor_from_k_vector(k_vector, k2, (0, iky, ikz)));
        }
    }

    // k-vectors with `ikx = 0` and `ikz = 0`
    for ikz in 1..kmax {
        let k_vector = cell.k_vector([0.0, 0.0, ikz as f64]);
        let k2 = k_vector.norm2();
        if k2 > parameters.kmax2 {
            continue;
        }

        factors.push(factor_from_k_vector(k_vector, k2, (0, 0, ikz)));
    }
}

/// Ewald summation for coulombic interactions.
///
/// The Ewald summation is based on a separation of the coulombic potential `U`
/// in two parts, using the trivial identity:
///
/// $$ U(x) = U(x) \times (f(x) + 1) - U(x) \times f(x) $$
///
/// where `f` is the `erf` function. This leads to a separation of the
/// conditionally convergent coulombic sum into two absolutely convergent sums:
/// one in real space, and the other in Fourier or k-space. For more information
/// about this algorithm see [FS2002].
///
/// [FS2002] Frenkel, D. & Smith, B. Understanding molecular simulation. (Academic press, 2002).
///
/// # Examples
///
/// ```
/// # use lumol_core::sys::{Particle, Molecule, UnitCell, System};
/// # use lumol_core::energy::{Ewald, SharedEwald};
/// # use lumol_core::types::Vector3D;
/// let ewald = SharedEwald::new(
///     Ewald::new(/* cutoff */ 12.0, /* kmax */ 7, /* alpha */ None)
/// );
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
/// // Use Ewald summation for electrostatic interactions
/// system.set_coulomb_potential(Box::new(ewald));
///
/// println!("energy is {}", system.potential_energy());
/// ```
pub struct Ewald {
    /// Various Ewald parameters
    parameters: EwaldParameters,
    /// Ewald pre-factors, only depending on the system unit cell
    factors: EwaldFactorVec,
    /// Restriction scheme
    restriction: PairRestriction,
    /// Cached phase factors (e^{i k r})
    eikr: Ewald3DArray,
    /// Fourier transform of the electrostatic density (\sum q_i e^{i k r})
    ///
    /// The vector contain the terms corresponding to the k-vectors in
    /// `self.factors.index`
    rho: Vec<Complex>,
    /// Caching the allocation for electric field calculation
    ///
    /// This will contain the electric field at each atom
    field: Vec<Vector3D>,
    /// Guard for cache invalidation of `self.factors`
    previous_cell: Option<UnitCell>,
    /// Update the cached quantities
    updater: Option<Box<dyn Fn(&mut Ewald) + Sync + Send>>,
}

impl Clone for Ewald {
    fn clone(&self) -> Ewald {
        Ewald {
            parameters: self.parameters.clone(),
            factors: self.factors.clone(),
            restriction: self.restriction,
            eikr: self.eikr.clone(),
            rho: self.rho.clone(),
            field: self.field.clone(),
            previous_cell: self.previous_cell,
            updater: None,
        }
    }
}

impl Deref for Ewald {
    type Target = EwaldParameters;
    fn deref(&self) -> &EwaldParameters {
        &self.parameters
    }
}

impl Ewald {
    /// Create an Ewald summation using the given `cutoff` radius in real space,
    /// and `kmax` points in k-space (Fourier space). If `alpha` is None, then
    /// the default value of `π / cutoff` is used.
    pub fn new<I: Into<Option<f64>>>(cutoff: f64, kmax: usize, alpha: I) -> Ewald {
        let alpha = alpha.into().unwrap_or(PI / cutoff);
        if cutoff < 0.0 {
            panic!("the cutoff can not be negative in Ewald");
        } else if alpha < 0.0 {
            panic!("alpha can not be negative in Ewald");
        } else if kmax == 0 {
            panic!("kmax can not be 0 in Ewald");
        }

        let parameters = EwaldParameters {
            alpha: alpha,
            rc: cutoff,
            kmax: kmax as isize,
            kmax2: 0.0,
        };

        Ewald {
            parameters: parameters,
            restriction: PairRestriction::None,
            factors: EwaldFactorVec::new(),
            eikr: Ewald3DArray::zeros((0..0, 0, 0)),
            rho: Vec::new(),
            field: Vec::new(),
            previous_cell: None,
            updater: None,
        }
    }

    /// Create an Ewald solver with the given real space `cutoff`, setting
    /// `alpha` and `kmax` to ensure that the energy is computed with the
    /// specified relative `accuracy`. The optimal parameter depends on the
    /// exact `configuration` used: both the total number of charges, and the
    /// unit cell.
    pub fn with_accuracy(cutoff: f64, accuracy: f64, configuration: &Configuration) -> Ewald {
        if cutoff < 0.0 {
            panic!("the cutoff can not be negative in Ewald");
        } else if accuracy < 0.0 {
            panic!("accuracy can not be negative in Ewald");
        } else if accuracy > 1.0 {
            warn!("accuracy is bigger than 1 in Ewald::with_precision");
        }

        // Compute squared total charge
        let mut q2 = 0.0;
        for charge in configuration.particles().charge {
            q2 += charge * charge;
        }
        q2 /= FOUR_PI_EPSILON_0;

        let natoms = configuration.size() as f64;
        let lengths = configuration.cell.lengths();
        let alpha = accuracy * f64::sqrt(natoms * cutoff * lengths[0] * lengths[1] * lengths[2]) / (2.0 * q2);
        let alpha = if alpha >= 1.0 {
            (1.35 - 0.15 * f64::ln(accuracy)) / cutoff
        } else {
            f64::sqrt(-f64::ln(alpha)) / cutoff
        };

        let min_length = lengths.min();
        let error = |kmax| {
            let arg: f64 = PI * kmax / (alpha * min_length);
            FRAC_2_SQRT_PI * q2 * alpha / min_length / f64::sqrt(kmax * natoms) * f64::exp(-arg * arg)
        };

        let mut kmax = 1;
        while error(kmax as f64) > accuracy {
            kmax += 1;
        }

        info!("Setting Ewald summation parameters: cutoff = {}, alpha = {}, kmax = {}", cutoff, alpha, kmax);

        Ewald::new(cutoff, kmax, alpha)
    }

    fn prepare(&mut self, cell: &UnitCell) {
        if let Some(ref prev_cell) = self.previous_cell {
            if cell == prev_cell {
                // Do not recompute
                return;
            }
        }
        self.previous_cell = Some(*cell);

        let max = cell.k_vector([1.0, 1.0, 1.0]).max() * self.parameters.kmax as f64;
        self.parameters.kmax2 = 1.0001 * max * max;

        let half_min_length = cell.lengths().min() / 2.0;
        if self.parameters.rc > half_min_length {
            warn_once!("The Ewald cutoff is too high for this unit cell, energy and forces might be wrong.");
        }
        if f64::exp(- self.parameters.alpha * half_min_length) > 0.05 {
            warn_once!(
"Ewald alpha parameter is too low for this unit cell, energy and forces might be wrong.
You can manually set alpha to a slightly higher value (current alpha is {})",
                self.parameters.alpha
            );
        }

        compute_ewald_factors(&mut self.factors, cell, &self.parameters);
    }
}

/// Real space part of the summation
impl Ewald {
    /// Get the real-space energy for one pair at distance `r` with charges `qi`
    /// and `qj` ; and with restriction information for this pair in `info`.
    #[allow(clippy::float_cmp)]  // checking info.scaling
    #[inline]
    fn real_space_energy_pair(&self, info: RestrictionInfo, qiqj: f64, r: f64) -> f64 {
        assert_eq!(info.scaling, 1.0, "Scaling restriction scheme using Ewald are not implemented");
        debug_assert!(!(r > self.rc && info.excluded), "excluded atoms are too far apart");
        if r > self.rc {
            return 0.0;
        }

        if info.excluded {
            // use a correction for excluded interaction, removing the energy
            // from k-space
            - qiqj / FOUR_PI_EPSILON_0 * erf(self.alpha * r) / r
        } else {
            qiqj / FOUR_PI_EPSILON_0 * erfc(self.alpha * r) / r
        }
    }

    /// Get the real-space force for one pair at distance `r` with charges
    /// `qi` and `qj` ; and with restriction information for this pair in
    /// `info`.
    #[allow(clippy::float_cmp)]  // checking info.scaling
    #[inline]
    fn real_space_force_pair(&self, info: RestrictionInfo, qiqj: f64, r: f64) -> f64 {
        assert_eq!(info.scaling, 1.0, "Scaling restriction scheme using Ewald are not implemented");
        debug_assert!(!(r > self.rc && info.excluded), "excluded atoms are too far apart");
        if r > self.rc {
            return 0.0;
        }

        if info.excluded {
            // use a correction for excluded interaction, removing the force
            // from k-space
            qiqj / (FOUR_PI_EPSILON_0 * r * r) * (
                self.alpha * FRAC_2_SQRT_PI * f64::exp(-self.alpha * self.alpha * r * r)
                - erf(self.alpha * r) / r
            )
        } else {
            qiqj / (FOUR_PI_EPSILON_0 * r * r) * (
                self.alpha * FRAC_2_SQRT_PI * f64::exp(-self.alpha * self.alpha * r * r)
                + erfc(self.alpha * r) / r
            )
        }
    }

    /// Real space contribution to the energy
    fn real_space_energy(&self, configuration: &Configuration) -> f64 {
        let natoms = configuration.size();
        let charges = configuration.particles().charge;

        let energies = (0..natoms).into_par_iter().map(|i| {
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

                let r = configuration.distance(i, j);
                local_energy += self.real_space_energy_pair(info, qi * qj, r);
            }

            local_energy
        });
        return energies.sum();
    }

    /// Real space contribution to the forces
    fn real_space_forces(&self, configuration: &Configuration, forces: &mut [Vector3D]) {
        assert_eq!(forces.len(), configuration.size());

        let natoms = configuration.size();
        let charges = configuration.particles().charge;
        // Each thread (and not each iteration of the loop below) get its own
        // storage in a `ThreadLocalVec`.
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

                let rij = configuration.nearest_image(i, j);
                let force = self.real_space_force_pair(info, qi * qj, rij.norm()) * rij;
                force_i += force;
                forces[j] -= force;
            }
            forces[i] += force_i;
        });

        // At this point all the forces are computed, but the results are
        // scattered across all thread local Vec, here we gather them.
        thread_local_forces.sum_into(forces);
    }

    /// Real space contribution to the atomic virial
    fn real_space_atomic_virial(&self, configuration: &Configuration) -> Matrix3 {
        let natoms = configuration.size();
        let charges = configuration.particles().charge;

        let virial = (0..natoms).into_par_iter().map(|i| {
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
                let force = self.real_space_force_pair(info, qi * qj, rij.norm()) * rij;
                local_virial += force.tensorial(&rij);
            }
            local_virial
        });
        return virial.sum();
    }

    /// Real space contribution to the molecular virial
    fn real_space_molecular_virial(&self, configuration: &Configuration) -> Matrix3 {
        let charges = configuration.particles().charge;
        let virial = configuration.molecules().enumerate().par_bridge().map(|(i, molecule_i)| {
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
                        let force = self.real_space_force_pair(info, q_a * q_b, r_ab.norm()) * r_ab;
                        let w_ab = force.tensorial(&r_ab);
                        local_virial += w_ab * (r_ab * r_ij) / r_ab.norm2();
                    }
                }
            }
            return local_virial;
        });

        return virial.sum();
     }

     fn real_space_move_molecule_cost(
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

                    let old_r = configuration.distance(part_i, part_j);
                    let new_r = configuration.cell.distance(&new_positions[i], &positions[part_j]);

                    let path = configuration.bond_path(part_i, part_j);
                    let info = self.restriction.information(path);

                    old_energy += self.real_space_energy_pair(info, qi * qj, old_r);
                    new_energy += self.real_space_energy_pair(info, qi * qj, new_r);
                }
            }
        }

        return new_energy - old_energy;
    }
}

/// Self-interaction correction
impl Ewald {
    /// Self-interaction contribution to the energy
    fn self_energy(&self, configuration: &Configuration) -> f64 {
        let q2 = configuration.particles()
                              .charge
                              .iter()
                              .map(|q| q * q)
                              .sum::<f64>();
        return -self.alpha / f64::sqrt(PI) * q2 / FOUR_PI_EPSILON_0;
    }
}


/// k-space part of the summation
impl Ewald {
    /// Compute the Fourier transform of the electrostatic density
    fn eik_dot_r(&mut self, configuration: &Configuration) {
        let natoms = configuration.size();
        let range = -self.kmax..(self.kmax + 1);
        self.eikr.resize_if_different((range, 3, natoms));
        self.rho.clear();

        let positions = configuration.particles().position;
        let charges = configuration.particles().charge;

        // do the k = -1, 0, 1 cases first
        for spatial in 0..3 {
            let mut k_idx = [0.0, 0.0, 0.0];
            k_idx[spatial] = 1.0;
            let k_vector = configuration.cell.k_vector(k_idx);
            for i in 0..natoms {
                self.eikr[(0, spatial, i)] = Complex::cartesian(1.0, 0.0);
                self.eikr[(1, spatial, i)] = Complex::polar(1.0, k_vector * positions[i]);
                self.eikr[(-1, spatial, i)] = self.eikr[(1, spatial, i)].conj();
            }
        }

        // compute the other values of k by recursion
        for spatial in 0..3 {
            for k in 2..(self.kmax + 1) {
                for i in 0..natoms {
                    self.eikr[(k, spatial, i)] = self.eikr[(k - 1, spatial, i)] * self.eikr[(1, spatial, i)];
                    self.eikr[(-k, spatial, i)] = self.eikr[(k, spatial, i)].conj();
                }
            }
        }

        for &(ikx, iky, ikz) in &self.factors.index {
            let mut partial = Complex::zero();
            for i in 0..natoms {
                let phi = self.eikr[(ikx, 0, i)] *
                          self.eikr[(iky, 1, i)] *
                          self.eikr[(ikz, 2, i)];
                partial += charges[i] * phi;
            }
            self.rho.push(partial);
        }
    }

    /// k-space contribution to the energy
    fn k_space_energy(&mut self, configuration: &Configuration) -> f64 {
        self.eik_dot_r(configuration);

        let energy = self.factors.energy
            .par_iter()
            .zip_eq(&self.rho)
            .map(|(factor, rho)| factor * rho.norm2())
            .sum::<f64>();

        return energy / FOUR_PI_EPSILON_0;
    }

    /// k-space contribution to the forces
    fn k_space_forces(&mut self, configuration: &Configuration, forces: &mut [Vector3D]) {
        assert_eq!(forces.len(), configuration.size());
        self.eik_dot_r(configuration);

        let natoms = configuration.size();
        self.field.clear();
        self.field.resize(natoms, Vector3D::zero());

        let thread_local_field = ThreadLocalVec::with_size(natoms);
        self.factors.index
            .par_iter()
            .zip_eq(&self.factors.field)
            .zip_eq(&self.rho)
            .for_each(|((&(ikx, iky, ikz), factor), rho)| {
                let mut field = thread_local_field.borrow_mut();
                for i in 0..natoms {
                    let eikr = self.eikr[(ikx, 0, i)] *
                               self.eikr[(iky, 1, i)] *
                               self.eikr[(ikz, 2, i)];
                    let partial = eikr * rho.conj();
                    field[i] += partial.imag() * factor;
                }
            });

            thread_local_field.sum_into(&mut self.field);

        let charges = configuration.particles().charge;
        for (force, &charge, field) in zip!(&mut *forces, charges, &self.field) {
            *force += charge * field / FOUR_PI_EPSILON_0;
        }
    }

    /// k-space contribution to the atomic virial
    fn k_space_atomic_virial(&mut self, configuration: &Configuration) -> Matrix3 {
        self.eik_dot_r(configuration);

        let virial = self.factors.virial
            .par_iter()
            .zip_eq(&self.rho)
            .map(|(factor, rho)| rho.norm2() * factor)
            .sum::<Matrix3>();

        return virial / FOUR_PI_EPSILON_0;
    }

    /// k-space contribution to the molecular virial
    fn k_space_molecular_virial(&mut self, configuration: &Configuration) -> Matrix3 {
        let atomic = self.k_space_atomic_virial(configuration);

        let mut forces = vec![Vector3D::zero(); configuration.size()];
        self.k_space_forces(configuration, &mut forces);

        let positions = configuration.particles().position;
        let mut correction = Matrix3::zero();
        for molecule in configuration.molecules() {
            let com = molecule.center_of_mass();
            for i in molecule.indexes() {
                let di = positions[i] - com;
                correction += forces[i].tensorial(&di);
            }
        }

        return atomic - correction;
    }

    /// Compute the Fourier transform of the electrostatic density changes
    /// while moving the molecule with the given `molecule_id` to
    /// `new_positions`
    fn delta_rho_move_rigid_molecules(
        &mut self,
        configuration: &Configuration,
        molecule_id: usize,
        new_positions: &[Vector3D],
    ) -> Vec<Complex> {
        let molecule = configuration.molecule(molecule_id);
        let mut new_energy_ikr = Ewald3DArray::zeros((-self.kmax..(self.kmax + 1), 3, molecule.size()));

        // Do the k=0, 1 cases first
        for spatial in 0..3 {
            let mut k_idx = [0.0, 0.0, 0.0];
            k_idx[spatial] = 1.0;
            let k_vector = configuration.cell.k_vector(k_idx);
            for i in 0..molecule.size() {
                new_energy_ikr[(0, spatial, i)] = Complex::cartesian(1.0, 0.0);
                new_energy_ikr[(1, spatial, i)] = Complex::polar(1.0, k_vector * new_positions[i]);
                new_energy_ikr[(-1, spatial, i)] = new_energy_ikr[(1, spatial, i)].conj();
            }
        }

        // Use recursive definition for computing the factor for all the other values of k.
        for spatial in 0..3 {
            for k in 2..(self.kmax + 1) {
                for i in 0..molecule.size() {
                    new_energy_ikr[(k, spatial, i)] = new_energy_ikr[(k - 1, spatial, i)] * new_energy_ikr[(1, spatial, i)];
                    new_energy_ikr[(-k, spatial, i)] = new_energy_ikr[(k, spatial, i)].conj();
                }
            }
        }

        let mut delta = Vec::new();
        let charges = configuration.particles().charge;
        for &(ikx, iky, ikz) in &self.factors.index {
            let mut partial = Complex::zero();
            for (i, part_i) in molecule.indexes().enumerate() {
                let old_phi = self.eikr[(ikx, 0, part_i)] *
                              self.eikr[(iky, 1, part_i)] *
                              self.eikr[(ikz, 2, part_i)];

                let new_phi = new_energy_ikr[(ikx, 0, i)] *
                              new_energy_ikr[(iky, 1, i)] *
                              new_energy_ikr[(ikz, 2, i)];

                partial += charges[part_i] * (new_phi - old_phi);
            }
            delta.push(partial);
        }

        return delta;
    }

    fn k_space_move_molecule_cost(
        &mut self,
        configuration: &Configuration,
        molecule_id: usize,
        new_positions: &[Vector3D],
    ) -> f64 {
        let mut old_energy = 0.0;
        for (factor, &rho) in zip!(&self.factors.energy, &self.rho) {
            old_energy += factor * rho.norm2();
        }
        old_energy /= FOUR_PI_EPSILON_0;

        let delta_rho = self.delta_rho_move_rigid_molecules(
            configuration, molecule_id, new_positions
        );

        let mut new_energy = 0.0;
        for (factor, &rho, &delta) in zip!(&self.factors.energy, &self.rho, &delta_rho) {
            new_energy += factor * (rho + delta).norm2();
        }
        new_energy /= FOUR_PI_EPSILON_0;

        self.updater = Some(Box::new(move |ewald: &mut Ewald| {
            for (rho, &delta) in zip!(&mut ewald.rho, &delta_rho) {
                *rho += delta;
            }
        }));

        return new_energy - old_energy;
    }
}

/// Thread-sade wrapper around Ewald implementing `CoulombicPotential`.
///
/// This wrapper allow to share a Ewald solver between threads (make it `Send
/// + Sync`) while still using caching in Monte Carlo simulations (with
/// interior mutability).
pub struct SharedEwald(RwLock<Ewald>);

impl SharedEwald {
    /// Wrap `ewald` in a thread-safe structure.
    ///
    /// # Example
    /// ```
    /// # use lumol_core::energy::{Ewald, SharedEwald, CoulombicPotential};
    /// let ewald = SharedEwald::new(Ewald::new(12.5, 10, None));
    /// let boxed: Box<dyn CoulombicPotential> = Box::new(ewald);
    /// ```
    pub fn new(ewald: Ewald) -> SharedEwald {
        SharedEwald(RwLock::new(ewald))
    }

    /// Get read access to the underlying Ewald solver
    fn read(&self) -> RwLockReadGuard<'_, Ewald> {
        // The lock should never be poisoned, because any panic will unwind
        // and finish the simulation.
        self.0.read().expect("Ewald lock is poisoned")
    }

    /// Get write access to the underlying Ewald solver
    fn write(&self) -> RwLockWriteGuard<'_, Ewald> {
        // The lock should never be poisoned, because any panic will unwind
        // and finish the simulation.
        self.0.write().expect("Ewald lock is poisoned")
    }
}

impl Clone for SharedEwald {
    fn clone(&self) -> SharedEwald {
        SharedEwald::new(self.read().clone())
    }
}

impl GlobalPotential for SharedEwald {
    fn cutoff(&self) -> Option<f64> {
        Some(self.read().rc)
    }

    fn energy(&self, configuration: &Configuration) -> f64 {
        let mut ewald = self.write();
        ewald.prepare(&configuration.cell);
        let real = ewald.real_space_energy(configuration);
        let self_e = ewald.self_energy(configuration);
        let k_space = ewald.k_space_energy(configuration);
        return real + self_e + k_space;
    }

    fn forces(&self, configuration: &Configuration, forces: &mut [Vector3D])  {
        assert_eq!(forces.len(), configuration.size());
        let mut ewald = self.write();
        ewald.prepare(&configuration.cell);

        ewald.real_space_forces(configuration, forces);
        // No self force
        ewald.k_space_forces(configuration, forces);
    }

    fn atomic_virial(&self, configuration: &Configuration) -> Matrix3 {
        let mut ewald = self.write();
        ewald.prepare(&configuration.cell);
        let real = ewald.real_space_atomic_virial(configuration);
        // No self virial
        let k_space = ewald.k_space_atomic_virial(configuration);
        return real + k_space;
    }

    fn molecular_virial(&self, configuration: &Configuration) -> Matrix3 {
        let mut ewald = self.write();
        ewald.prepare(&configuration.cell);
        let real = ewald.real_space_molecular_virial(configuration);
        // No self virial
        let k_space = ewald.k_space_molecular_virial(configuration);
        return real + k_space;
    }
}

impl CoulombicPotential for SharedEwald {
    fn set_restriction(&mut self, restriction: PairRestriction) {
        self.write().restriction = restriction;
    }
}

impl GlobalCache for SharedEwald {
    fn move_molecule_cost(
        &self,
        configuration: &Configuration,
        molecule_id: usize,
        new_positions: &[Vector3D]
    ) -> f64 {
        let mut ewald = self.write();
        ewald.prepare(&configuration.cell);
        let real = ewald.real_space_move_molecule_cost(configuration, molecule_id, new_positions);
        /* No self cost */
        let k_space = ewald.k_space_move_molecule_cost(configuration, molecule_id, new_positions);
        return real + k_space;
    }

    fn update(&self) {
        let mut ewald = self.write();
        if ewald.updater.is_some() {
            let mut updater = None;
            ::std::mem::swap(&mut updater, &mut ewald.updater);
            let updater = updater.unwrap();
            updater(&mut ewald);
        }
    }
}

#[cfg(test)]
mod tests {
    pub use super::*;
    use crate::System;
    use crate::utils::system_from_xyz;

    pub fn nacl_pair() -> System {
        let mut system = system_from_xyz("2
        cell: 20.0
        Cl 0.0 0.0 0.0
        Na 1.5 0.0 0.0
        ");
        system.particles_mut().charge[0] = -1.0;
        system.particles_mut().charge[1] = 1.0;
        return system;
    }

    pub fn water() -> System {
        let mut system = system_from_xyz("3
        cell: 20.0
        O  0.0  0.0  0.0
        H -0.7 -0.7  0.3
        H  0.3 -0.3 -0.8
        ");
        assert!(system.add_bond(0, 1).is_empty());
        assert!(system.add_bond(0, 2).is_empty());
        assert!(system.molecules().count() == 1);

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
    fn with_accuracy() {
        let ewald = Ewald::with_accuracy(8.5, 1e-6, &water());
        assert!(f64::abs(ewald.alpha - 0.2998) < 1e-4);
        assert_eq!(ewald.kmax, 5);
    }

    mod errors {
        use super::*;
        use crate::GlobalPotential;
        use crate::UnitCell;

        #[test]
        #[should_panic(expected="Ewald is not defined with infinite unit cell")]
        fn infinite_cell() {
            let mut system = nacl_pair();
            system.cell = UnitCell::infinite();
            let ewald = SharedEwald::new(Ewald::new(8.0, 10, None));
            let _ = ewald.energy(&system);
        }

        #[test]
        #[should_panic(expected="the cutoff can not be negative in Ewald")]
        fn negative_cutoff() {
            let _ = Ewald::new(-8.0, 10, None);
        }

        #[test]
        #[should_panic(expected="alpha can not be negative in Ewald")]
        fn negative_alpha() {
            let _ = Ewald::new(8.0, 10, -45.2);
        }

        #[test]
        #[should_panic(expected="kmax can not be 0 in Ewald")]
        fn kmax_null() {
            let _ = Ewald::new(8.0, 0, None);
        }
    }

    mod pairs {
        use super::*;
        use crate::GlobalPotential;
        use approx::{assert_ulps_eq, assert_relative_eq};

        #[test]
        #[allow(clippy::unreadable_literal)]
        fn energy() {
            let system = nacl_pair();
            let ewald = SharedEwald::new(Ewald::new(8.0, 10, None));

            let energy = ewald.energy(&system);
            // This was computed by hand
            let energy_brute_force = -0.09262397663346732;
            assert_ulps_eq!(energy, energy_brute_force, epsilon=1e-4);

            // Just checking that this does not crash
            let ewald = SharedEwald::new(Ewald::new(8.0, 1, None));
            let _ = ewald.energy(&system);
        }

        #[test]
        fn real_forces_finite_differences() {
            let mut system = nacl_pair();
            let mut ewald = Ewald::new(8.0, 10, None);
            ewald.prepare(&system.cell);

            let e = ewald.real_space_energy(&system);
            let eps = 1e-9;
            system.particles_mut().position[0][0] += eps;

            let e1 = ewald.real_space_energy(&system);
            let mut forces = vec![Vector3D::zero(); 2];
            ewald.real_space_forces(&system, &mut forces);
            assert_relative_eq!((e - e1) / eps, forces[0][0], epsilon=1e-6);
        }

        #[test]
        fn k_space_forces_finite_differences() {
            let mut system = nacl_pair();
            // Using a small cutoff to increase the weight of k-space
            let mut ewald = Ewald::new(2.0, 10, None);
            ewald.prepare(&system.cell);

            let e = ewald.k_space_energy(&system);
            let eps = 1e-9;
            system.particles_mut().position[0][0] += eps;

            let e1 = ewald.k_space_energy(&system);
            let mut forces = vec![Vector3D::zero(); 2];
            ewald.k_space_forces(&system, &mut forces);
            assert_relative_eq!((e - e1) / eps, forces[0][0], epsilon=1e-6);
        }

        #[test]
        fn total_forces() {
            let mut system = nacl_pair();
            let ewald = SharedEwald::new(Ewald::new(8.0, 10, None));

            let mut forces = vec![Vector3D::zero(); 2];
            ewald.forces(&system, &mut forces);
            let norm = (forces[0] + forces[1]).norm();
            // Total force should be null
            assert_ulps_eq!(norm, 0.0);

            // Force is attractive
            assert!(forces[0][0] > 0.0);
            assert!(forces[1][0] < 0.0);

            // Finite difference computation of the force
            let e = ewald.energy(&system);
            let eps = 1e-9;
            system.particles_mut().position[0][0] += eps;

            let e1 = ewald.energy(&system);
            let mut forces = vec![Vector3D::zero(); 2];
            ewald.forces(&system, &mut forces);
            assert_relative_eq!((e - e1) / eps, forces[0][0], epsilon=1e-6);
        }
    }

    mod molecules {
        use super::*;
        use crate::Vector3D;
        use crate::{GlobalPotential, PairRestriction, CoulombicPotential};

        use approx::{assert_ulps_eq, assert_relative_eq};

        #[test]
        #[allow(clippy::unreadable_literal)]
        fn energy() {
            let system = water();
            let mut ewald = SharedEwald::new(Ewald::new(8.0, 10, None));
            ewald.set_restriction(PairRestriction::InterMolecular);

            let energy = ewald.energy(&system);
            let expected = -0.000009243868813825495;
            assert_ulps_eq!(energy, expected);
        }

        #[test]
        fn real_space_forces_finite_differences() {
            let mut system = water();
            let mut ewald = Ewald::new(8.0, 10, None);
            ewald.restriction = PairRestriction::InterMolecular;
            ewald.prepare(&system.cell);

            let mut forces = vec![Vector3D::zero(); 3];
            ewald.real_space_forces(&system, &mut forces);
            let force = forces[0][0];

            let eps = 1e-9;
            let e = ewald.real_space_energy(&system);
            system.particles_mut().position[0][0] += eps;
            let e1 = ewald.real_space_energy(&system);

            assert_relative_eq!((e - e1) / eps, force, epsilon = 1e-6);
        }

        #[test]
        fn k_space_forces_finite_differences() {
            let mut system = water();
            let mut ewald = Ewald::new(8.0, 10, None);
            ewald.restriction = PairRestriction::InterMolecular;
            ewald.prepare(&system.cell);

            let mut forces = vec![Vector3D::zero(); 3];
            ewald.k_space_forces(&system, &mut forces);
            let force = forces[0][0];

            let eps = 1e-9;
            let e = ewald.k_space_energy(&system);
            system.particles_mut().position[0][0] += eps;
            let e1 = ewald.k_space_energy(&system);

            assert_relative_eq!((e - e1) / eps, force, epsilon = 1e-6);
        }

        #[test]
        fn total_forces() {
            let mut system = water();
            let mut ewald = SharedEwald::new(Ewald::new(8.0, 10, None));
            ewald.set_restriction(PairRestriction::InterMolecular);

            let mut forces = vec![Vector3D::zero(); 3];
            ewald.forces(&system, &mut forces);
            let norm = (forces[0] + forces[1]).norm();
            // Total force should be null
            assert!(norm.abs() < 1e-3);

            let eps = 1e-9;
            let e = ewald.energy(&system);
            system.particles_mut().position[0][0] += eps;
            let e1 = ewald.energy(&system);

            let force = forces[0][0];
            assert_relative_eq!((e - e1) / eps, force, epsilon = 1e-6);
        }
    }

    mod atomic_virial {
        use super::*;
        use crate::Matrix3;
        use crate::PairRestriction;

        use approx::assert_relative_eq;

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

        #[test]
        fn virial_is_energy() {
            // A nice property of Ewald summation for point charge systems
            let system = water();
            let ewald = SharedEwald::new(Ewald::new(8.0, 10, None));

            let energy = ewald.energy(&system);
            let virial = ewald.atomic_virial(&system).trace();
            assert_relative_eq!(energy, virial, max_relative = 1e-3);

            let system = nacl_pair();
            let ewald = SharedEwald::new(Ewald::new(8.0, 10, None));

            let energy = ewald.energy(&system);
            let virial = ewald.atomic_virial(&system).trace();
            assert_relative_eq!(energy, virial, max_relative = 1e-3);
        }

        #[test]
        fn real_space_finite_differences() {
            let mut system = water();
            let mut ewald = Ewald::new(8.0, 10, None);
            ewald.restriction = PairRestriction::InterMolecular;
            ewald.prepare(&system.cell);

            let eps = 1e-9;
            let virial = ewald.real_space_atomic_virial(&system);
            let mut finite_diff = Matrix3::zero();

            for i in 0..3 {
                for j in 0..3 {
                    ewald.prepare(&system.cell);
                    let e = ewald.real_space_energy(&system);
                    scale(&mut system, i, j, eps);
                    ewald.prepare(&system.cell);
                    let e1 = ewald.real_space_energy(&system);
                    finite_diff[i][j] = (e - e1) / eps;
                }
            }

            assert_relative_eq!(virial, finite_diff, epsilon = 1e-6);
        }

        #[test]
        fn k_space_finite_differences() {
            let mut system = water();
            let mut ewald = Ewald::new(2.0, 10, None);
            ewald.restriction = PairRestriction::InterMolecular;
            ewald.prepare(&system.cell);

            let eps = 1e-9;
            let virial = ewald.k_space_atomic_virial(&system);
            let mut finite_diff = Matrix3::zero();

            for i in 0..3 {
                for j in 0..3 {
                    ewald.prepare(&system.cell);
                    let e = ewald.k_space_energy(&system);
                    scale(&mut system, i, j, eps);
                    ewald.prepare(&system.cell);
                    let e1 = ewald.k_space_energy(&system);
                    finite_diff[i][j] = (e - e1) / eps;
                }
            }

            // Make sure the finite_diff matrix is symmetric
            finite_diff = (finite_diff + finite_diff.transposed()) / 2.0;
            assert_relative_eq!(virial, finite_diff, epsilon = 1e-6);
        }
    }

    use approx::assert_relative_eq;

    #[test]
    #[allow(clippy::items_after_statements)]
    fn move_molecule() {
        let mut system = system_from_xyz("6
        cell: 20.0
        H  0.3 -0.3 -0.8
        O  0.0  0.0  0.0
        H -0.7 -0.7  0.3
        H  2.3  1.7 -0.8
        O  2.0  2.0  0.0
        H  1.3  1.3  0.3
        ");
        assert!(system.add_bond(0, 1).is_empty());
        assert!(system.add_bond(1, 2).is_empty());
        assert!(system.add_bond(3, 4).is_empty());
        assert!(system.add_bond(4, 5).is_empty());
        assert!(system.molecules().count() == 2);

        for particle in system.particles_mut() {
            if particle.name == "O" {
                *particle.charge = -0.8476;
            } else if particle.name == "H" {
                *particle.charge = 0.4238;
            }
        }

        type CostCompute = fn (ewald: &SharedEwald, system: &System, molecule: usize, positions: &[Vector3D]) -> f64;
        type EnergyCompute = fn (ewald: &SharedEwald, system: &System) -> f64;

        #[allow(clippy::unreadable_literal)]
        let check_cache = |mut system: System, ewald: Ewald, compute_energy: EnergyCompute, compute_cost: CostCompute| {
            let mut ewald = SharedEwald::new(ewald);
            ewald.set_restriction(PairRestriction::InterMolecular);
            ewald.write().prepare(&system.cell);

            let check = ewald.clone();
            // Initialize cached values
            let _ = compute_energy(&ewald, &system);
            let old_energy = compute_energy(&check, &system);

            let new_positions = &[
                Vector3D::new(0.41727, 2.29401, -0.0558),
                Vector3D::new(0.5097743599026461, 3.194114034722624, -0.020364564697826326),
                Vector3D::new(-0.2501317777731211, 3.562366060753896, -0.6178033542374419),
            ];
            let cost = compute_cost(&ewald, &system, 0, new_positions);

            system.particles_mut().position[0] = new_positions[0];
            system.particles_mut().position[1] = new_positions[1];
            system.particles_mut().position[2] = new_positions[2];
            let new_energy = compute_energy(&check, &system);
            assert_relative_eq!(cost, new_energy - old_energy, max_relative = 1e-12);
        };

        // Real space energy
        check_cache(
            system.clone(),
            Ewald::new(8.0, 10, None),
            |ewald, system| {
                ewald.read().real_space_energy(system)
            },
            |ewald, system, molecule, positions| {
                ewald.read().real_space_move_molecule_cost(system, molecule, positions)
            }
        );

        // k-space energy
        check_cache(
            system.clone(),
            Ewald::new(2.0, 10, None),
            |ewald, system| {
                ewald.write().k_space_energy(system)
            },
            |ewald, system, molecule, positions| {
                ewald.write().k_space_move_molecule_cost(system, molecule, positions)
            }
        );

        // Whole energy at once
        check_cache(
            system,
            Ewald::new(8.0, 10, None),
            |ewald, system| {
                ewald.energy(system)
            },
            |ewald, system, molecule, positions| {
                ewald.move_molecule_cost(system, molecule, positions)
            }
        );
    }

    // Comparing the value for each component of Ewald energy with the NIST
    // reference. See `tests/nist-spce.rs` for more information. These tests
    // check values that are not accessible from the outside of lumol-core.
    mod nist {
        #![allow(clippy::unreadable_literal)]
        use super::*;

        use std::path::Path;
        use std::fs::File;
        use std::io::Read;

        pub fn get_system(path: &str) -> System {
            let path = Path::new(env!("CARGO_MANIFEST_DIR"))
                .join("..")
                .join("tests")
                .join("data")
                .join("nist-spce")
                .join(path);

            let mut file = File::open(path).unwrap();
            let mut buffer = String::new();
            let _ = file.read_to_string(&mut buffer).unwrap();

            let mut system = system_from_xyz(&buffer);

            for i in 0..system.size() {
                if i % 3 == 0 {
                    let _ = system.add_bond(i, i + 1);
                    let _ = system.add_bond(i, i + 2);
                }
            }

            for particle in system.particles_mut() {
                match particle.name.as_ref() {
                    "H" => *particle.charge = 0.42380,
                    "O" => *particle.charge = -2.0 * 0.42380,
                    other => panic!("Unknown particle name: {other}"),
                }
            }

            return system;
        }

        #[test]
        fn cells() {
            let system = get_system("spce-1.xyz");
            assert_eq!(system.cell.a(), 20.0);
            assert_eq!(system.cell.b(), 20.0);
            assert_eq!(system.cell.c(), 20.0);

            let system = get_system("spce-2.xyz");
            assert_eq!(system.cell.a(), 20.0);
            assert_eq!(system.cell.b(), 20.0);
            assert_eq!(system.cell.c(), 20.0);

            let system = get_system("spce-3.xyz");
            assert_eq!(system.cell.a(), 20.0);
            assert_eq!(system.cell.b(), 20.0);
            assert_eq!(system.cell.c(), 20.0);

            let system = get_system("spce-4.xyz");
            assert_eq!(system.cell.a(), 30.0);
            assert_eq!(system.cell.b(), 30.0);
            assert_eq!(system.cell.c(), 30.0);
        }

        #[test]
        fn with_accuracy() {
            let ewald = Ewald::with_accuracy(9.0, 1e-5, &get_system("spce-1.xyz"));
            assert!(f64::abs(ewald.alpha - 0.2826) < 1e-4);
            assert_eq!(ewald.kmax, 5);

            let ewald = Ewald::with_accuracy(9.0, 1e-5, &get_system("spce-2.xyz"));
            assert!(f64::abs(ewald.alpha - 0.2900) < 1e-4);
            assert_eq!(ewald.kmax, 5);

            let ewald = Ewald::with_accuracy(9.0, 1e-5, &get_system("spce-3.xyz"));
            assert!(f64::abs(ewald.alpha - 0.2943) < 1e-4);
            assert_eq!(ewald.kmax, 5);

            let ewald = Ewald::with_accuracy(9.0, 1e-5, &get_system("spce-4.xyz"));
            assert!(f64::abs(ewald.alpha - 0.2912) < 1e-4);
            assert_eq!(ewald.kmax, 8);
        }

        #[allow(clippy::unreadable_literal)]
        mod cutoff_9 {
            use super::*;
            use crate::consts::K_BOLTZMANN;
            use crate::units;
            const CUTOFF: f64 = 9.0;

            #[test]
            fn nist1() {
                let system = get_system("spce-1.xyz");

                let alpha = 5.6 / system.cell.a();
                let mut ewald = Ewald::new(CUTOFF, 5, alpha);
                ewald.restriction = PairRestriction::InterMolecular;
                ewald.prepare(&system.cell);

                let energy = ewald.real_space_energy(&system) / K_BOLTZMANN;
                let expected = 2.251086e6;
                assert_relative_eq!(energy, expected, max_relative = 1e-4);

                let energy = ewald.k_space_energy(&system) / K_BOLTZMANN;
                let expected = 6.27009e3;
                assert_relative_eq!(energy, expected, max_relative = 2e-3);

                let energy = ewald.self_energy(&system) / K_BOLTZMANN;
                let expected = -2.84469e6;
                assert_relative_eq!(energy, expected, max_relative = 1e-4);
            }

            #[test]
            fn nist1_virial() {
                let system = get_system("spce-1.xyz");
                let mut ewald = Ewald::new(CUTOFF, 8, 0.364209);
                ewald.restriction = PairRestriction::InterMolecular;
                ewald.prepare(&system.cell);

                let convert = units::from(1.0, "atm").unwrap() * system.volume();

                let virial = ewald.real_space_atomic_virial(&system);
                let expected = convert * Matrix3::new([
                    [-3852.8846,   264.92131,  15.263331],
                    [ 264.92131,  -3778.4993, -108.79484],
                    [ 15.263331,  -108.79484,  -1885.924],
                ]);
                assert_relative_eq!(virial, expected, max_relative = 1e-3);

                let virial = ewald.k_space_atomic_virial(&system);
                let expected = convert * Matrix3::new([
                    [-124.54878,  5.3873858,   15.093282],
                    [ 5.3873858, -111.92329,  -36.333371],
                    [ 15.093282, -36.333371, -248.12515 ],
                ]);
                assert_relative_eq!(virial, expected, max_relative = 1e-3);

                let ewald = SharedEwald::new(ewald);
                let energy = ewald.energy(&system);
                let virial = ewald.atomic_virial(&system).trace();
                assert_relative_eq!(energy, virial, max_relative = 1e-3);
            }

            #[test]
            fn nist2() {
                let system = get_system("spce-2.xyz");

                let alpha = 5.6 / system.cell.a();
                let mut ewald = Ewald::new(CUTOFF, 5, alpha);
                ewald.restriction = PairRestriction::InterMolecular;
                ewald.prepare(&system.cell);

                let energy = ewald.real_space_energy(&system) / K_BOLTZMANN;
                let expected = 4.4269e6;
                assert_relative_eq!(energy, expected, max_relative = 1e-4);

                let energy = ewald.k_space_energy(&system) / K_BOLTZMANN;
                let expected = 6.03495e3;
                assert_relative_eq!(energy, expected, max_relative = 5e-3);

                let energy = ewald.self_energy(&system) / K_BOLTZMANN;
                let expected = -5.68938e6;
                assert_relative_eq!(energy, expected, max_relative = 1e-4);
            }

            #[test]
            fn nist2_virial() {
                let system = get_system("spce-2.xyz");
                let mut ewald = Ewald::new(CUTOFF, 8, 0.370036);
                ewald.restriction = PairRestriction::InterMolecular;
                ewald.prepare(&system.cell);

                let convert = units::from(1.0, "atm").unwrap() * system.volume();

                let virial = ewald.real_space_atomic_virial(&system);
                let expected = convert * Matrix3::new([
                    [-6388.8957, -158.39688, -344.71965],
                    [-158.39688, -6988.116,   443.24826],
                    [-344.71965,  443.24826, -7169.9086],
                ]);
                assert_relative_eq!(virial, expected, max_relative = 1e-3);

                let virial = ewald.k_space_atomic_virial(&system);
                let expected = convert * Matrix3::new([
                    [-296.86013, -10.296217, -3.0451399],
                    [-10.296217, -250.90827, -7.9310376],
                    [-3.0451399, -7.9310376, -299.97078],
                ]);
                assert_relative_eq!(virial, expected, max_relative = 1e-3);

                let ewald = SharedEwald::new(ewald);
                let energy = ewald.energy(&system);
                let virial = ewald.atomic_virial(&system).trace();
                assert_relative_eq!(energy, virial, max_relative = 1e-3);
            }

            #[test]
            fn nist3() {
                let system = get_system("spce-3.xyz");

                let alpha = 5.6 / system.cell.a();
                let mut ewald = Ewald::new(CUTOFF, 5, alpha);
                ewald.restriction = PairRestriction::InterMolecular;
                ewald.prepare(&system.cell);

                let energy = ewald.real_space_energy(&system) / K_BOLTZMANN;
                let expected = 6.46678e6;
                assert_relative_eq!(energy, expected, max_relative = 1e-4);

                let energy = ewald.k_space_energy(&system) / K_BOLTZMANN;
                let expected = 5.24461e3;
                assert_relative_eq!(energy, expected, max_relative = 2e-3);

                let energy = ewald.self_energy(&system) / K_BOLTZMANN;
                let expected = -8.53407e6;
                assert_relative_eq!(energy, expected, max_relative = 1e-4);
            }

            #[test]
            fn nist3_virial() {
                let system = get_system("spce-3.xyz");
                let mut ewald = Ewald::new(CUTOFF, 8, 0.373403);
                ewald.restriction = PairRestriction::InterMolecular;
                ewald.prepare(&system.cell);

                let convert = units::from(1.0, "atm").unwrap() * system.volume();

                let virial = ewald.real_space_atomic_virial(&system);
                let expected = convert * Matrix3::new([
                    [-10038.765, -1366.2912,  95.727168],
                    [-1366.2912, -12683.884, -213.42417],
                    [ 95.727168, -213.42417, -11624.416],
                ]);
                assert_relative_eq!(virial, expected, max_relative = 1e-3);

                let virial = ewald.k_space_atomic_virial(&system);
                let expected = convert * Matrix3::new([
                    [-255.20669,  18.680095, -37.1178  ],
                    [ 18.680095, -275.89451, -5.1062841],
                    [  -37.1178, -5.1062841, -232.9918 ],
                ]);
                assert_relative_eq!(virial, expected, max_relative = 1e-3);

                let ewald = SharedEwald::new(ewald);
                let energy = ewald.energy(&system);
                let virial = ewald.atomic_virial(&system).trace();
                assert_relative_eq!(energy, virial, max_relative = 1e-3);
            }

            #[test]
            fn nist4() {
                let system = get_system("spce-4.xyz");

                let alpha = 5.6 / system.cell.a();
                let mut ewald = Ewald::new(CUTOFF, 5, alpha);
                ewald.restriction = PairRestriction::InterMolecular;
                ewald.prepare(&system.cell);

                let energy = ewald.real_space_energy(&system) / K_BOLTZMANN;
                let expected = 1.07011e7;
                assert_relative_eq!(energy, expected, max_relative = 1e-4);

                let energy = ewald.k_space_energy(&system) / K_BOLTZMANN;
                let expected = 7.58785e3;
                assert_relative_eq!(energy, expected, max_relative = 2e-3);

                let energy = ewald.self_energy(&system) / K_BOLTZMANN;
                let expected = -1.42235e7;
                assert_relative_eq!(energy, expected, max_relative = 1e-4);
            }

            #[test]
            fn nist4_virial() {
                let system = get_system("spce-4.xyz");
                let mut ewald = Ewald::new(CUTOFF, 12, 0.370914);
                ewald.restriction = PairRestriction::InterMolecular;
                ewald.prepare(&system.cell);

                let convert = units::from(1.0, "atm").unwrap() * system.volume();

                let virial = ewald.real_space_atomic_virial(&system);
                let expected = convert * Matrix3::new([
                    [-5586.0358, -136.58999, -30.075826],
                    [-136.58999, -5846.9018, -121.54568],
                    [-30.075826, -121.54568, -4988.8238],
                ]);
                assert_relative_eq!(virial, expected, max_relative = 5e-3);

                let virial = ewald.k_space_atomic_virial(&system);
                let expected = convert * Matrix3::new([
                    [ -482.8217, -17.877365,  11.908942],
                    [-17.877365, -465.33957,  3.5274012],
                    [ 11.908942,  3.5274012, -536.23714],
                ]);
                assert_relative_eq!(virial, expected, max_relative = 1e-3);

                let ewald = SharedEwald::new(ewald);
                let energy = ewald.energy(&system);
                let virial = ewald.atomic_virial(&system).trace();
                assert_relative_eq!(energy, virial, max_relative = 1e-3);
            }
        }

        #[allow(clippy::unreadable_literal)]
        mod cutoff_10 {
            use super::*;
            use crate::consts::K_BOLTZMANN;
            use crate::units;
            const CUTOFF: f64 = 10.0;

            #[test]
            fn nist1() {
                let system = get_system("spce-1.xyz");

                let alpha = 5.6 / system.cell.a();
                let mut ewald = Ewald::new(CUTOFF, 5, alpha);
                ewald.restriction = PairRestriction::InterMolecular;
                ewald.prepare(&system.cell);

                let energy = ewald.real_space_energy(&system) / K_BOLTZMANN;
                let expected = 2.251101e6;
                assert_relative_eq!(energy, expected, max_relative = 1e-4);

                let energy = ewald.k_space_energy(&system) / K_BOLTZMANN;
                let expected = 6.27009e3;
                assert_relative_eq!(energy, expected, max_relative = 2e-3);

                let energy = ewald.self_energy(&system) / K_BOLTZMANN;
                let expected = -2.84469e6;
                assert_relative_eq!(energy, expected, max_relative = 1e-4);
            }

            #[test]
            fn nist1_virial() {
                let system = get_system("spce-1.xyz");
                let mut ewald = Ewald::new(CUTOFF, 7, 0.326983);
                ewald.restriction = PairRestriction::InterMolecular;
                ewald.prepare(&system.cell);

                let convert = units::from(1.0, "atm").unwrap() * system.volume();

                let virial = ewald.real_space_atomic_virial(&system);
                let expected = convert * Matrix3::new([
                    [-3916.4836,  268.02381,  23.686442],
                    [ 268.02381, -3839.2663, -118.84454],
                    [ 23.686442, -118.84454, -1945.1106],
                ]);
                assert_relative_eq!(virial, expected, max_relative = 1e-3);

                let virial = ewald.k_space_atomic_virial(&system);
                let expected = convert * Matrix3::new([
                    [-61.298853,  2.0228732,    6.70617],
                    [ 2.0228732, -51.134442,   -26.2246],
                    [   6.70617,   -26.2246, -188.90856],
                ]);
                assert_relative_eq!(virial, expected, max_relative = 1e-3);

                let ewald = SharedEwald::new(ewald);
                let energy = ewald.energy(&system);
                let virial = ewald.atomic_virial(&system).trace();
                assert_relative_eq!(energy, virial, max_relative = 1e-3);
            }

            #[test]
            fn nist2() {
                let system = get_system("spce-2.xyz");

                let alpha = 5.6 / system.cell.a();
                let mut ewald = Ewald::new(CUTOFF, 5, alpha);
                ewald.restriction = PairRestriction::InterMolecular;
                ewald.prepare(&system.cell);

                let energy = ewald.real_space_energy(&system) / K_BOLTZMANN;
                let expected = 4.42703e6;
                assert_relative_eq!(energy, expected, max_relative = 1e-4);

                let energy = ewald.k_space_energy(&system) / K_BOLTZMANN;
                let expected = 6.03495e3;
                assert_relative_eq!(energy, expected, max_relative = 5e-3);

                let energy = ewald.self_energy(&system) / K_BOLTZMANN;
                let expected = -5.68938e6;
                assert_relative_eq!(energy, expected, max_relative = 1e-4);
            }

            #[test]
            fn nist2_virial() {
                let system = get_system("spce-2.xyz");
                let mut ewald = Ewald::new(CUTOFF, 8, 0.332241);
                ewald.restriction = PairRestriction::InterMolecular;
                ewald.prepare(&system.cell);

                let convert = units::from(1.0, "atm").unwrap() * system.volume();

                let virial = ewald.real_space_atomic_virial(&system);
                let expected = convert * Matrix3::new([
                    [-6515.4904, -161.25038, -340.79656],
                    [-161.25038,  -7099.813,  440.80231],
                    [-340.79656,  440.80231, -7285.3392],
                ]);
                assert_relative_eq!(virial, expected, max_relative = 1e-3);

                let virial = ewald.k_space_atomic_virial(&system);
                let expected = convert * Matrix3::new([
                    [-171.16013, -7.5200295, -6.7910143],
                    [-7.5200295, -140.58412, -5.5476543],
                    [-6.7910143, -5.5476543, -185.27155],
                ]);
                assert_relative_eq!(virial, expected, max_relative = 1e-3);

                let ewald = SharedEwald::new(ewald);
                let energy = ewald.energy(&system);
                let virial = ewald.atomic_virial(&system).trace();
                assert_relative_eq!(energy, virial, max_relative = 1e-3);
            }

            #[test]
            fn nist3() {
                let system = get_system("spce-3.xyz");

                let alpha = 5.6 / system.cell.a();
                let mut ewald = Ewald::new(CUTOFF, 5, alpha);
                ewald.restriction = PairRestriction::InterMolecular;
                ewald.prepare(&system.cell);

                let energy = ewald.real_space_energy(&system) / K_BOLTZMANN;
                let expected = 6.46701e6;
                assert_relative_eq!(energy, expected, max_relative = 1e-4);

                let energy = ewald.k_space_energy(&system) / K_BOLTZMANN;
                let expected = 5.24461e3;
                assert_relative_eq!(energy, expected, max_relative = 2e-3);

                let energy = ewald.self_energy(&system) / K_BOLTZMANN;
                let expected = -8.53407e6;
                assert_relative_eq!(energy, expected, max_relative = 1e-4);
            }

            #[test]
            fn nist3_virial() {
                let system = get_system("spce-3.xyz");
                let mut ewald = Ewald::new(CUTOFF, 8, 0.335278);
                ewald.restriction = PairRestriction::InterMolecular;
                ewald.prepare(&system.cell);

                let convert = units::from(1.0, "atm").unwrap() * system.volume();

                let virial = ewald.real_space_atomic_virial(&system);
                let expected = convert * Matrix3::new([
                    [-10159.552, -1365.6921,  85.837259],
                    [-1365.6921, -12806.406, -202.68293],
                    [ 85.837259, -202.68293, -11744.811],
                ]);
                assert_relative_eq!(virial, expected, max_relative = 1e-3);

                let virial = ewald.k_space_atomic_virial(&system);
                let expected = convert * Matrix3::new([
                    [-136.28835,   18.49649, -27.386553],
                    [  18.49649, -155.56804, -15.667645],
                    [-27.386553, -15.667645, -114.48139],
                ]);
                assert_relative_eq!(virial, expected, max_relative = 1e-3);

                let ewald = SharedEwald::new(ewald);
                let energy = ewald.energy(&system);
                let virial = ewald.atomic_virial(&system).trace();
                assert_relative_eq!(energy, virial, max_relative = 1e-3);
            }

            #[test]
            fn nist4() {
                let system = get_system("spce-4.xyz");

                let alpha = 5.6 / system.cell.a();
                let mut ewald = Ewald::new(CUTOFF, 5, alpha);
                ewald.restriction = PairRestriction::InterMolecular;
                ewald.prepare(&system.cell);

                let energy = ewald.real_space_energy(&system) / K_BOLTZMANN;
                let expected = 1.057604e7;
                assert_relative_eq!(energy, expected, max_relative = 1e-4);

                let energy = ewald.k_space_energy(&system) / K_BOLTZMANN;
                let expected = 7.58785e3;
                assert_relative_eq!(energy, expected, max_relative = 2e-3);

                let energy = ewald.self_energy(&system) / K_BOLTZMANN;
                let expected = -1.42235e7;
                assert_relative_eq!(energy, expected, max_relative = 1e-4);
            }

            #[test]
            fn nist4_virial() {
                let system = get_system("spce-4.xyz");
                let mut ewald = Ewald::new(CUTOFF, 11, 0.333033);
                ewald.restriction = PairRestriction::InterMolecular;
                ewald.prepare(&system.cell);

                let convert = units::from(1.0, "atm").unwrap() * system.volume();

                let virial = ewald.real_space_atomic_virial(&system);
                let expected = convert * Matrix3::new([
                    [-5782.9498, -144.67187, -31.932323],
                    [-144.67187, -6052.8821, -121.95614],
                    [-31.932323, -121.95614, -5205.1646],
                ]);
                assert_relative_eq!(virial, expected, max_relative = 1e-3);

                let virial = ewald.k_space_atomic_virial(&system);
                let expected = convert * Matrix3::new([
                    [-288.2612,  -9.676199,  13.775332],
                    [-9.676199, -261.38081,  3.8066323],
                    [13.775332,  3.8066323, -322.02403],
                ]);
                assert_relative_eq!(virial, expected, max_relative = 1e-3);

                let ewald = SharedEwald::new(ewald);
                let energy = ewald.energy(&system);
                let virial = ewald.atomic_virial(&system).trace();
                assert_relative_eq!(energy, virial, max_relative = 1e-3);
            }
        }
    }
}
