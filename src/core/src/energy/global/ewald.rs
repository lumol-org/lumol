// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license
#![cfg_attr(rustfmt, rustfmt_skip)]

use std::ops::{Index, IndexMut, Deref, Range};
use std::sync::{RwLock, RwLockReadGuard, RwLockWriteGuard};
use std::f64::consts::{PI, FRAC_2_SQRT_PI};
use std::f64;

use ndarray::Zip;

use math::*;
use sys::{Configuration, UnitCell, CellShape};
use types::{Matrix3, Vector3D, Array3, Complex, Zero, One};
use consts::ELCC;
use energy::{PairRestriction, RestrictionInfo};
use parallel::ThreadLocalStore;
use parallel::prelude::*;

use super::{GlobalPotential, CoulombicPotential, GlobalCache};

/// 3D array with negative indexing on the first dimmension, for use in Ewald
/// phase factors.
///
/// # Examples
///
///  ```no-run
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
///
/// All the vectors contains the term corresponding to the k-vector indexes in
/// `self.kvecs`.
#[derive(Clone, Debug)]
struct EwaldFactors {
    /// Energetic pre-factor: `4 π / V exp(- k² / (4 α²)) / k²`
    energy: Vec<f64>,
    /// Electric field/force pre-factor: `8 π / V exp(- k² / (4 α²)) / k² \vec k / k`
    efield: Vec<Vector3D>,
    /// Virial pre-factor: `𝟙 - 2 (1 / k² + 1 / (4 α²)) \vec k ⊗ \vec k / k²`
    virial: Vec<Matrix3>,
    /// Indexes in k-space
    kvecs: Vec<(isize, isize, isize)>,
}

impl EwaldFactors {
    /// Create a new empty EwaldFactors
    pub fn new() -> EwaldFactors {
        EwaldFactors {
            energy: Vec::new(),
            efield: Vec::new(),
            virial: Vec::new(),
            kvecs: Vec::new(),
        }
    }

    /// Remove any data in the factors
    fn clear(&mut self) {
        self.energy.clear();
        self.efield.clear();
        self.virial.clear();
        self.kvecs.clear();
    }

    /// Reserve memory for at leats `size` items
    fn reserve(&mut self, size: usize) {
        self.energy.reserve(size);
        self.efield.reserve(size);
        self.virial.reserve(size);
        self.kvecs.reserve(size);
    }

    /// Compute the factors for the given `cell` and Ewald `parameters`
    pub fn compute(&mut self, cell: &UnitCell, parameters: &EwaldParameters) {
        self.clear();
        let kmax = parameters.kmax;
        let kmax3d = 4 * kmax * kmax * kmax + 6 * kmax * kmax + 3 * kmax;
        self.reserve(kmax3d as usize);

        match cell.shape() {
            CellShape::Infinite => panic!("Ewald is not defined with infinite unit cell"),
            CellShape::Orthorhombic => self.compute_ortho(cell, parameters),
            CellShape::Triclinic => self.compute_triclinic(cell, parameters),
        }
    }

    fn compute_ortho(&mut self, cell: &UnitCell, parameters: &EwaldParameters) {
        // TODO: there is a faster algorithm for orthorhombic cell
        self.compute_triclinic(cell, parameters);
    }

    fn compute_triclinic(&mut self, cell: &UnitCell, parameters: &EwaldParameters) {
        let alpha_sq_inv_fourth = 0.25 / (parameters.alpha * parameters.alpha);
        let four_pi_v = 4.0 * PI / cell.volume();

        let kmax = parameters.kmax;
        // k-vectors with a positive `ikx`
        for ikx in 1..kmax {
            for iky in -kmax..kmax {
                for ikz in -kmax..kmax {
                    let kvec = cell.k_vector([ikx as f64, iky as f64, ikz as f64]);
                    let k2 = kvec.norm2();
                    if k2 > parameters.kmax2 {
                        continue;
                    }

                    self.kvecs.push((ikx, iky, ikz));
                    let efact = four_pi_v * f64::exp(- k2 * alpha_sq_inv_fourth) / k2;
                    self.energy.push(efact);
                    self.efield.push(2.0 * efact * kvec);
                    let vfact = -2.0 * (1.0 / k2 + alpha_sq_inv_fourth);
                    let virial = Matrix3::one() + vfact * kvec.tensorial(&kvec);
                    self.virial.push(efact * virial);
                }
            }
        }

        // k-vectors with `ikx = 0`
        for iky in 1..kmax {
            for ikz in -kmax..kmax {
                let kvec = cell.k_vector([0.0, iky as f64, ikz as f64]);
                let k2 = kvec.norm2();
                if k2 > parameters.kmax2 {
                    continue;
                }

                self.kvecs.push((0, iky, ikz));
                let efact = four_pi_v * f64::exp(- k2 * alpha_sq_inv_fourth) / k2;
                self.energy.push(efact);
                self.efield.push(2.0 * efact * kvec);
                let vfact = -2.0 * (1.0 / k2 + alpha_sq_inv_fourth);
                let virial = Matrix3::one() + vfact * kvec.tensorial(&kvec);
                self.virial.push(efact * virial);
            }
        }

        // k-vectors with `ikx = 0` and `ikz = 0`
        for ikz in 1..kmax {
            let kvec = cell.k_vector([0.0, 0.0, ikz as f64]);
            let k2 = kvec.norm2();
            if k2 > parameters.kmax2 {
                continue;
            }

            self.kvecs.push((0, 0, ikz));
            let efact = four_pi_v * f64::exp(- k2 * alpha_sq_inv_fourth) / k2;
            self.energy.push(efact);
            self.efield.push(2.0 * efact * kvec);
            let vfact = -2.0 * (1.0 / k2 + alpha_sq_inv_fourth);
            let virial = Matrix3::one() + vfact * kvec.tensorial(&kvec);
            self.virial.push(efact * virial);
        }
    }
}

/// Ewald summation for coulombic interactions.
///
/// The Ewald summation is based on a separation of the coulombic potential `U`
/// in two parts, using the trivial identity: `U(x) = U(x) * (f(x) + 1) - U(x) *
/// f(x)` where `f` is the `erf` function. This leads to a separation of the
/// conditionally convergent coulombic sum into two absolutely convergent sums:
/// one in real space, and the other in Fourier or k-space. For more information
/// about this algorithm see [FS2002].
///
/// # Examples
///
/// ```
/// use lumol_core::energy::{Ewald, SharedEwald};
/// use lumol_core::units;
///
/// let ewald = SharedEwald::new(
///     Ewald::new(/* cutoff */ 12.0, /* kmax */ 7, /* alpha */ None)
/// );
///
/// use lumol_core::sys::System;
/// use lumol_core::sys::Particle;
/// use lumol_core::sys::UnitCell;
/// use lumol_core::types::Vector3D;
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
/// system.add_particle(na);
/// system.add_particle(cl);
///
/// // Use Ewald summation for electrostatic interactions
/// system.set_coulomb_potential(Box::new(ewald));
///
/// assert_eq!(system.potential_energy(), -0.0695110962119286);
/// ```
///
/// [FS2002] Frenkel, D. & Smith, B. Understanding molecular simulation. (Academic press, 2002).
pub struct Ewald {
    /// Various Ewald parameters
    parameters: EwaldParameters,
    /// Ewald pre-factors, only depending on the system unit cell
    factors: EwaldFactors,
    /// Restriction scheme
    restriction: PairRestriction,
    /// Cached phase factors (e^{i k r})
    eikr: Ewald3DArray,
    /// Fourier transform of the electrostatic density (\sum q_i e^{i k r})
    ///
    /// The vector contain the terms corresponding to the k-vectors in
    /// `self.factors.kvecs`
    rho: Vec<Complex>,
    /// Caching the allocation for electric field calculation
    ///
    /// This will contain the electric field at each atom
    efield: Vec<Vector3D>,
    /// Guard for cache invalidation of `self.factors`
    previous_cell: Option<UnitCell>,
    /// Update the cached quantities
    updater: Option<Box<Fn(&mut Ewald) + Sync + Send>>,
}

impl Clone for Ewald {
    fn clone(&self) -> Ewald {
        Ewald {
            parameters: self.parameters.clone(),
            factors: self.factors.clone(),
            restriction: self.restriction.clone(),
            eikr: self.eikr.clone(),
            rho: self.rho.clone(),
            efield: self.efield.clone(),
            previous_cell: self.previous_cell.clone(),
            updater: None,
        }
    }
}

// direct acces to parameters as `self.<xxx>`
impl Deref for Ewald {
    type Target = EwaldParameters;
    fn deref(&self) -> &EwaldParameters {
        &self.parameters
    }
}

impl Ewald {
    /// Create an Ewald summation using the given `cutoff` radius in real
    /// space, and `kmax` points in k-space (Fourier space). If `alpha` is None,
    /// then the default value of `3 * π / (4 * cutoff)` is used.
    pub fn new<I: Into<Option<f64>>>(cutoff: f64, kmax: usize, alpha: I) -> Ewald {
        let alpha = alpha.into().unwrap_or(3.0 * PI / (cutoff * 4.0));
        if cutoff < 0.0 {
            fatal_error!("the cutoff can not be negative in Ewald");
        } else if alpha < 0.0 {
            fatal_error!("alpha can not be negative in Ewald");
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
            factors: EwaldFactors::new(),
            eikr: Ewald3DArray::zeros((0..0, 0, 0)),
            rho: Vec::new(),
            efield: Vec::new(),
            previous_cell: None,
            updater: None,
        }
    }

    fn precompute(&mut self, cell: &UnitCell) {
        if let Some(ref prev_cell) = self.previous_cell {
            if cell == prev_cell {
                // Do not recompute
                return;
            }
        }
        self.previous_cell = Some(*cell);

        let max = cell.k_vector([1.0, 1.0, 1.0]).max() * self.parameters.kmax as f64;
        self.parameters.kmax2 = 1.0001 * max * max;

        if self.parameters.rc > cell.lengths().min() / 2.0 {
            warn!("The Ewald cutoff is too high for this unit cell, energy might be wrong.");
        }

        self.factors.compute(cell, &self.parameters);
    }
}

/// Real space part of the summation
impl Ewald {
    /// Get the real-space energy for one pair at distance `r` with charges `qi`
    /// and `qj` ; and with restriction information for this pair in `info`.
    #[inline]
    fn real_space_energy_pair(&self, info: RestrictionInfo, qi: f64, qj: f64, r: f64) -> f64 {
        if r > self.rc || info.excluded {
            return 0.0
        }
        assert_eq!(info.scaling, 1.0, "Scaling restriction scheme using Ewald are not implemented");
        return qi * qj * erfc(self.alpha * r) / r / ELCC;
    }

    /// Get the real-space force for one pair at distance `rij` with charges
    /// `qi` and `qj` ; and with restriction information for this pair in
    /// `info`.
    #[inline]
    fn real_space_force_pair(&self, info: RestrictionInfo, qi: f64, qj: f64, rij: &Vector3D) -> Vector3D {
        let r = rij.norm();
        if r > self.rc || info.excluded {
            return Vector3D::new(0.0, 0.0, 0.0)
        }
        assert_eq!(info.scaling, 1.0, "Scaling restriction scheme using Ewald are not implemented");
        let mut factor = erfc(self.alpha * r) / r;
        factor += self.alpha * FRAC_2_SQRT_PI * exp(-self.alpha * self.alpha * r * r);
        factor *= qi * qj / (r * r) / ELCC;
        return factor * rij;
    }

    /// Real space contribution to the energy
    fn real_space_energy(&self, configuration: &Configuration) -> f64 {
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

                let distance = configuration.bond_distance(i, j);
                let info = self.restriction.information(distance);

                let r = configuration.distance(i, j);
                local_energy += self.real_space_energy_pair(info, qi, qj, r);
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

        let thread_forces = ThreadLocalStore::new(|| vec![Vector3D::zero(); natoms]);

        (0..natoms).into_par_iter().for_each(|i| {
            // Get the thread local forces Vec
            let mut thread_forces = thread_forces.borrow_mut();

            let qi = charges[i];
            if qi == 0.0 {
                return;
            }

            for j in i + 1..natoms {
                let qj = charges[j];
                if qj == 0.0 {
                    continue;
                }

                let distance = configuration.bond_distance(i, j);
                let info = self.restriction.information(distance);

                let rij = configuration.nearest_image(i, j);
                let force = self.real_space_force_pair(info, qi, qj, &rij);
                thread_forces[i] += force;
                thread_forces[j] -= force;
            }
        });

        // reduce the thread local values
        thread_forces.sum_local_values(forces);
    }

    /// Real space contribution to the virial
    fn real_space_virial(&self, configuration: &Configuration) -> Matrix3 {
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

                let distance = configuration.bond_distance(i, j);
                let info = self.restriction.information(distance);

                let rij = configuration.nearest_image(i, j);
                let force = self.real_space_force_pair(info, qi, qj, &rij);
                local_virial += force.tensorial(&rij);
            }
            local_virial
        });
        return virials.sum();
    }

    fn real_space_move_particles_cost(&self, configuration: &Configuration, idxes: &[usize], newpos: &[Vector3D]) -> f64 {
        let charges = configuration.particles().charge;
        let positions = configuration.particles().position;

        let mut e_old = 0.0;
        let mut e_new = 0.0;

        // Iterate over all interactions between a moved particle and a
        // particle not moved
        for (idx, &i) in idxes.iter().enumerate() {
            let qi = charges[i];
            if qi == 0.0 {continue}
            for j in (0..configuration.size()).filter(|x| !idxes.contains(x)) {
                let qj = charges[j];
                if qi == 0.0 {continue}

                let r_old = configuration.distance(i, j);
                let r_new = configuration.cell.distance(&newpos[idx], &positions[j]);

                let distance = configuration.bond_distance(i, j);
                let info = self.restriction.information(distance);

                e_old += self.real_space_energy_pair(info, qi, qj, r_old);
                e_new += self.real_space_energy_pair(info, qi, qj, r_new);
            }
        }

        // Iterate over all interactions between two moved particles
        for (idx, &i) in idxes.iter().enumerate() {
            let qi = charges[i];
            if qi == 0.0 {continue}
            for (jdx, &j) in idxes.iter().enumerate().skip(i + 1) {
                let qj = charges[j];
                if qj == 0.0 {continue}

                let r_old = configuration.distance(i, j);
                let r_new = configuration.cell.distance(&newpos[idx], &newpos[jdx]);

                let distance = configuration.bond_distance(i, j);
                let info = self.restriction.information(distance);

                e_old += self.real_space_energy_pair(info, qi, qj, r_old);
                e_new += self.real_space_energy_pair(info, qi, qj, r_new);
            }
        }

        return e_new - e_old;
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
        return -self.alpha / sqrt(PI) * q2 / ELCC;
    }
}


/// k-space part of the summation
impl Ewald {
    /// Compute the Fourier transform of the electrostatic density
    fn eik_dot_r(&mut self, configuration: &Configuration) {
        let natoms = configuration.size();
        let range = -self.kmax..self.kmax;
        self.eikr.resize_if_different((range, 3, natoms));
        self.rho.clear();

        let positions = configuration.particles().position;
        let charges = configuration.particles().charge;

        // do the k = -1, 0, 1 cases first
        for spatial in 0..3 {
            let mut k_idx = [0.0, 0.0, 0.0];
            k_idx[spatial] = 1.0;
            let kvec = configuration.cell.k_vector(k_idx);
            for i in 0..natoms {
                self.eikr[(0, spatial, i)] = Complex::cartesian(1.0, 0.0);
                self.eikr[(1, spatial, i)] = Complex::polar(1.0, kvec * positions[i]);
                self.eikr[(-1, spatial, i)] = self.eikr[(1, spatial, i)].conj();
            }
        }

        // compute the other values of k by recursion
        for spatial in 0..3 {
            for k in 2..self.kmax {
                for i in 0..natoms {
                    self.eikr[(k, spatial, i)] = self.eikr[(k - 1, spatial, i)] * self.eikr[(1, spatial, i)];
                    self.eikr[(-k, spatial, i)] = self.eikr[(k, spatial, i)].conj();
                }
            }
        }

        for &(ikx, iky, ikz) in &self.factors.kvecs {
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
    fn kspace_energy(&mut self, configuration: &Configuration) -> f64 {
        self.eik_dot_r(configuration);
        let energy: f64 = Zip::from(&self.factors.energy)
            .and(&self.rho)
            .par_map(|(factor, rho)| factor * rho.norm2())
            .sum();
        return energy / ELCC;
    }

    /// k-space contribution to the forces
    fn kspace_forces(&mut self, configuration: &Configuration, forces: &mut [Vector3D]) {
        assert_eq!(forces.len(), configuration.size());
        self.eik_dot_r(configuration);

        let natoms = configuration.size();
        self.efield.clear();
        self.efield.resize(natoms, Vector3D::zero());

        let thread_efield = ThreadLocalStore::new(|| vec![Vector3D::zero(); natoms]);
        Zip::from(&self.factors.kvecs)
            .and(&self.factors.efield)
            .and(&self.rho)
            .par_apply(|&(ikx, iky, ikz), factor, rho| {
                let mut thread_efield = thread_efield.borrow_mut();
                for i in 0..natoms {
                    let eikr = self.eikr[(ikx, 0, i)] *
                               self.eikr[(iky, 1, i)] *
                               self.eikr[(ikz, 2, i)];
                    let partial = eikr * rho.conj();
                    thread_efield[i] += partial.imag() * factor;
                }
        });

        thread_efield.sum_local_values(&mut self.efield);

        let charges = configuration.particles().charge;
        for (force, &charge, field) in izip!(&mut *forces, charges, &self.efield) {
            *force += charge * field / ELCC;
        }
    }

    /// k-space contribution to the virial
    fn kspace_virial(&mut self, configuration: &Configuration) -> Matrix3 {
        self.eik_dot_r(configuration);

        let virial: Matrix3 = Zip::from(&self.factors.virial)
            .and(&self.rho)
            .par_map(|(factor, rho)| rho.norm2() * factor)
            .sum();

        return virial / ELCC;
    }

    /// Compute the Fourier transform of the electrostatic density changes
    /// while moving the particles in `idxes` from there curent position in
    /// the configuration to `newpos`;
    fn delta_rho_move_particles(&mut self, configuration: &Configuration, idxes: &[usize], newpos: &[Vector3D]) -> Vec<Complex> {
        let natoms = idxes.len();
        let mut new_eikr = Ewald3DArray::zeros((-self.kmax..self.kmax, 3, natoms));

        // Do the k=0, 1 cases first
        for spatial in 0..3 {
            let mut k_idx = [0.0, 0.0, 0.0];
            k_idx[spatial] = 1.0;
            let kvec = configuration.cell.k_vector(k_idx);
            for i in 0..natoms {
                new_eikr[(0, spatial, i)] = Complex::cartesian(1.0, 0.0);
                new_eikr[(1, spatial, i)] = Complex::polar(1.0, kvec * newpos[i]);
                new_eikr[(-1, spatial, i)] = new_eikr[(1, spatial, i)].conj();
            }
        }

        // Use recursive definition for computing the factor for all the other values of k.
        for spatial in 0..3 {
            for k in 2..self.kmax {
                for i in 0..natoms {
                    new_eikr[(k, spatial, i)] = new_eikr[(k - 1, spatial, i)] * new_eikr[(1, spatial, i)];
                    new_eikr[(-k, spatial, i)] = new_eikr[(k, spatial, i)].conj();
                }
            }
        }

        let mut delta = Vec::new();
        let charges = configuration.particles().charge;
        for &(ikx, iky, ikz) in &self.factors.kvecs {
            let mut partial = Complex::zero();
            for (idx, &i) in idxes.iter().enumerate() {
                let old_phi = self.eikr[(ikx, 0, i)] *
                              self.eikr[(iky, 1, i)] *
                              self.eikr[(ikz, 2, i)];

                let new_phi = new_eikr[(ikx, 0, idx)] *
                              new_eikr[(iky, 1, idx)] *
                              new_eikr[(ikz, 2, idx)];

                partial += charges[i] * (new_phi - old_phi);
            }
            delta.push(partial);
        }

        return delta;
    }

    fn kspace_move_particles_cost(&mut self, configuration: &Configuration, idxes: &[usize], newpos: &[Vector3D]) -> f64 {
        let mut e_old = 0.0;
        for (factor, &rho) in izip!(&self.factors.energy, &self.rho) {
            e_old += factor * rho.norm2();
        }
        e_old /= ELCC;

        let delta_rho = self.delta_rho_move_particles(configuration, idxes, newpos);

        let mut e_new = 0.0;
        for (factor, &rho, &delta) in izip!(&self.factors.energy, &self.rho, &delta_rho) {
            e_new += factor * (rho + delta).norm2();
        }
        e_new /= ELCC;

        self.updater = Some(Box::new(move |ewald: &mut Ewald| {
            for (rho, &delta) in izip!(&mut ewald.rho, &delta_rho) {
                *rho += delta;
            }
        }));

        return e_new - e_old;
    }
}

/// Molecular correction for Ewald summation
impl Ewald {
    /// Get the molecular correction energy for the pair with charges `qi` and
    /// `qj`, at distance `rij` and with restriction information in `info`.
    #[inline]
    fn molcorrect_energy_pair(&self, info: RestrictionInfo, qi: f64, qj: f64, r: f64) -> f64 {
        assert!(info.excluded, "Can not compute molecular correction for non-excluded pair");
        assert_eq!(info.scaling, 1.0, "Scaling restriction scheme using Ewald are not implemented");
        assert!(r < self.rc, "Atoms in molecule are separated by more than the cutoff radius of Ewald sum.");

        return - qi * qj / ELCC * erf(self.alpha * r) / r;
    }

    /// Get the molecular correction force for the pair with charges `qi` and
    /// `qj`, at distance `rij` and with restriction information in `info`.
    #[inline]
    fn molcorrect_force_pair(&self, info: RestrictionInfo, qi: f64, qj: f64, rij: &Vector3D) -> Vector3D {
        assert!(info.excluded, "Can not compute molecular correction for non-excluded pair");
        assert_eq!(info.scaling, 1.0, "Scaling restriction scheme using Ewald are not implemented");
        let r = rij.norm();
        assert!(r < self.rc, "Atoms in molecule are separated by more than the cutoff radius of Ewald sum.");

        let qiqj = qi * qj / (ELCC * r * r);
        let factor = qiqj * (2.0 * self.alpha / sqrt(PI) * exp(-self.alpha * self.alpha * r * r) - erf(self.alpha * r) / r);
        return factor * rij;
    }

    /// Molecular correction contribution to the energy
    fn molcorrect_energy(&self, configuration: &Configuration) -> f64 {
        let natoms = configuration.size();
        let charges = configuration.particles().charge;
        let mut energy = 0.0;

        for i in 0..natoms {
            let qi = charges[i];
            if qi == 0.0 {continue}
            // I can not manage to get this work with a loop from (i+1) to N. The finite
            // difference test (testing that the force is the same that the finite difference
            // of the energy) always fail. So let's use it that way for now.
            for j in i+1..natoms {
                // Only account for excluded pairs
                let distance = configuration.bond_distance(i, j);
                let info = self.restriction.information(distance);
                if !info.excluded {continue}

                let qj = charges[j];
                if qj == 0.0 {continue}

                let r = configuration.distance(i, j);
                energy += self.molcorrect_energy_pair(info, qi, qj, r);
            }
        }
        return energy;
    }

    /// Molecular correction contribution to the forces
    fn molcorrect_forces(&self, configuration: &Configuration, forces: &mut [Vector3D]) {
        let natoms = configuration.size();
        let charges = configuration.particles().charge;
        assert_eq!(forces.len(), natoms);

        for i in 0..natoms {
            let qi = charges[i];
            if qi == 0.0 {continue}
            for j in i+1..natoms {
                let distance = configuration.bond_distance(i, j);
                let info = self.restriction.information(distance);
                // Only account for excluded pairs
                if !info.excluded {continue}

                let qj = charges[j];
                if qj == 0.0 {continue}

                let rij = configuration.nearest_image(i, j);
                let force = self.molcorrect_force_pair(info, qi, qj, &rij);
                forces[i] += force;
                forces[j] -= force;
            }
        }
    }

    /// Molecular correction contribution to the virial
    fn molcorrect_virial(&self, configuration: &Configuration) -> Matrix3 {
        let natoms = configuration.size();
        let charges = configuration.particles().charge;
        let mut virial = Matrix3::zero();

        for i in 0..natoms {
            let qi = charges[i];
            if qi == 0.0 {continue}
            for j in i+1..natoms {
                let distance = configuration.bond_distance(i, j);
                let info = self.restriction.information(distance);
                // Only account for excluded pairs
                if !info.excluded {continue}

                let qj = charges[j];
                if qj == 0.0 {continue}

                let rij = configuration.nearest_image(i, j);
                let force = self.molcorrect_force_pair(info, qi, qj, &rij);
                virial += force.tensorial(&rij);
            }
        }
        return virial;
    }

    fn molcorrect_move_particles_cost(&mut self, configuration: &Configuration, idxes: &[usize], newpos: &[Vector3D]) -> f64 {
        let charges = configuration.particles().charge;
        let positions = configuration.particles().position;

        let mut e_old = 0.0;
        let mut e_new = 0.0;

        // Iterate over all interactions between a moved particle and a
        // particle not moved
        for (idx, &i) in idxes.iter().enumerate() {
            let qi = charges[i];
            if qi == 0.0 {continue}
            for j in (0..configuration.size()).filter(|x| !idxes.contains(x)) {
                let qj = charges[j];
                if qi == 0.0 {continue}

                let distance = configuration.bond_distance(i, j);
                let info = self.restriction.information(distance);
                if !info.excluded {continue}

                let r_old = configuration.distance(i, j);
                let r_new = configuration.cell.distance(&newpos[idx], &positions[j]);

                e_old += self.molcorrect_energy_pair(info, qi, qj, r_old);
                e_new += self.molcorrect_energy_pair(info, qi, qj, r_new);
            }
        }

        // Iterate over all interactions between two moved particles
        for (idx, &i) in idxes.iter().enumerate() {
            let qi = charges[i];
            if qi == 0.0 {continue}
            for (jdx, &j) in idxes.iter().enumerate().skip(i + 1) {
                let qj = charges[j];
                if qj == 0.0 {continue}

                let distance = configuration.bond_distance(i, j);
                let info = self.restriction.information(distance);
                if !info.excluded {continue}

                let r_old = configuration.distance(i, j);
                let r_new = configuration.cell.distance(&newpos[idx], &newpos[jdx]);

                e_old += self.molcorrect_energy_pair(info, qi, qj, r_old);
                e_new += self.molcorrect_energy_pair(info, qi, qj, r_new);
            }
        }

        return e_new - e_old;
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
    /// let boxed: Box<CoulombicPotential> = Box::new(ewald);
    /// ```
    pub fn new(ewald: Ewald) -> SharedEwald {
        SharedEwald(RwLock::new(ewald))
    }

    /// Get read access to the underlying Ewald solver
    fn read(&self) -> RwLockReadGuard<Ewald> {
        // The lock should never be poisonned, because any panic will unwind
        // and finish the simulation.
        self.0.read().expect("Ewald lock is poisonned")
    }

    /// Get write access to the underlying Ewald solver
    fn write(&self) -> RwLockWriteGuard<Ewald> {
        // The lock should never be poisonned, because any panic will unwind
        // and finish the simulation.
        self.0.write().expect("Ewald lock is poisonned")
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
        ewald.precompute(&configuration.cell);
        let real = ewald.real_space_energy(configuration);
        let self_e = ewald.self_energy(configuration);
        let kspace = ewald.kspace_energy(configuration);
        let molecular = ewald.molcorrect_energy(configuration);
        return real + self_e + kspace + molecular;
    }

    fn forces(&self, configuration: &Configuration, forces: &mut [Vector3D])  {
        assert_eq!(forces.len(), configuration.size());
        let mut ewald = self.write();
        ewald.precompute(&configuration.cell);

        ewald.real_space_forces(configuration, forces);
        // No self force
        ewald.kspace_forces(configuration, forces);
        ewald.molcorrect_forces(configuration, forces);
    }

    fn virial(&self, configuration: &Configuration) -> Matrix3 {
        let mut ewald = self.write();
        ewald.precompute(&configuration.cell);
        let real = ewald.real_space_virial(configuration);
        // No self virial
        let kspace = ewald.kspace_virial(configuration);
        let molecular = ewald.molcorrect_virial(configuration);
        return real + kspace + molecular;
    }
}

impl CoulombicPotential for SharedEwald {
    fn set_restriction(&mut self, restriction: PairRestriction) {
        self.write().restriction = restriction;
    }
}

impl GlobalCache for SharedEwald {
    fn move_particles_cost(&self, configuration: &Configuration, idxes: &[usize], newpos: &[Vector3D]) -> f64 {
        let mut ewald = self.write();
        ewald.precompute(&configuration.cell);
        let real = ewald.real_space_move_particles_cost(configuration, idxes, newpos);
        /* No self cost */
        let kspace = ewald.kspace_move_particles_cost(configuration, idxes, newpos);
        let molecular = ewald.molcorrect_move_particles_cost(configuration, idxes, newpos);
        return real + kspace + molecular;
    }

    fn update(&self) {
        use std::mem;

        let mut ewald = self.write();
        if ewald.updater.is_some() {
            let mut updater = None;
            mem::swap(&mut updater, &mut ewald.updater);
            let updater = updater.unwrap();
            updater(&mut *ewald);
        }
    }
}

#[cfg(test)]
mod tests {
    pub use super::*;
    use sys::System;
    use utils::system_from_xyz;

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
        use utils::system_from_xyz;
        let mut system = system_from_xyz("3
        bonds cell: 20.0
        O  0.0  0.0  0.0
        H -0.7 -0.7  0.3
        H  0.3 -0.3 -0.8
        ");
        assert!(system.molecules().len() == 1);

        for particle in system.particles_mut() {
            if particle.name == "O" {
                *particle.charge = -0.8476;
            } else if particle.name == "H" {
                *particle.charge = 0.4238;
            }
        }
        return system;
    }

    mod errors {
        use super::*;
        use energy::GlobalPotential;
        use sys::UnitCell;

        #[test]
        #[should_panic]
        fn infinite_cell() {
            let mut system = nacl_pair();
            system.cell = UnitCell::infinite();
            let ewald = SharedEwald::new(Ewald::new(8.0, 10, None));
            let _ = ewald.energy(&system);
        }

        #[test]
        #[should_panic]
        fn negative_cutoff() {
            let _ = Ewald::new(-8.0, 10, None);
        }

        #[test]
        #[should_panic]
        fn negative_alpha() {
            let _ = Ewald::new(8.0, 10, -45.2);
        }
    }

    mod pairs {
        use super::*;
        use energy::GlobalPotential;

        #[test]
        fn energy() {
            let system = nacl_pair();
            let ewald = SharedEwald::new(Ewald::new(8.0, 10, None));

            let energy = ewald.energy(&system);
            // This was computed by hand
            let energy_brute_force = -0.09262397663346732;
            assert_ulps_eq!(energy, energy_brute_force, epsilon=1e-4);
        }

        #[test]
        fn real_forces_finite_differences() {
            let mut system = nacl_pair();
            let mut ewald = Ewald::new(8.0, 10, None);
            ewald.precompute(&system.cell);

            let e = ewald.real_space_energy(&system);
            let eps = 1e-9;
            system.particles_mut().position[0][0] += eps;

            let e1 = ewald.real_space_energy(&system);
            let mut forces = vec![Vector3D::zero(); 2];
            ewald.real_space_forces(&system, &mut forces);
            assert_relative_eq!((e - e1) / eps, forces[0][0], epsilon=1e-6);
        }

        #[test]
        fn kspace_forces_finite_differences() {
            let mut system = nacl_pair();
            // Using a small cutoff to increase the weight of k-space
            let mut ewald = Ewald::new(2.0, 10, None);
            ewald.precompute(&system.cell);

            let e = ewald.kspace_energy(&system);
            let eps = 1e-9;
            system.particles_mut().position[0][0] += eps;

            let e1 = ewald.kspace_energy(&system);
            let mut forces = vec![Vector3D::zero(); 2];
            ewald.kspace_forces(&system, &mut forces);
            assert_relative_eq!((e - e1) / eps, forces[0][0], epsilon=1e-6);
        }

        #[test]
        fn molcorrect_forces_finite_differences() {
            let mut system = nacl_pair();
            let mut ewald = Ewald::new(8.0, 10, None);
            ewald.precompute(&system.cell);

            let e = ewald.molcorrect_energy(&system);
            let eps = 1e-9;
            system.particles_mut().position[0][0] += eps;

            let e1 = ewald.molcorrect_energy(&system);
            let mut forces = vec![Vector3D::zero(); 2];
            ewald.molcorrect_forces(&system, &mut forces);
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
        use types::{Vector3D, Zero};
        use energy::{GlobalPotential, PairRestriction, CoulombicPotential};

        #[test]
        fn energy() {
            let system = water();
            let mut ewald = SharedEwald::new(Ewald::new(8.0, 10, None));
            ewald.set_restriction(PairRestriction::InterMolecular);

            let energy = ewald.energy(&system);
            let expected = -0.000009243358925787454;
            assert_ulps_eq!(energy, expected);

            let molcorrect = ewald.read().molcorrect_energy(&system);
            let expected = 0.02452968743897957;
            assert_ulps_eq!(molcorrect, expected);
        }

        #[test]
        fn real_space_forces_finite_differences() {
            let mut system = water();
            let mut ewald = Ewald::new(8.0, 10, None);
            ewald.restriction = PairRestriction::InterMolecular;
            ewald.precompute(&system.cell);

            let mut forces = vec![Vector3D::zero(); 3];
            ewald.real_space_forces(&system, &mut forces);
            let force = forces[0][0];

            let eps = 1e-9;
            let e = ewald.real_space_energy(&system);
            system.particles_mut().position[0][0] += eps;
            let e1 = ewald.real_space_energy(&system);

            // No real space energetic contribution here, we only have one
            // molecule.
            assert_ulps_eq!(e, 0.0);
            assert_ulps_eq!(e1, 0.0);
            assert_ulps_eq!(force, 0.0);
        }

        #[test]
        fn kspace_forces_finite_differences() {
            let mut system = water();
            let mut ewald = Ewald::new(8.0, 10, None);
            ewald.restriction = PairRestriction::InterMolecular;
            ewald.precompute(&system.cell);

            let mut forces = vec![Vector3D::zero(); 3];
            ewald.kspace_forces(&system, &mut forces);
            let force = forces[0][0];

            let eps = 1e-9;
            let e = ewald.kspace_energy(&system);
            system.particles_mut().position[0][0] += eps;
            let e1 = ewald.kspace_energy(&system);

            assert_relative_eq!((e - e1) / eps, force, epsilon = 1e-6);
        }

        #[test]
        fn molcorrect_forces_finite_differences() {
            let mut system = water();
            let mut ewald = Ewald::new(8.0, 10, None);
            ewald.restriction = PairRestriction::InterMolecular;
            ewald.precompute(&system.cell);

            let mut forces = vec![Vector3D::zero(); 3];
            ewald.molcorrect_forces(&system, &mut forces);
            let force = forces[0][0];

            let eps = 1e-9;
            let e = ewald.molcorrect_energy(&system);
            system.particles_mut().position[0][0] += eps;
            let e1 = ewald.molcorrect_energy(&system);

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

    mod virial {
        use super::*;
        use types::{Matrix3, Zero, One};
        use energy::PairRestriction;

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

        #[test]
        fn virial_is_energy() {
            // A nice property of Ewald summation for point charge systems
            let system = water();
            let ewald = SharedEwald::new(Ewald::new(8.0, 10, None));

            let energy = ewald.energy(&system);
            let virial = ewald.virial(&system).trace();
            assert_relative_eq!(energy, virial, max_relative = 1e-3);

            let system = nacl_pair();
            let ewald = SharedEwald::new(Ewald::new(8.0, 10, None));

            let energy = ewald.energy(&system);
            let virial = ewald.virial(&system).trace();
            assert_relative_eq!(energy, virial, max_relative = 1e-3);
        }

        #[test]
        fn real_space_finite_differences() {
            let mut system = water();
            let mut ewald = Ewald::new(8.0, 10, None);
            ewald.restriction = PairRestriction::InterMolecular;
            ewald.precompute(&system.cell);

            let eps = 1e-9;
            let virial = ewald.real_space_virial(&system);
            let mut finite_diff = Matrix3::zero();

            for i in 0..3 {
                for j in 0..3 {
                    ewald.precompute(&system.cell);
                    let e = ewald.real_space_energy(&system);
                    scale(&mut system, i, j, eps);
                    ewald.precompute(&system.cell);
                    let e1 = ewald.real_space_energy(&system);
                    finite_diff[i][j] = (e - e1) / eps;
                }
            }

            assert_relative_eq!(virial, finite_diff, epsilon = 1e-6);
        }

        #[test]
        fn kspace_finite_differences() {
            let mut system = water();
            let mut ewald = Ewald::new(2.0, 10, None);
            ewald.restriction = PairRestriction::InterMolecular;
            ewald.precompute(&system.cell);

            let eps = 1e-9;
            let virial = ewald.kspace_virial(&system);
            let mut finite_diff = Matrix3::zero();

            for i in 0..3 {
                for j in 0..3 {
                    ewald.precompute(&system.cell);
                    let e = ewald.kspace_energy(&system);
                    scale(&mut system, i, j, eps);
                    ewald.precompute(&system.cell);
                    let e1 = ewald.kspace_energy(&system);
                    finite_diff[i][j] = (e - e1) / eps;
                }
            }

            // Make sure the finite_diff matrix is symetric
            finite_diff = (finite_diff + finite_diff.transposed()) / 2.0;
            assert_relative_eq!(virial, finite_diff, epsilon = 1e-6);
        }

        #[test]
        fn molcorrect_finite_differences() {
            let mut system = water();
            let mut ewald = Ewald::new(2.0, 10, None);
            ewald.restriction = PairRestriction::InterMolecular;
            ewald.precompute(&system.cell);

            let eps = 1e-9;
            let virial = ewald.molcorrect_virial(&system);
            let mut finite_diff = Matrix3::zero();

            for i in 0..3 {
                for j in 0..3 {
                    ewald.precompute(&system.cell);
                    let e = ewald.molcorrect_energy(&system);
                    scale(&mut system, i, j, eps);
                    ewald.precompute(&system.cell);
                    let e1 = ewald.molcorrect_energy(&system);
                    finite_diff[i][j] = (e - e1) / eps;
                }
            }

            assert_relative_eq!(virial, finite_diff, epsilon = 1e-6);
        }
    }

    mod cache {
        use super::*;
        use sys::System;
        use types::Vector3D;
        use energy::{GlobalPotential, PairRestriction, CoulombicPotential, GlobalCache};

        pub fn testing_system() -> System {
            use utils::system_from_xyz;
            let mut system = system_from_xyz("6
            bonds cell: 20.0
            O  0.0  0.0  0.0
            H -0.7 -0.7  0.3
            H  0.3 -0.3 -0.8
            O  2.0  2.0  0.0
            H  1.3  1.3  0.3
            H  2.3  1.7 -0.8
            ");
            assert!(system.molecules().len() == 2);

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
        fn move_atoms() {
            let mut system = testing_system();
            let mut ewald = SharedEwald::new(Ewald::new(8.0, 10, None));
            ewald.set_restriction(PairRestriction::InterMolecular);
            ewald.write().precompute(&system.cell);

            let check = ewald.clone();
            // Initialize cached values
            let _ = ewald.energy(&system);

            let old_e = check.energy(&system);
            let idxes = &[0, 1];
            let newpos = &[Vector3D::new(0.0, 0.0, 0.5), Vector3D::new(-0.7, 0.2, 1.5)];

            let cost = ewald.move_particles_cost(&system, idxes, newpos);

            system.particles_mut().position[0] = newpos[0];
            system.particles_mut().position[1] = newpos[1];
            let new_e = check.energy(&system);
            assert_ulps_eq!(cost, new_e - old_e);
        }

        #[test]
        fn move_atoms_real_space() {
            let mut system = testing_system();
            let mut ewald = Ewald::new(8.0, 10, None);
            ewald.restriction = PairRestriction::InterMolecular;
            ewald.precompute(&system.cell);

            let check = ewald.clone();
            // Initialize cached values
            let _ = ewald.real_space_energy(&system);

            let old_e = check.real_space_energy(&system);
            let idxes = &[0, 1];
            let newpos = &[Vector3D::new(0.0, 0.0, 0.5), Vector3D::new(-0.7, 0.2, 1.5)];

            let cost = ewald.real_space_move_particles_cost(&system, idxes, newpos);

            system.particles_mut().position[0] = newpos[0];
            system.particles_mut().position[1] = newpos[1];
            let new_e = check.real_space_energy(&system);
            assert_ulps_eq!(cost, new_e - old_e);
        }

        #[test]
        fn move_atoms_kspace() {
            let mut system = testing_system();
            let mut ewald = Ewald::new(2.0, 10, None);
            ewald.restriction = PairRestriction::InterMolecular;
            ewald.precompute(&system.cell);

            let mut check = ewald.clone();
            // Initialize cached values
            let _ = ewald.kspace_energy(&system);

            let old_e = check.kspace_energy(&system);
            let idxes = &[0, 1];
            let newpos = &[Vector3D::new(0.0, 0.0, 0.5), Vector3D::new(-0.7, 0.2, 1.5)];

            let cost = ewald.kspace_move_particles_cost(&system, idxes, newpos);

            system.particles_mut().position[0] = newpos[0];
            system.particles_mut().position[1] = newpos[1];
            let new_e = check.kspace_energy(&system);
            assert_ulps_eq!(cost, new_e - old_e);
        }

        #[test]
        fn move_atoms_molcorrect() {
            let mut system = testing_system();
            let mut ewald = Ewald::new(8.0, 10, None);
            ewald.restriction = PairRestriction::InterMolecular;
            ewald.precompute(&system.cell);

            let check = ewald.clone();
            // Initialize cached values
            let _ = ewald.molcorrect_energy(&system);

            let old_e = check.molcorrect_energy(&system);
            let idxes = &[0, 1];
            let newpos = &[Vector3D::new(0.0, 0.0, 0.5), Vector3D::new(-0.7, 0.2, 1.5)];

            let cost = ewald.molcorrect_move_particles_cost(&system, idxes, newpos);

            system.particles_mut().position[0] = newpos[0];
            system.particles_mut().position[1] = newpos[1];
            let new_e = check.molcorrect_energy(&system);
            assert_ulps_eq!(cost, new_e - old_e);
        }
    }

    // Comparing the value for each component of Ewald energy with the NIST
    // reference. See `tests/nist-spce.rs` for more information. These tests
    // check values that are not accessible from the outside of lumol-core.
    mod nist {
        use super::*;

        use std::path::Path;
        use std::fs::File;
        use std::io::Read;

        pub fn get_system(path: &str) -> System {
            let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
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
                    other => panic!("Unknown particle name: {}", other),
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

        mod cutoff_9 {
            use super::*;
            use consts::K_BOLTZMANN;
            const CUTOFF: f64 = 9.0;

            #[test]
            fn nist1() {
                let system = get_system("spce-1.xyz");

                let alpha = 5.6 / system.cell.a();
                let mut ewald = Ewald::new(CUTOFF, 5, alpha);
                ewald.restriction = PairRestriction::InterMolecular;
                ewald.precompute(&system.cell);

                let energy = ewald.real_space_energy(&system) / K_BOLTZMANN;
                let expected = -5.58904e5;
                assert_relative_eq!(energy, expected, max_relative = 1e-4);

                let energy = ewald.kspace_energy(&system) / K_BOLTZMANN;
                let expected = 6.27009e3;
                assert_relative_eq!(energy, expected, max_relative = 2e-3);

                let energy = ewald.self_energy(&system) / K_BOLTZMANN;
                let expected = -2.84469e6;
                assert_relative_eq!(energy, expected, max_relative = 1e-4);

                let energy = ewald.molcorrect_energy(&system) / K_BOLTZMANN;
                let expected = 2.80999e6;
                assert_relative_eq!(energy, expected, max_relative = 1e-4);
            }

            #[test]
            fn nist2() {
                let system = get_system("spce-2.xyz");

                let alpha = 5.6 / system.cell.a();
                let mut ewald = Ewald::new(CUTOFF, 5, alpha);
                ewald.restriction = PairRestriction::InterMolecular;
                ewald.precompute(&system.cell);

                let energy = ewald.real_space_energy(&system) / K_BOLTZMANN;
                let expected = -1.19308e6;
                assert_relative_eq!(energy, expected, max_relative = 1e-4);

                let energy = ewald.kspace_energy(&system) / K_BOLTZMANN;
                let expected = 6.03495e3;
                assert_relative_eq!(energy, expected, max_relative = 5e-3);

                let energy = ewald.self_energy(&system) / K_BOLTZMANN;
                let expected = -5.68938e6;
                assert_relative_eq!(energy, expected, max_relative = 1e-4);

                let energy = ewald.molcorrect_energy(&system) / K_BOLTZMANN;
                let expected = 5.61998e6;
                assert_relative_eq!(energy, expected, max_relative = 1e-4);
            }

            #[test]
            fn nist3() {
                let system = get_system("spce-3.xyz");

                let alpha = 5.6 / system.cell.a();
                let mut ewald = Ewald::new(CUTOFF, 5, alpha);
                ewald.restriction = PairRestriction::InterMolecular;
                ewald.precompute(&system.cell);

                let energy = ewald.real_space_energy(&system) / K_BOLTZMANN;
                let expected = -1.96320e6;
                assert_relative_eq!(energy, expected, max_relative = 1e-4);

                let energy = ewald.kspace_energy(&system) / K_BOLTZMANN;
                let expected = 5.24461e3;
                assert_relative_eq!(energy, expected, max_relative = 2e-3);

                let energy = ewald.self_energy(&system) / K_BOLTZMANN;
                let expected = -8.53407e6;
                assert_relative_eq!(energy, expected, max_relative = 1e-4);

                let energy = ewald.molcorrect_energy(&system) / K_BOLTZMANN;
                let expected = 8.42998e6;
                assert_relative_eq!(energy, expected, max_relative = 1e-4);
            }

            #[test]
            fn nist4() {
                let system = get_system("spce-4.xyz");

                let alpha = 5.6 / system.cell.a();
                let mut ewald = Ewald::new(CUTOFF, 5, alpha);
                ewald.restriction = PairRestriction::InterMolecular;
                ewald.precompute(&system.cell);

                let energy = ewald.real_space_energy(&system) / K_BOLTZMANN;
                let expected = -3.44720e6;
                assert_relative_eq!(energy, expected, max_relative = 1e-4);

                let energy = ewald.kspace_energy(&system) / K_BOLTZMANN;
                let expected = 7.58785e3;
                assert_relative_eq!(energy, expected, max_relative = 2e-3);

                let energy = ewald.self_energy(&system) / K_BOLTZMANN;
                let expected = -1.42235e7;
                assert_relative_eq!(energy, expected, max_relative = 1e-4);

                let energy = ewald.molcorrect_energy(&system) / K_BOLTZMANN;
                let expected = 1.41483e7;
                assert_relative_eq!(energy, expected, max_relative = 1e-4);
            }
        }

        mod cutoff_10 {
            use super::*;
            use consts::K_BOLTZMANN;
            const CUTOFF: f64 = 10.0;

            #[test]
            fn nist1() {
                let system = get_system("spce-1.xyz");

                let alpha = 5.6 / system.cell.a();
                let mut ewald = Ewald::new(CUTOFF, 5, alpha);
                ewald.restriction = PairRestriction::InterMolecular;
                ewald.precompute(&system.cell);

                let energy = ewald.real_space_energy(&system) / K_BOLTZMANN;
                let expected = -5.58889e5;
                assert_relative_eq!(energy, expected, max_relative = 1e-4);

                let energy = ewald.kspace_energy(&system) / K_BOLTZMANN;
                let expected = 6.27009e3;
                assert_relative_eq!(energy, expected, max_relative = 2e-3);

                let energy = ewald.self_energy(&system) / K_BOLTZMANN;
                let expected = -2.84469e6;
                assert_relative_eq!(energy, expected, max_relative = 1e-4);

                let energy = ewald.molcorrect_energy(&system) / K_BOLTZMANN;
                let expected = 2.80999e6;
                assert_relative_eq!(energy, expected, max_relative = 1e-4);
            }

            #[test]
            fn nist2() {
                let system = get_system("spce-2.xyz");

                let alpha = 5.6 / system.cell.a();
                let mut ewald = Ewald::new(CUTOFF, 5, alpha);
                ewald.restriction = PairRestriction::InterMolecular;
                ewald.precompute(&system.cell);

                let energy = ewald.real_space_energy(&system) / K_BOLTZMANN;
                let expected = -1.19295e6;
                assert_relative_eq!(energy, expected, max_relative = 1e-4);

                let energy = ewald.kspace_energy(&system) / K_BOLTZMANN;
                let expected = 6.03495e3;
                assert_relative_eq!(energy, expected, max_relative = 5e-3);

                let energy = ewald.self_energy(&system) / K_BOLTZMANN;
                let expected = -5.68938e6;
                assert_relative_eq!(energy, expected, max_relative = 1e-4);

                let energy = ewald.molcorrect_energy(&system) / K_BOLTZMANN;
                let expected = 5.61998e6;
                assert_relative_eq!(energy, expected, max_relative = 1e-4);
            }

            #[test]
            fn nist3() {
                let system = get_system("spce-3.xyz");

                let alpha = 5.6 / system.cell.a();
                let mut ewald = Ewald::new(CUTOFF, 5, alpha);
                ewald.restriction = PairRestriction::InterMolecular;
                ewald.precompute(&system.cell);

                let energy = ewald.real_space_energy(&system) / K_BOLTZMANN;
                let expected = -1.96297e6;
                assert_relative_eq!(energy, expected, max_relative = 1e-4);

                let energy = ewald.kspace_energy(&system) / K_BOLTZMANN;
                let expected = 5.24461e3;
                assert_relative_eq!(energy, expected, max_relative = 2e-3);

                let energy = ewald.self_energy(&system) / K_BOLTZMANN;
                let expected = -8.53407e6;
                assert_relative_eq!(energy, expected, max_relative = 1e-4);

                let energy = ewald.molcorrect_energy(&system) / K_BOLTZMANN;
                let expected = 8.42998e6;
                assert_relative_eq!(energy, expected, max_relative = 1e-4);
            }

            #[test]
            fn nist4() {
                let system = get_system("spce-4.xyz");

                let alpha = 5.6 / system.cell.a();
                let mut ewald = Ewald::new(CUTOFF, 5, alpha);
                ewald.restriction = PairRestriction::InterMolecular;
                ewald.precompute(&system.cell);

                let energy = ewald.real_space_energy(&system) / K_BOLTZMANN;
                let expected = -3.57226e6;
                assert_relative_eq!(energy, expected, max_relative = 1e-4);

                let energy = ewald.kspace_energy(&system) / K_BOLTZMANN;
                let expected = 7.58785e3;
                assert_relative_eq!(energy, expected, max_relative = 2e-3);

                let energy = ewald.self_energy(&system) / K_BOLTZMANN;
                let expected = -1.42235e7;
                assert_relative_eq!(energy, expected, max_relative = 1e-4);

                let energy = ewald.molcorrect_energy(&system) / K_BOLTZMANN;
                let expected = 1.41483e7;
                assert_relative_eq!(energy, expected, max_relative = 1e-4);
            }
        }
    }
}
