// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license
#![cfg_attr(rustfmt, rustfmt_skip)]

use std::sync::{RwLock, RwLockReadGuard, RwLockWriteGuard};
use std::f64::consts::{PI, FRAC_2_SQRT_PI};
use std::f64;

use ndarray::Zip;

use math::*;
use sys::{Configuration, UnitCell, CellShape};
use types::{Matrix3, Vector3D, Array3, Complex, Zero};
use consts::ELCC;
use energy::{PairRestriction, RestrictionInfo};
use parallel::ThreadLocalStore;
use parallel::prelude::*;

use super::{GlobalPotential, CoulombicPotential, GlobalCache};

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
///     Ewald::new(/* cutoff */ 12.0, /* kmax */ 7)
/// );
///
/// use lumol_core::sys::System;
/// use lumol_core::sys::Particle;
/// use lumol_core::sys::UnitCell;
/// use lumol_core::types::Vector3D;
///
/// // Setup a system containing a NaCl pair
/// let mut system = System::with_cell(UnitCell::cubic(10.0));
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
/// assert_eq!(system.potential_energy(), -0.07042996180522723);
/// ```
///
/// [FS2002] Frenkel, D. & Smith, B. Understanding molecular simulation. (Academic press, 2002).
#[derive(Clone, Debug)]
pub struct Ewald {
    /// Splitting parameter between k-space and real space
    alpha: f64,
    /// Cutoff radius in real space
    rc: f64,
    /// Number of points to use in k-space
    kmax: usize,
    /// Cutoff in k-space
    kmax2: f64,
    /// Restriction scheme
    restriction: PairRestriction,
    /// Caching exponential factors exp(-k^2 / (4 alpha^2)) / k^2
    expfactors: Array3<f64>,
    /// Phases for the Fourier transform, cached allocation
    fourier_phases: Array3<Complex>,
    /// Fourier transform of the electrostatic density
    rho: Array3<Complex>,
    /// Fourier transform of the electrostatic density modifications, cached
    /// allocation and for updating `self.rho`
    delta_rho: Array3<Complex>,
    /// Guard for cache invalidation of `expfactors`
    previous_cell: Option<UnitCell>,
}

impl Ewald {
    /// Create an Ewald summation using the given `cutoff` radius in real
    /// space, and `kmax` points in k-space (Fourier space).
    pub fn new(cutoff: f64, kmax: usize) -> Ewald {
        let expfactors = Array3::zeros((kmax, kmax, kmax));
        let rho = Array3::zeros((kmax, kmax, kmax));
        Ewald {
            alpha: 3.0 * PI / (cutoff * 4.0),
            rc: cutoff,
            kmax: kmax,
            kmax2: 0.0,
            restriction: PairRestriction::None,
            expfactors: expfactors,
            fourier_phases: Array3::zeros((0, 0, 0)),
            rho: rho.clone(),
            delta_rho: rho,
            previous_cell: None,
        }
    }

    /// Set the value of the alpha parameter for ewald computation. The default is to use
    /// `alpha = 3 * π / (4 * rc)`.
    pub fn set_alpha(&mut self, alpha: f64) {
        assert!(alpha > 0.0, "Ewald parameter alpha must be positive");
        self.alpha = alpha;
    }

    fn precompute(&mut self, cell: &UnitCell) {
        if let Some(ref prev_cell) = self.previous_cell {
            if cell == prev_cell {
                // Do not recompute
                return;
            }
        }
        match cell.shape() {
            CellShape::Infinite => {
                fatal_error!("Can not use Ewald sum with Infinite cell.");
            },
            CellShape::Triclinic => {
                fatal_error!("Can not (yet) use Ewald sum with Triclinic cell.");
            },
            CellShape::Orthorhombic => {
                // All good!
            },
        }
        self.previous_cell = Some(*cell);

        // Because we do a spherical truncation in k space, we have to transform
        // kmax into a spherical cutoff 'radius'
        let lenghts = cell.lengths();
        let max_lenght = f64::max(f64::max(lenghts[0], lenghts[1]), lenghts[2]);
        let min_lenght = f64::min(f64::min(lenghts[0], lenghts[1]), lenghts[2]);
        let k_rc = self.kmax as f64 * (2.0 * PI / max_lenght);
        self.kmax2 = k_rc * k_rc;

        if self.rc > min_lenght / 2.0 {
            warn!("The Ewald cutoff is too high for this unit cell, energy might be wrong.");
        }

        // Now, we precompute the exp(-k^2 / (4 a^2)) / k^2 terms. We use the
        // symmetry to only store (ikx >= 0 && iky >= 0  && ikz >= 0 ) terms
        let (rec_vx, rec_vy, rec_vz) = cell.reciprocal_vectors();
        for ikx in 0..self.kmax {
            let kx = (ikx as f64) * rec_vx;
            for iky in 0..self.kmax {
                let ky = kx + (iky as f64) * rec_vy;
                for ikz in 0..self.kmax {
                    let k = ky + (ikz as f64) * rec_vz;
                    let k2 = k.norm2();
                    if k2 > self.kmax2 {
                        self.expfactors[(ikx, iky, ikz)] = 0.0;
                        continue;
                    }
                    self.expfactors[(ikx, iky, ikz)] = exp(-k2 / (4.0 * self.alpha * self.alpha)) / k2;
                    if ikx != 0 {self.expfactors[(ikx, iky, ikz)] *= 2.0;}
                    if iky != 0 {self.expfactors[(ikx, iky, ikz)] *= 2.0;}
                    if ikz != 0 {self.expfactors[(ikx, iky, ikz)] *= 2.0;}
                }
            }
        }
        self.expfactors[(0, 0, 0)] = 0.0;
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
        let mut factor = erfc(self.alpha * r) / r;
        factor += self.alpha * FRAC_2_SQRT_PI * exp(-self.alpha * self.alpha * r * r);
        factor *= qi * qj / (r * r) / ELCC;
        return factor * rij;
    }

    /// Real space contribution to the energy
    fn real_space_energy(&self, configuration: &Configuration) -> f64 {
        let natoms = configuration.size();
        let charges = configuration.particles().charge;

        let mut energy = 0.0;
        for i in 0..natoms {
            let qi = charges[i];
            if qi == 0.0 {continue}
            for j in i+1..natoms {
                let qj = charges[j];
                if qj == 0.0 {continue}

                let distance = configuration.bond_distance(i, j);
                let info = self.restriction.information(distance);

                let r = configuration.distance(i, j);
                energy += self.real_space_energy_pair(info, qi, qj, r);
            }
        }
        return energy;
    }

    /// Real space contribution to the forces
    fn real_space_forces(&self, configuration: &Configuration, forces: &mut [Vector3D]) {
        assert_eq!(forces.len(), configuration.size());

        let natoms = configuration.size();
        let charges = configuration.particles().charge;

        for i in 0..natoms {
            let qi = charges[i];
            if qi == 0.0 {continue}
            for j in i+1..natoms {
                let qj = charges[j];
                if qj == 0.0 {continue}

                let distance = configuration.bond_distance(i, j);
                let info = self.restriction.information(distance);

                let rij = configuration.nearest_image(i, j);
                let force = self.real_space_force_pair(info, qi, qj, &rij);
                forces[i] += force;
                forces[j] -= force;
            }
        }
    }

    /// Real space contribution to the virial
    fn real_space_virial(&self, configuration: &Configuration) -> Matrix3 {
        let natoms = configuration.size();
        let charges = configuration.particles().charge;

        let mut virial = Matrix3::zero();
        for i in 0..natoms {
            let qi = charges[i];
            if qi == 0.0 {continue}
            for j in i+1..natoms {
                let qj = charges[j];
                if qj == 0.0 {continue}

                let distance = configuration.bond_distance(i, j);
                let info = self.restriction.information(distance);

                let rij = configuration.nearest_image(i, j);
                let force = self.real_space_force_pair(info, qi, qj, &rij);
                virial -= force.tensorial(&rij);
            }
        }
        return virial;
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
    fn density_fft(&mut self, configuration: &Configuration) {
        let natoms = configuration.size();
        let charges = configuration.particles().charge;
        let positions = configuration.particles().position;
        self.fourier_phases.resize_if_different((self.kmax, natoms, 3));

        // Do the k=0, 1 cases first
        for i in 0..natoms {
            let ri = configuration.cell.fractional(&positions[i]);
            for j in 0..3 {
                self.fourier_phases[(0, i, j)] = Complex::polar(1.0, 0.0);
                self.fourier_phases[(1, i, j)] = Complex::polar(1.0, -2.0 * PI * ri[j]);
            }
        }

        // Use recursive definition for computing the factor for all the other values of k.
        for k in 2..self.kmax {
            for i in 0..natoms {
                for j in 0..3 {
                    self.fourier_phases[(k, i, j)] = self.fourier_phases[(k - 1, i, j)]
                                                   * self.fourier_phases[(1, i, j)];
                }
            }
        }

        for ikx in 0..self.kmax {
            for iky in 0..self.kmax {
                for ikz in 0..self.kmax {
                    let mut rho = Complex::zero();
                    for j in 0..natoms {
                        let phi = self.fourier_phases[(ikx, j, 0)] * self.fourier_phases[(iky, j, 1)] * self.fourier_phases[(ikz, j, 2)];
                        rho += charges[j] * phi;
                    }
                    self.rho[(ikx, iky, ikz)] = rho;
                }
            }
        }
    }

    /// k-space contribution to the energy
    fn kspace_energy(&mut self, configuration: &Configuration) -> f64 {
        self.density_fft(configuration);
        let mut energy = 0.0;

        for ikx in 0..self.kmax {
            for iky in 0..self.kmax {
                for ikz in 0..self.kmax {
                    // The k = 0 case and the cutoff in k-space are already
                    // handled in `expfactors`
                    if self.expfactors[(ikx, iky, ikz)].abs() < f64::EPSILON {continue}
                    let density = self.rho[(ikx, iky, ikz)].norm();
                    energy += self.expfactors[(ikx, iky, ikz)] * density * density;
                }
            }
        }
        energy *= 2.0 * PI / (configuration.cell.volume() * ELCC);
        return energy;
    }

    /// k-space contribution to the forces
    fn kspace_forces(&mut self, configuration: &Configuration, forces: &mut [Vector3D]) {
        assert_eq!(forces.len(), configuration.size());
        self.density_fft(configuration);

        let charges = configuration.particles().charge;
        let factor = 4.0 * PI / (configuration.cell.volume() * ELCC);
        let (rec_kx, rec_ky, rec_kz) = configuration.cell.reciprocal_vectors();
        let thread_forces_store = ThreadLocalStore::new(|| vec![Vector3D::zero(); configuration.size()]);

        Zip::indexed(&*self.expfactors).par_apply(|(ikx, iky, ikz), &expfactor| {
            // The k = 0 and the cutoff in k-space are already handled
            // in `expfactors`.
            if expfactor < f64::EPSILON { return; }

            let f = expfactor * factor;
            let k = (ikx as f64) * rec_kx + (iky as f64) * rec_ky + (ikz as f64) * rec_kz;

            for i in 0..configuration.size() {
                let qi = charges[i];

                let fourier_i = self.fourier_phases[(ikx, i, 0)] *
                                self.fourier_phases[(iky, i, 1)] *
                                self.fourier_phases[(ikz, i, 2)];
                let fourier_i = fourier_i.imag();

                let mut thread_forces = thread_forces_store.borrow_mut();
                let mut force_i = Vector3D::zero();

                for j in (i + 1)..configuration.size() {
                    let qj = charges[j];
                    let force = f * self.kspace_force_factor(j, (ikx, iky, ikz), qi * qj, fourier_i) * k;
                    force_i -= force;
                    thread_forces[j] += force;
                }

                thread_forces[i] += force_i;
            }
        });

        thread_forces_store.sum_local_values(forces);
    }

    /// Get the force factor for particles `i` and `j` with charges `qi` and
    /// `qj`, at k point  `(ikx, iky, ikz)`
    #[inline]
    fn kspace_force_factor(&self, j: usize, (ikx, iky, ikz): (usize, usize, usize), qiqj: f64, fourier_i: f64) -> f64 {
        // Here the compiler is smart enough to optimize away
        // the useless computation of the last real part.
        let fourier_j = self.fourier_phases[(ikx, j, 0)] *
                        self.fourier_phases[(iky, j, 1)] *
                        self.fourier_phases[(ikz, j, 2)];
        let fourier_j = fourier_j.imag();

        return qiqj * (fourier_i - fourier_j);
    }

    /// k-space contribution to the virial
    fn kspace_virial(&mut self, configuration: &Configuration) -> Matrix3 {
        self.density_fft(configuration);

        let charges = configuration.particles().charge;
        let factor = 4.0 * PI / (configuration.cell.volume() * ELCC);
        let (rec_kx, rec_ky, rec_kz) = configuration.cell.reciprocal_vectors();

        Zip::indexed(&*self.expfactors).par_map(|((ikx, iky, ikz), &expfactor)| {
            if expfactor < f64::EPSILON { return Matrix3::zero(); }

            let f = expfactor * factor;
            let k = (ikx as f64) * rec_kx + (iky as f64) * rec_ky + (ikz as f64) * rec_kz;

            let mut local_virial = Matrix3::zero();
            for i in 0..configuration.size() {
                let qi = charges[i];

                let fourier_i = self.fourier_phases[(ikx, i, 0)] *
                                self.fourier_phases[(iky, i, 1)] *
                                self.fourier_phases[(ikz, i, 2)];
                let fourier_i = fourier_i.imag();

                for j in (i + 1)..configuration.size() {
                    let qj = charges[j];
                    let force = f * self.kspace_force_factor(j, (ikx, iky, ikz), qi * qj, fourier_i) * k;
                    let rij = configuration.nearest_image(i, j);
                    local_virial += force.tensorial(&rij);
                }
            }
            return local_virial;
        }).sum()
    }

    fn compute_delta_rho_move_particles(&mut self, configuration: &Configuration, idxes: &[usize], newpos: &[Vector3D]) {
        let natoms = idxes.len();
        let positions = configuration.particles().position;
        let charges = configuration.particles().charge;

        let mut new_fourier_phases = Array3::zeros((self.kmax, natoms, 3));
        let mut old_fourier_phases = Array3::zeros((self.kmax, natoms, 3));

        // Do the k=0, 1 cases first
        for (idx, &i) in idxes.iter().enumerate() {
            let old_ri = configuration.cell.fractional(&positions[i]);
            let new_ri = configuration.cell.fractional(&newpos[idx]);
            for j in 0..3 {
                old_fourier_phases[(0, idx, j)] = Complex::polar(1.0, 0.0);
                old_fourier_phases[(1, idx, j)] = Complex::polar(1.0, -2.0 * PI * old_ri[j]);

                new_fourier_phases[(0, idx, j)] = Complex::polar(1.0, 0.0);
                new_fourier_phases[(1, idx, j)] = Complex::polar(1.0, -2.0 * PI * new_ri[j]);
            }
        }

        // Use recursive definition for computing the factor for all the other values of k.
        for k in 2..self.kmax {
            for idx in 0..natoms {
                for j in 0..3 {
                    old_fourier_phases[(k, idx, j)] = old_fourier_phases[(k - 1, idx, j)]
                                                    * old_fourier_phases[(1, idx, j)];

                    new_fourier_phases[(k, idx, j)] = new_fourier_phases[(k - 1, idx, j)]
                                                    * new_fourier_phases[(1, idx, j)];
                }
            }
        }

        for ikx in 0..self.kmax {
            for iky in 0..self.kmax {
                for ikz in 0..self.kmax {
                    self.delta_rho[(ikx, iky, ikz)] = Complex::polar(0.0, 0.0);
                    for (idx, &i) in idxes.iter().enumerate() {
                        let old_phi = old_fourier_phases[(ikx, idx, 0)] *
                                      old_fourier_phases[(iky, idx, 1)] *
                                      old_fourier_phases[(ikz, idx, 2)];
                        let new_phi = new_fourier_phases[(ikx, idx, 0)] *
                                      new_fourier_phases[(iky, idx, 1)] *
                                      new_fourier_phases[(ikz, idx, 2)];

                        self.delta_rho[(ikx, iky, ikz)] -= charges[i] * old_phi;
                        self.delta_rho[(ikx, iky, ikz)] += charges[i] * new_phi;
                    }
                }
            }
        }
    }

    fn kspace_move_particles_cost(&mut self, configuration: &Configuration, idxes: &[usize], newpos: &[Vector3D]) -> f64 {
        let e_old = self.kspace_energy(configuration);

        let mut e_new = 0.0;
        self.compute_delta_rho_move_particles(configuration, idxes, newpos);
        for ikx in 0..self.kmax {
            for iky in 0..self.kmax {
                for ikz in 0..self.kmax {
                    // The k = 0 case and the cutoff in k-space are already
                    // handled in `expfactors`.
                    if self.expfactors[(ikx, iky, ikz)].abs() < f64::EPSILON {continue}
                    let rho = self.rho[(ikx, iky, ikz)] + self.delta_rho[(ikx, iky, ikz)];
                    let density = rho.norm();
                    e_new += self.expfactors[(ikx, iky, ikz)] * density * density;
                }
            }
        }
        e_new *= 2.0 * PI / (configuration.cell.volume() * ELCC);

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
                virial -= force.tensorial(&rij);
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
    /// let ewald = SharedEwald::new(Ewald::new(12.5, 10));
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
        /* No self force */
        ewald.kspace_forces(configuration, forces);
        ewald.molcorrect_forces(configuration, forces);
    }

    fn virial(&self, configuration: &Configuration) -> Matrix3 {
        let mut ewald = self.write();
        ewald.precompute(&configuration.cell);
        let real = ewald.real_space_virial(configuration);
        /* No self virial */
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
        let mut ewald = self.write();
        for ikx in 0..ewald.kmax {
            for iky in 0..ewald.kmax {
                for ikz in 0..ewald.kmax {
                    ewald.rho[(ikx, iky, ikz)] += ewald.delta_rho[(ikx, iky, ikz)];
                }
            }
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
            system.cell = UnitCell::new();
            let ewald = SharedEwald::new(Ewald::new(8.0, 10));
            let _ = ewald.energy(&system);
        }

        #[test]
        #[should_panic]
        fn triclinic_cell() {
            let mut system = nacl_pair();
            system.cell = UnitCell::triclinic(10.0, 10.0, 10.0, 90.0, 90.0, 90.0);
            let ewald = SharedEwald::new(Ewald::new(8.0, 10));
            let _ = ewald.energy(&system);
        }

        #[test]
        #[should_panic]
        fn negative_alpha() {
            let mut ewald = Ewald::new(8.0, 10);
            ewald.set_alpha(-45.2);
        }
    }

    mod pairs {
        use super::*;
        use energy::GlobalPotential;

        #[test]
        fn energy() {
            let system = nacl_pair();
            let ewald = SharedEwald::new(Ewald::new(8.0, 10));

            let energy = ewald.energy(&system);
            // This was computed by hand
            let energy_brute_force = -0.09262397663346732;
            assert_ulps_eq!(energy, energy_brute_force, epsilon=1e-4);
        }

        #[test]
        fn forces() {
            let mut system = nacl_pair();
            let ewald = SharedEwald::new(Ewald::new(8.0, 10));

            let mut forces = vec![Vector3D::zero(); 2];
            ewald.forces(&system, &mut forces);
            let norm = (forces[0] + forces[1]).norm();
            // Total force should be null
            assert_ulps_eq!(norm, 0.0);

            // Force is attractive
            for i in 0..3 {
                assert!(forces[0][i] > 0.0);
                assert!(forces[1][i] < 0.0);
            }

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
            let mut ewald = SharedEwald::new(Ewald::new(8.0, 10));
            ewald.set_restriction(PairRestriction::InterMolecular);

            let energy = ewald.energy(&system);
            let expected = 0.0002257554843856993;
            assert_ulps_eq!(energy, expected);

            let molcorrect = ewald.read().molcorrect_energy(&system);
            let expected = 0.02452968743897957;
            assert_ulps_eq!(molcorrect, expected);
        }

        #[test]
        fn forces() {
            let mut system = water();
            let mut ewald = SharedEwald::new(Ewald::new(8.0, 10));
            ewald.set_restriction(PairRestriction::InterMolecular);

            let mut forces = vec![Vector3D::zero(); 3];
            ewald.forces(&system, &mut forces);
            let norm = (forces[0] + forces[1]).norm();
            // Total force should be null
            assert!(norm.abs() < 1e-3);

            // Finite difference computation of all the force components
            let energy = ewald.energy(&system);
            let real_energy = ewald.read().real_space_energy(&system);
            let kspace_energy = ewald.write().kspace_energy(&system);
            let molcorrect_energy = ewald.read().molcorrect_energy(&system);

            let eps = 1e-9;
            system.particles_mut().position[0][0] += eps;

            let energy_1 = ewald.energy(&system);
            let real_energy_1 = ewald.read().real_space_energy(&system);
            let kspace_energy_1 = ewald.write().kspace_energy(&system);
            let molcorrect_energy_1 = ewald.read().molcorrect_energy(&system);

            let mut forces = vec![Vector3D::zero(); 3];
            ewald.forces(&system, &mut forces);
            let force = forces[0][0];

            let mut forces_buffer = vec![Vector3D::zero(); system.size()];
            ewald.read().real_space_forces(&system, &mut forces_buffer);
            let real_force = forces_buffer[0][0];

            let mut forces_buffer = vec![Vector3D::zero(); system.size()];
            ewald.write().kspace_forces(&system, &mut forces_buffer);
            let kspace_force = forces_buffer[0][0];

            let mut forces_buffer = vec![Vector3D::zero(); system.size()];
            ewald.read().molcorrect_forces(&system, &mut forces_buffer);
            let molcorrect_force = forces_buffer[0][0];

            let force_fda = (energy - energy_1) / eps;
            assert!(abs((force_fda - force) / force) < 1e-4);

            // No real space energetic contribution here, we only have one
            // molecule.
            assert_ulps_eq!(real_energy, 0.0);
            assert_ulps_eq!(real_energy_1, 0.0);
            assert_ulps_eq!(real_force, 0.0);

            let kspace_force_fda = (kspace_energy - kspace_energy_1) / eps;
            assert!(abs((kspace_force_fda - kspace_force) / kspace_force) < 1e-4);

            let molcorrect_force_fda = (molcorrect_energy - molcorrect_energy_1) / eps;
            assert!(abs((molcorrect_force_fda - molcorrect_force) / molcorrect_force) < 1e-4);
        }
    }

    mod virial {
        use super::*;
        use types::{Vector3D, Zero};
        use energy::{GlobalPotential, PairRestriction, CoulombicPotential};

        #[test]
        fn real_space() {
            let system = nacl_pair();
            let ewald = Ewald::new(8.0, 10);

            // real space
            let virial = ewald.real_space_virial(&system);
            let mut forces = vec![Vector3D::zero(); 2];
            ewald.real_space_forces(&system, &mut forces);
            let expected = forces[0].tensorial(&Vector3D::new(1.5, 0.0, 0.0));
            assert_ulps_eq!(virial, expected);
        }

        #[test]
        fn kspace() {
            let system = nacl_pair();
            let mut ewald = Ewald::new(8.0, 10);

            let virial = ewald.kspace_virial(&system);
            let mut forces = vec![Vector3D::zero(); 2];
            ewald.kspace_forces(&system, &mut forces);
            let expected = forces[0].tensorial(&Vector3D::new(1.5, 0.0, 0.0));
            assert_ulps_eq!(virial, expected);
        }

        #[test]
        fn molcorrect() {
            let mut system = nacl_pair();
            let _ = system.add_bond(0, 1);

            let mut ewald = SharedEwald::new(Ewald::new(8.0, 10));
            ewald.set_restriction(PairRestriction::InterMolecular);

            let virial = ewald.read().molcorrect_virial(&system);
            let mut forces = vec![Vector3D::zero(); 2];
            ewald.read().molcorrect_forces(&system, &mut forces);
            let expected = forces[0].tensorial(&Vector3D::new(1.5, 0.0, 0.0));
            assert_ulps_eq!(virial, expected);
        }

        #[test]
        fn total() {
            let system = nacl_pair();
            let ewald = SharedEwald::new(Ewald::new(8.0, 10));
            let mut forces = vec![Vector3D::zero(); 2];

            let virial = ewald.virial(&system);
            ewald.forces(&system, &mut forces);
            let expected = forces[0].tensorial(&Vector3D::new(1.5, 0.0, 0.0));
            assert_ulps_eq!(virial, expected, max_ulps=25);
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
            let mut ewald = SharedEwald::new(Ewald::new(8.0, 10));
            ewald.set_restriction(PairRestriction::InterMolecular);

            let ewald_check = ewald.clone();

            let old_e = ewald_check.energy(&system);
            let idxes = &[0, 1];
            let newpos = &[Vector3D::new(0.0, 0.0, 0.5), Vector3D::new(-0.7, 0.2, 1.5)];

            let cost = ewald.move_particles_cost(&system, idxes, newpos);

            system.particles_mut().position[0] = newpos[0];
            system.particles_mut().position[1] = newpos[1];
            let new_e = ewald_check.energy(&system);
            assert_ulps_eq!(cost, new_e - old_e);
        }

        #[test]
        fn move_atoms_real_space() {
            let mut system = testing_system();
            let mut ewald = SharedEwald::new(Ewald::new(8.0, 10));
            ewald.set_restriction(PairRestriction::InterMolecular);

            let ewald_check = ewald.clone();

            let old_e = ewald_check.read().real_space_energy(&system);
            let idxes = &[0, 1];
            let newpos = &[Vector3D::new(0.0, 0.0, 0.5), Vector3D::new(-0.7, 0.2, 1.5)];

            let cost = ewald.read().real_space_move_particles_cost(&system, idxes, newpos);

            system.particles_mut().position[0] = newpos[0];
            system.particles_mut().position[1] = newpos[1];
            let new_e = ewald_check.read().real_space_energy(&system);
            assert_ulps_eq!(cost, new_e - old_e);
        }

        #[test]
        fn move_atoms_kspace() {
            let mut system = testing_system();
            let mut ewald = SharedEwald::new(Ewald::new(8.0, 10));
            ewald.set_restriction(PairRestriction::InterMolecular);

            let ewald_check = ewald.clone();

            let old_e = ewald_check.write().kspace_energy(&system);
            let idxes = &[0, 1];
            let newpos = &[Vector3D::new(0.0, 0.0, 0.5), Vector3D::new(-0.7, 0.2, 1.5)];

            let cost = ewald.write().kspace_move_particles_cost(&system, idxes, newpos);

            system.particles_mut().position[0] = newpos[0];
            system.particles_mut().position[1] = newpos[1];
            let new_e = ewald_check.write().kspace_energy(&system);
            assert_ulps_eq!(cost, new_e - old_e);
        }

        #[test]
        fn move_atoms_molcorrect() {
            let mut system = testing_system();
            let mut ewald = SharedEwald::new(Ewald::new(8.0, 10));
            ewald.set_restriction(PairRestriction::InterMolecular);

            let ewald_check = ewald.clone();

            let old_e = ewald_check.read().molcorrect_energy(&system);
            let idxes = &[0, 1];
            let newpos = &[Vector3D::new(0.0, 0.0, 0.5), Vector3D::new(-0.7, 0.2, 1.5)];

            let cost = ewald.write().molcorrect_move_particles_cost(&system, idxes, newpos);

            system.particles_mut().position[0] = newpos[0];
            system.particles_mut().position[1] = newpos[1];
            let new_e = ewald_check.read().molcorrect_energy(&system);
            assert_ulps_eq!(cost, new_e - old_e);
        }
    }
}
