// Cymbalum, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

//! Ewald summation of coulombic potential
use special::{erfc, erf};
use std::f64::consts::{PI, FRAC_2_SQRT_PI};
use std::f64;

use system::{System, UnitCell, CellType};
use types::{Matrix3, Vector3D, Array3, Complex, Zero};
use constants::ELCC;
use potentials::{PairRestriction, RestrictionInfo};

use super::{GlobalPotential, CoulombicPotential, GlobalCache};

/// Ewald summation of the coulombic interactions. The Ewald summation is based
/// on a separation of the coulombic potential U in two parts, using the trivial
/// identity:
///
///     U(x) = U(x) * (f(x) + 1) - U(x) * f(x)
///
/// where `f(x)` is the `erf` function. This leads to a separation of the
/// conditionally convergent coulombic sum into two absolutly convergent sums:
/// one in real space, and the other in Fourier or k-space.
///
/// For more informations about this algorithm see [FS2002].
///
/// [FS2002] Frenkel, D. & Smit, B. Understanding molecular simulation. (Academic press, 2002).
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
    fourier_phases: Array3<f64>,
    /// Fourier transform of the electrostatic density
    rho: Array3<Complex>,
    /// Fourier transform of the electrostatic density modifications, cached
    /// allocation and for updating `self.rho`
    delta_rho: Array3<Complex>,
    /// Guard for cache invalidation of expfactors
    previous_cell: Option<UnitCell>,
}

impl Ewald {
    /// Create an Ewald summation using the `rc` cutoff radius in real space,
    /// and `kmax` points in k-space (Fourier space).
    pub fn new(rc: f64, kmax: usize) -> Ewald {
        let kmax = kmax + 1;
        let expfactors = Array3::zeros((kmax, kmax, kmax));
        let rho = Array3::zeros((kmax, kmax, kmax));
        Ewald {
            alpha: 3.0 * PI / (rc * 4.0),
            rc: rc,
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

    fn precompute(&mut self, cell: &UnitCell) {
        if let Some(ref prev_cell) = self.previous_cell {
            if cell == prev_cell {
                // Do not recompute
                return;
            }
        }
        match cell.celltype() {
            CellType::Infinite => {
                fatal_error!("Can not use Ewald sum with Infinite cell.");
            },
            CellType::Triclinic => {
                fatal_error!("Can not (yet) use Ewald sum with Triclinic cell.");
            },
            CellType::Orthorombic => {
                // All good!
            },
        }
        self.previous_cell = Some(cell.clone());

        // Because we do a spherical truncation in k space, we have to transform
        // kmax into a spherical cutoff 'radius'
        let lenghts = cell.lengths();
        let max_lenght = f64::max(f64::max(lenghts.0, lenghts.1), lenghts.2);
        let min_lenght = f64::min(f64::min(lenghts.0, lenghts.1), lenghts.2);
        let k_rc = self.kmax as f64 * (2.0 * PI / max_lenght);
        self.kmax2 = k_rc * k_rc;

        if self.rc > min_lenght / 2.0 {
            warn!("The Ewald cutoff is too high for this unit cell, energy might be wrong.");
        }

        // Now, we precompute the exp(-k^2/4a^2)/k^2 terms. We use the symmetry to
        // only store (ikx >= 0 && iky >= 0  && ikz >= 0 ) terms
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
                    self.expfactors[(ikx, iky, ikz)] = f64::exp(-k2 / (4.0 * self.alpha * self.alpha)) / k2;
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
    /// and `qj` ; and with restriction informations for this pair in `info`.
    #[inline(always)]
    fn real_space_energy_pair(&self, info: RestrictionInfo, qi: f64, qj: f64, r: f64) -> f64 {
        if r > self.rc || info.excluded {
            return 0.0
        }
        assert!(info.scaling  == 1.0, "Scaling restriction scheme using Ewald are not implemented");
        return qi * qj * erfc(self.alpha * r) / r / ELCC;
    }

    /// Get the real-space force for one pair at distance `rij` with charges
    /// `qi` and `qj` ; and with restriction informations for this pair in
    /// `info`.
    #[inline(always)]
    fn real_space_force_pair(&self, info: RestrictionInfo, qi: f64, qj: f64, rij: &Vector3D) -> Vector3D {
        let r = rij.norm();
        if r > self.rc || info.excluded {
            return Vector3D::new(0.0, 0.0, 0.0)
        }
        let mut factor = erfc(self.alpha * r) / r;
        factor += self.alpha * FRAC_2_SQRT_PI * f64::exp(-self.alpha * self.alpha * r * r);
        factor *= qi * qj / (r * r) / ELCC;
        return factor * rij;
    }

    /// Real space contribution to the energy
    fn real_space_energy(&self, system: &System) -> f64 {
        let natoms = system.size();
        let mut energy = 0.0;
        for i in 0..natoms {
            let qi = system[i].charge;
            if qi == 0.0 {continue}
            for j in i+1..natoms {
                let qj = system[j].charge;
                if qj == 0.0 {continue}

                let r = system.distance(i, j);
                let info = self.restriction.informations(system, i, j);
                energy += self.real_space_energy_pair(info, qi, qj, r);
            }
        }
        return energy;
    }

    /// Real space contribution to the forces
    fn real_space_forces(&self, system: &System, res: &mut Vec<Vector3D>) {
        let natoms = system.size();
        assert!(res.len() == system.size());

        for i in 0..natoms {
            let qi = system[i].charge;
            if qi == 0.0 {continue}
            for j in i+1..natoms {
                let qj = system[j].charge;
                if qj == 0.0 {continue}

                let rij = system.wraped_vector(i, j);
                let info = self.restriction.informations(system, i, j);

                let force = self.real_space_force_pair(info, qi, qj, &rij);
                res[i] = res[i] + force;
                res[j] = res[j] - force;
            }
        }
    }

    /// Real space contribution to the virial
    fn real_space_virial(&self, system: &System) -> Matrix3 {
        let natoms = system.size();
        let mut res = Matrix3::zero();
        for i in 0..natoms {
            let qi = system[i].charge;
            if qi == 0.0 {continue}
            for j in i+1..natoms {
                let qj = system[j].charge;
                if qj == 0.0 {continue}

                let rij = system.wraped_vector(i, j);
                let info = self.restriction.informations(system, i, j);

                let force = self.real_space_force_pair(info, qi, qj, &rij);
                res = res - force.tensorial(&rij);
            }
        }
        return res;
    }

    fn real_space_move_particles_cost(&self, system: &System, idxes: &[usize], newpos: &[Vector3D]) -> f64 {
        let mut e_old = 0.0;
        let mut e_new = 0.0;

        // Iterate over all interactions between a moved particle and a
        // particle not moved
        for (idx, &i) in idxes.iter().enumerate() {
            let qi = system[i].charge;
            if qi == 0.0 {continue}
            for j in (0..system.size()).filter(|x| !idxes.contains(x)) {
                let qj = system[j].charge;
                if qi == 0.0 {continue}

                let r_old = system.distance(i, j);
                let r_new = system.cell().distance(&newpos[idx], &system[j].position);
                let info = self.restriction.informations(system, i, j);

                e_old += self.real_space_energy_pair(info, qi, qj, r_old);
                e_new += self.real_space_energy_pair(info, qi, qj, r_new);
            }
        }

        // Iterate over all interactions between two moved particles
        for (idx, &i) in idxes.iter().enumerate() {
            let qi = system[i].charge;
            if qi == 0.0 {continue}
            for (jdx, &j) in idxes.iter().enumerate().skip(i + 1) {
                let qj = system[j].charge;
                if qj == 0.0 {continue}

                let r_old = system.distance(i, j);
                let r_new = system.cell().distance(&newpos[idx], &newpos[jdx]);

                let info = self.restriction.informations(system, i, j);

                e_old += self.real_space_energy_pair(info, qi, qj, r_old);
                e_new += self.real_space_energy_pair(info, qi, qj, r_new);
            }
        }

        return e_new - e_old;
    }
}

/// Self-interaction corection
impl Ewald {
    /// Self-interaction contribution to the energy
    fn self_energy(&self, system: &System) -> f64 {
        let mut q2 = 0.0;
        for i in 0..system.size() {
            q2 += system[i].charge * system[i].charge;
        }
        return -self.alpha / f64::sqrt(PI) * q2 / ELCC;
    }
}

/// k-space part of the summation
impl Ewald {
    /// Compute the Fourier transform of the electrostatic density
    fn density_fft(&mut self, system: &System) {
        let natoms = system.size();
        self.fourier_phases.resize_if_different((self.kmax, natoms, 3));

        // Do the k=0, 1 cases first
        for i in 0..natoms {
            let ri = system.cell().fractional(&system[i].position);
            for j in 0..3 {
                self.fourier_phases[(0, i, j)] = 0.0;
                self.fourier_phases[(1, i, j)] = -2.0 * PI * ri[j];
            }
        }

        // Use recursive definition for computing the factor for all the other values of k.
        for k in 2..self.kmax {
            for i in 0..natoms {
                for j in 0..3 {
                    self.fourier_phases[(k, i, j)] = self.fourier_phases[(k - 1, i, j)]
                                                   + self.fourier_phases[(1, i, j)];
                }
            }
        }

        for ikx in 0..self.kmax {
            for iky in 0..self.kmax {
                for ikz in 0..self.kmax {
                    self.rho[(ikx, iky, ikz)] = Complex::polar(0.0, 0.0);
                    for j in 0..natoms {
                        let phi = self.fourier_phases[(ikx, j, 0)] + self.fourier_phases[(iky, j, 1)] + self.fourier_phases[(ikz, j, 2)];
                        self.rho[(ikx, iky, ikz)] = self.rho[(ikx, iky, ikz)] + Complex::polar(system[j].charge, phi);
                    }
                }
            }
        }
    }

    /// k-space contribution to the energy
    fn kspace_energy(&mut self, system: &System) -> f64 {
        self.density_fft(system);
        let mut energy = 0.0;

        for ikx in 0..self.kmax {
            for iky in 0..self.kmax {
                for ikz in 0..self.kmax {
                    // The k = 0 case and the cutoff in k-space are already
                    // handled in expfactors
                    if self.expfactors[(ikx, iky, ikz)].abs() < f64::MIN {continue}
                    let density = self.rho[(ikx, iky, ikz)].norm();
                    energy += self.expfactors[(ikx, iky, ikz)] * density * density;
                }
            }
        }
        energy *= 2.0 * PI / (system.cell().volume() * ELCC);
        return energy;
    }

    /// k-space contribution to the forces
    fn kspace_forces(&mut self, system: &System, res: &mut Vec<Vector3D>) {
        assert!(res.len() == system.size());
        self.density_fft(system);

        let factor = 4.0 * PI / (system.cell().volume() * ELCC);
        let (rec_kx, rec_ky, rec_kz) = system.cell().reciprocal_vectors();

        for ikx in 0..self.kmax {
            for iky in 0..self.kmax {
                for ikz in 0..self.kmax {
                    // The k = 0 and the cutoff in k-space are already handled in
                    // expfactors.
                    if self.expfactors[(ikx, iky, ikz)].abs() < f64::MIN {continue}
                    let k = (ikx as f64) * rec_kx + (iky as f64) * rec_ky + (ikz as f64) * rec_kz;
                    for i in 0..system.size() {
                        let qi = system[i].charge;
                        for j in (i + 1)..system.size() {
                            let qj = system[j].charge;
                            let force = factor * self.kspace_force_factor(i, j, ikx, iky, ikz, qi, qj) * k;

                            res[i] = res[i] - force;
                            res[j] = res[j] + force;
                        }
                    }
                }
            }
        }
    }

    /// Get the force factor for particles `i` and `j` with charges `qi` and
    /// `qj`, at k point  `(ikx, iky, ikz)`
    #[inline(always)]
    fn kspace_force_factor(&self, i: usize, j: usize, ikx: usize, iky: usize, ikz: usize, qi: f64, qj: f64) -> f64 {
        let fourier_i = self.fourier_phases[(ikx, i, 0)]
                      + self.fourier_phases[(iky, i, 1)]
                      + self.fourier_phases[(ikz, i, 2)];
        let fourier_j = self.fourier_phases[(ikx, j, 0)]
                      + self.fourier_phases[(iky, j, 1)]
                      + self.fourier_phases[(ikz, j, 2)];
        let sin_kr = fast_sin(fourier_i - fourier_j);

        return qi * qj * self.expfactors[(ikx, iky, ikz)] * sin_kr;
    }

    /// k-space contribution to the virial
    fn kspace_virial(&mut self, system: &System) -> Matrix3 {
        self.density_fft(system);
        let mut res = Matrix3::zero();

        let factor = 4.0 * PI / (system.cell().volume() * ELCC);
        let (rec_kx, rec_ky, rec_kz) = system.cell().reciprocal_vectors();

        for ikx in 0..self.kmax {
            for iky in 0..self.kmax {
                for ikz in 0..self.kmax {
                    // The k = 0 and the cutoff in k-space are already handled in
                    // expfactors.
                    if self.expfactors[(ikx, iky, ikz)].abs() < f64::MIN {continue}
                    let k = (ikx as f64) * rec_kx + (iky as f64) * rec_ky + (ikz as f64) * rec_kz;
                    for i in 0..system.size() {
                        let qi = system[i].charge;
                        for j in (i + 1)..system.size() {
                            let qj = system[j].charge;
                            let force = factor * self.kspace_force_factor(i, j, ikx, iky, ikz, qi, qj) * k;
                            let rij = system.wraped_vector(i, j);

                            res = res + force.tensorial(&rij);
                        }
                    }
                }
            }
        }
        return res;
    }

    fn compute_delta_rho_move_particles(&mut self, system: &System, idxes: &[usize], newpos: &[Vector3D]) {
        let natoms = idxes.len();
        let mut new_fourier_phases = Array3::zeros((self.kmax, natoms, 3));
        let mut old_fourier_phases = Array3::zeros((self.kmax, natoms, 3));

        // Do the k=0, 1 cases first
        for (idx, &i) in idxes.iter().enumerate() {
            let old_ri = system.cell().fractional(&system[i].position);
            let new_ri = system.cell().fractional(&newpos[idx]);
            for j in 0..3 {
                old_fourier_phases[(0, i, j)] = 0.0;
                old_fourier_phases[(1, i, j)] = -2.0 * PI * old_ri[j];

                new_fourier_phases[(0, i, j)] = 0.0;
                new_fourier_phases[(1, i, j)] = -2.0 * PI * new_ri[j];
            }
        }

        // Use recursive definition for computing the factor for all the other values of k.
        for k in 2..self.kmax {
            for i in 0..natoms {
                for j in 0..3 {
                    old_fourier_phases[(k, i, j)] = old_fourier_phases[(k - 1, i, j)]
                                                   + old_fourier_phases[(1, i, j)];

                    new_fourier_phases[(k, i, j)] = new_fourier_phases[(k - 1, i, j)]
                                                   + new_fourier_phases[(1, i, j)];
                }
            }
        }

        for ikx in 0..self.kmax {
            for iky in 0..self.kmax {
                for ikz in 0..self.kmax {
                    self.delta_rho[(ikx, iky, ikz)] = Complex::polar(0.0, 0.0);
                    for j in 0..natoms {
                        let old_phi = old_fourier_phases[(ikx, j, 0)] + old_fourier_phases[(iky, j, 1)] + old_fourier_phases[(ikz, j, 2)];
                        let new_phi = new_fourier_phases[(ikx, j, 0)] + new_fourier_phases[(iky, j, 1)] + new_fourier_phases[(ikz, j, 2)];

                        self.delta_rho[(ikx, iky, ikz)] = self.delta_rho[(ikx, iky, ikz)] - Complex::polar(system[j].charge, old_phi);
                        self.delta_rho[(ikx, iky, ikz)] = self.delta_rho[(ikx, iky, ikz)] + Complex::polar(system[j].charge, new_phi);
                    }
                }
            }
        }
    }

    fn kspace_move_particles_cost(&mut self, system: &System, idxes: &[usize], newpos: &[Vector3D]) -> f64 {
        let e_old = self.kspace_energy(system);

        let mut e_new = 0.0;
        self.compute_delta_rho_move_particles(system, idxes, newpos);
        for ikx in 0..self.kmax {
            for iky in 0..self.kmax {
                for ikz in 0..self.kmax {
                    // The k = 0 case and the cutoff in k-space are already
                    // handled in expfactors
                    if self.expfactors[(ikx, iky, ikz)].abs() < f64::MIN {continue}
                    let rho = self.rho[(ikx, iky, ikz)] + self.delta_rho[(ikx, iky, ikz)];
                    let density = rho.norm();
                    e_new += self.expfactors[(ikx, iky, ikz)] * density * density;
                }
            }
        }
        e_new *= 2.0 * PI / (system.cell().volume() * ELCC);

        return e_new - e_old;
    }
}

/// This an implementation of the sin function which is faster than the libm
/// sin function on my i7 processor. I'll have to benchmarck this on other
/// architectures and OS.
///
/// Using this function in `kspace_force_factor` gives me a 40% speedup of the
/// overall simulation time.
///
/// This code comes from http://forum.devmaster.net/t/fast-and-accurate-sine-cosine/9648/85
fn fast_sin(mut x: f64) -> f64 {
    const FRAC_1_TWO_PI: f64 = 1.0 / (2.0 * PI);
    const A: f64 = 7.58946638440411;
    const B: f64 = 1.6338434577536627;

    x = x * FRAC_1_TWO_PI;
    x = x - f64::floor(x + 0.5);
    x = A * x * (0.5 - f64::abs(x));
    x = x * (B + f64::abs(x));

    return x;
}

/// Molecular correction for Ewald summation
impl Ewald {
    /// Get the molecular correction energy for the pair with charges `qi` and
    /// `qj`, at distance `rij` and with restriction informations in `info`.
    #[inline(always)]
    fn molcorrect_energy_pair(&self, info: RestrictionInfo, qi: f64, qj: f64, r: f64) -> f64 {
        assert!(info.excluded, "Can not compute molecular correction for non-excluded pair");
        assert!(info.scaling  == 1.0, "Scaling restriction scheme using Ewald are not implemented");
        assert!(r < self.rc, "Atoms in molecule are separated by more than the cutoff radius of Ewald sum.");

        return 0.5 * qi * qj / ELCC * erf(self.alpha * r)/r;
    }

    /// Get the molecular correction force for the pair with charges `qi` and
    /// `qj`, at distance `rij` and with restriction informations in `info`.
    #[inline(always)]
    fn molcorrect_force_pair(&self, info: RestrictionInfo, qi: f64, qj: f64, rij: &Vector3D) -> Vector3D {
        assert!(info.excluded, "Can not compute molecular correction for non-excluded pair");
        assert!(info.scaling  == 1.0, "Scaling restriction scheme using Ewald are not implemented");
        let r = rij.norm();
        assert!(r < self.rc, "Atoms in molecule are separated by more than the cutoff radius of Ewald sum.");

        let qiqj = qi * qj / (ELCC * r * r);
        let factor = qiqj * (2.0 * self.alpha / f64::sqrt(PI) * f64::exp(-self.alpha * self.alpha * r * r) - erf(self.alpha * r) / r);
        return factor * rij;
    }

    /// Molecular correction contribution to the energy
    fn molcorrect_energy(&self, system: &System) -> f64 {
        let natoms = system.size();
        let mut energy = 0.0;

        for i in 0..natoms {
            let qi = system[i].charge;
            if qi == 0.0 {continue}
            // I can not manage to get this work with a loop from (i+1) to N. The finite
            // difference test (testing that the force is the same that the finite difference
            // of the energy) always fail. So let's use it that way for now.
            for j in 0..natoms {
                if i == j {continue}
                // Only account for excluded pairs
                let info = self.restriction.informations(system, i, j);
                if !info.excluded {continue}

                let qj = system[j].charge;
                if qj == 0.0 {continue}

                let r = system.distance(i, j);
                energy += self.molcorrect_energy_pair(info, qi, qj, r);
            }
        }
        return energy;
    }

    /// Molecular correction contribution to the forces
    fn molcorrect_forces(&self, system: &System, res: &mut Vec<Vector3D>) {
        let natoms = system.size();
        assert!(res.len() == natoms);

        for i in 0..natoms {
            let qi = system[i].charge;
            if qi == 0.0 {continue}
            for j in 0..natoms {
                if i == j {continue}
                let info = self.restriction.informations(system, i, j);
                // Only account for excluded pairs
                if !info.excluded {continue}

                let qj = system[j].charge;
                if qi == 0.0 {continue}

                let rij = system.wraped_vector(i, j);
                let force = self.molcorrect_force_pair(info, qi, qj, &rij);
                res[i] = res[i] - force;
            }
        }
    }

    /// Molecular correction contribution to the virial
    fn molcorrect_virial(&self, system: &System) -> Matrix3 {
        let natoms = system.size();
        let mut res = Matrix3::zero();

        for i in 0..natoms {
            let qi = system[i].charge;
            if qi == 0.0 {continue}
            for j in 0..natoms {
                if i == j {continue}
                let info = self.restriction.informations(system, i, j);
                // Only account for excluded pairs
                if !info.excluded {continue}

                let qj = system[j].charge;
                if qi == 0.0 {continue}

                let rij = system.wraped_vector(i, j);
                let force = self.molcorrect_force_pair(info, qi, qj, &rij);
                res = res + force.tensorial(&rij);
            }
        }
        return res;
    }

    fn molcorrect_move_particles_cost(&mut self, system: &System, idxes: &[usize], newpos: &[Vector3D]) -> f64 {
        let mut e_old = 0.0;
        let mut e_new = 0.0;

        // Iterate over all interactions between a moved particle and a
        // particle not moved
        for (idx, &i) in idxes.iter().enumerate() {
            let qi = system[i].charge;
            if qi == 0.0 {continue}
            for j in (0..system.size()).filter(|x| !idxes.contains(x)) {
                let qj = system[j].charge;
                if qi == 0.0 {continue}

                let info = self.restriction.informations(system, i, j);
                if !info.excluded {continue}

                let r_old = system.distance(i, j);
                let r_new = system.cell().distance(&newpos[idx], &system[j].position);

                e_old += self.molcorrect_energy_pair(info, qi, qj, r_old);
                e_new += self.molcorrect_energy_pair(info, qi, qj, r_new);
            }
        }

        // Iterate over all interactions between two moved particles
        for (idx, &i) in idxes.iter().enumerate() {
            let qi = system[i].charge;
            if qi == 0.0 {continue}
            for (jdx, &j) in idxes.iter().enumerate().skip(i + 1) {
                let qj = system[j].charge;
                if qj == 0.0 {continue}

                let info = self.restriction.informations(system, i, j);
                if !info.excluded {continue}

                let r_old = system.distance(i, j);
                let r_new = system.cell().distance(&newpos[idx], &newpos[jdx]);

                e_old += self.molcorrect_energy_pair(info, qi, qj, r_old);
                e_new += self.molcorrect_energy_pair(info, qi, qj, r_new);
            }
        }

        return e_new - e_old;
    }
}

impl GlobalPotential for Ewald {
    fn energy(&mut self, system: &System) -> f64 {
        self.precompute(system.cell());
        let real = self.real_space_energy(system);
        let self_e = self.self_energy(system);
        let kspace = self.kspace_energy(system);
        let molecular = self.molcorrect_energy(system);
        return real + self_e + kspace + molecular;
    }

    fn forces(&mut self, system: &System) -> Vec<Vector3D> {
        self.precompute(system.cell());
        let mut res = vec![Vector3D::zero(); system.size()];
        self.real_space_forces(system, &mut res);
        /* No self force */
        self.kspace_forces(system, &mut res);
        self.molcorrect_forces(system, &mut res);
        return res;
    }

    fn virial(&mut self, system: &System) -> Matrix3 {
        self.precompute(system.cell());
        let real = self.real_space_virial(system);
        /* No self virial */
        let kspace = self.kspace_virial(system);
        let molecular = self.molcorrect_virial(system);
        return real + kspace + molecular;
    }
}

impl CoulombicPotential for Ewald {
    fn set_restriction(&mut self, restriction: PairRestriction) {
        self.restriction = restriction;
    }
}

impl GlobalCache for Ewald {
    fn move_particles_cost(&mut self, system: &System, idxes: &[usize], newpos: &[Vector3D]) -> f64 {
        self.precompute(system.cell());
        let real = self.real_space_move_particles_cost(system, idxes, newpos);
        /* No self cost */
        let kspace = self.kspace_move_particles_cost(system, idxes, newpos);
        let molecular = self.molcorrect_move_particles_cost(system, idxes, newpos);
        return real + kspace + molecular;
    }

    fn update(&mut self) {
        for ikx in 0..self.kmax {
            for iky in 0..self.kmax {
                for ikz in 0..self.kmax {
                    self.rho[(ikx, iky, ikz)] = self.rho[(ikx, iky, ikz)] + self.delta_rho[(ikx, iky, ikz)];
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    pub use super::*;
    use system::{System, UnitCell, Particle};
    use types::{Vector3D, Zero};
    use potentials::GlobalPotential;

    const E_BRUTE_FORCE: f64 = -0.09262397663346732;

    pub fn testing_system() -> System {
        let mut system = System::from_cell(UnitCell::cubic(20.0));

        system.add_particle(Particle::new("Cl"));
        system[0].charge = -1.0;
        system[0].position = Vector3D::zero();

        system.add_particle(Particle::new("Na"));
        system[1].charge = 1.0;
        system[1].position = Vector3D::new(1.5, 0.0, 0.0);

        return system;
    }

    #[test]
    #[should_panic]
    fn infinite_cell() {
        let mut system = testing_system();
        system.set_cell(UnitCell::new());
        let mut ewald = Ewald::new(8.0, 10);
        ewald.energy(&system);
    }

    #[test]
    #[should_panic]
    fn triclinic_cell() {
        let mut system = testing_system();
        system.set_cell(UnitCell::triclinic(0.0, 0.0, 0.0, 0.0, 0.0, 0.0));
        let mut ewald = Ewald::new(8.0, 10);
        ewald.energy(&system);
    }

    #[test]
    fn energy() {
        let system = testing_system();
        let mut ewald = Ewald::new(8.0, 10);

        let e = ewald.energy(&system);
        assert_approx_eq!(e, E_BRUTE_FORCE, 1e-4);
    }

    #[test]
    fn forces() {
        let mut system = testing_system();
        let mut ewald = Ewald::new(8.0, 10);

        let forces = ewald.forces(&system);
        let norm = (forces[0] + forces[1]).norm();
        // Total force should be null
        assert_approx_eq!(norm, 0.0, 1e-9);

        // Force is attractive
        for i in 0..3 {
            assert!(forces[0][i] > 0.0);
            assert!(forces[1][i] < 0.0);
        }

        // Finite difference computation of the force
        let e = ewald.energy(&system);
        let eps = 1e-9;
        system[0].position[0] += eps;

        let e1 = ewald.energy(&system);
        let force = ewald.forces(&system)[0][0];
        assert_approx_eq!((e - e1)/eps, force, 1e-6);
    }

    mod virial {
        use super::*;
        use types::{Vector3D, Zero};
        use potentials::GlobalPotential;

        #[test]
        fn real_space() {
            let system = testing_system();
            let ewald = Ewald::new(8.0, 10);

            // real space
            let virial = ewald.real_space_virial(&system);
            let mut forces = vec![Vector3D::zero(); 2];
            ewald.real_space_forces(&system, &mut forces);
            let w = forces[0].tensorial(&Vector3D::new(1.5, 0.0, 0.0));

            for i in 0..3 {
                for j in 0..3 {
                    assert_approx_eq!(virial[(i, j)], w[(i, j)]);
                }
            }
        }

        #[test]
        fn kspace() {
            let system = testing_system();
            let mut ewald = Ewald::new(8.0, 10);

            let virial = ewald.kspace_virial(&system);
            let mut forces = vec![Vector3D::zero(); 2];
            ewald.kspace_forces(&system, &mut forces);
            let w = forces[0].tensorial(&Vector3D::new(1.5, 0.0, 0.0));

            for i in 0..3 {
                for j in 0..3 {
                    assert_approx_eq!(virial[(i, j)], w[(i, j)]);
                }
            }
        }

        #[test]
        fn molcorrect() {
            let system = testing_system();
            let ewald = Ewald::new(8.0, 10);

            let virial = ewald.molcorrect_virial(&system);
            let mut forces = vec![Vector3D::zero(); 2];
            ewald.molcorrect_forces(&system, &mut forces);
            let w = forces[0].tensorial(&Vector3D::new(1.5, 0.0, 0.0));

            for i in 0..3 {
                for j in 0..3 {
                    assert_approx_eq!(virial[(i, j)], w[(i, j)]);
                }
            }
        }

        #[test]
        fn total() {
            let system = testing_system();
            let mut ewald = Ewald::new(8.0, 10);

            let virial = ewald.virial(&system);
            let force = ewald.forces(&system)[0];
            let w = force.tensorial(&Vector3D::new(1.5, 0.0, 0.0));

            for i in 0..3 {
                for j in 0..3 {
                    assert_approx_eq!(virial[(i, j)], w[(i, j)]);
                }
            }
        }
    }
}
