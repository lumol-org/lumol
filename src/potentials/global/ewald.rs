/* Cymbalum, Molecular Simulation in Rust - Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
 */
//! Ewald summation of coulombic potential
extern crate special;
use self::special::{erfc, erf};

use std::f64::consts::{PI, FRAC_2_SQRT_PI};
use std::f64;
use std::cell::{Cell, RefCell};

use universe::{Universe, UnitCell, CellType};
use types::{Matrix3, Vector3D, Array3, Complex};
use constants::ELCC;
use potentials::PairRestriction;

use super::{GlobalPotential, CoulombicPotential};

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
#[derive(Debug)]
pub struct Ewald {
    /// Splitting parameter between k-space and real space
    alpha: f64,
    /// Cutoff radius in real space
    rc: f64,
    /// Number of points to use in k-space
    kmax: usize,
    /// Cutoff in k-space
    kmax2: Cell<f64>,
    /// Restriction scheme
    restriction: PairRestriction,
    /// Caching exponential factors exp(-k^2 / (4 alpha^2)) / k^2
    expfactors: RefCell<Array3<f64>>,
    /// Phases for the Fourier transform
    fourier_phases: RefCell<Array3<f64>>,
    /// Fourier transform of the electrostatic density
    rho: RefCell<Array3<Complex>>,
    /// Guard for cache invalidation of expfactors
    previous_cell: RefCell<UnitCell>,
}

impl Ewald {
    /// Create an Ewald summation using the `rc` cutoff radius in real space,
    /// and `kmax` points in k-space (Fourier space).
    pub fn new(rc: f64, kmax: usize) -> Ewald {
        let kmax = kmax + 1;
        let expfactors = Array3::with_size((kmax, kmax, kmax));
        let rho = Array3::with_size((kmax, kmax, kmax));
        Ewald {
            alpha: 3.0 * PI / (rc * 4.0),
            rc: rc,
            kmax: kmax,
            kmax2: Cell::new(0.0),
            restriction: PairRestriction::None,
            expfactors: RefCell::new(expfactors),
            fourier_phases: RefCell::new(Array3::new()),
            rho: RefCell::new(rho),
            previous_cell: RefCell::new(UnitCell::new()),
        }
    }

    fn precompute(&self, cell: &UnitCell) {
        if *cell == *self.previous_cell.borrow() {
            // Do not recompute
            return;
        }
        match cell.celltype() {
            CellType::Infinite => {
                error!("Can not use Ewald sum with Infinite cell.");
                panic!();
            },
            CellType::Triclinic => {
                error!("Can not (yet) use Ewald sum with Triclinic cell.");
                unimplemented!();
            },
            CellType::Orthorombic => {
                // All good!
            },
        }
        *self.previous_cell.borrow_mut() = *cell;
        let mut expfactors = self.expfactors.borrow_mut();

        // Because we do a spherical truncation in k space, we have to transform
        // kmax into a spherical cutoff 'radius'
        let lenghts = cell.lengths();
        let max_lenght = f64::max(f64::max(lenghts.0, lenghts.1), lenghts.2);
        let min_lenght = f64::min(f64::min(lenghts.0, lenghts.1), lenghts.2);
        let k_rc = self.kmax as f64 * (2.0 * PI / max_lenght);
        self.kmax2.set(k_rc * k_rc);

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
                    if k2 > self.kmax2.get() {
                        expfactors[(ikx, iky, ikz)] = 0.0;
                        continue;
                    }
                    expfactors[(ikx, iky, ikz)] = f64::exp(-k2 / (4.0 * self.alpha * self.alpha)) / k2;
                    if ikx != 0 {expfactors[(ikx, iky, ikz)] *= 2.0;}
                    if iky != 0 {expfactors[(ikx, iky, ikz)] *= 2.0;}
                    if ikz != 0 {expfactors[(ikx, iky, ikz)] *= 2.0;}
                }
            }
        }
        expfactors[(0, 0, 0)] = 0.0;
    }
}

/// Real space part of the summation
impl Ewald {
    /// Real space contribution to the energy
    fn real_space_energy(&self, universe: &Universe) -> f64 {
        let natoms = universe.size();
        let mut energy = 0.0;
        for i in 0..natoms {
            for j in (i + 1)..natoms {
                if self.restriction.is_excluded_pair(universe, i, j) {continue}
                let s = self.restriction.scaling(universe, i, j);
                assert!(s == 1.0, "Scaling restriction scheme using Ewald are not implemented");

                let r = universe.distance(i, j);
                if r > self.rc {continue};

                energy += s * universe[i].charge * universe[j].charge * erfc(self.alpha * r) / r;
            }
        }
        return energy / ELCC;
    }

    /// Real space contribution to the forces
    fn real_space_forces(&self, universe: &Universe, res: &mut Vec<Vector3D>) {
        let natoms = universe.size();
        assert!(res.len() == universe.size());

        for i in 0..natoms {
            for j in (i + 1)..natoms {
                if self.restriction.is_excluded_pair(universe, i, j) {continue}
                let s = self.restriction.scaling(universe, i, j);
                assert!(s == 1.0, "Scaling restriction scheme using Ewald are not implemented");

                let rij = universe.wraped_vector(i, j);
                let r = rij.norm();
                if r > self.rc {continue};

                let factor = s * self.real_space_force_factor(r, universe[i].charge, universe[j].charge);
                let force = factor * rij;
                res[i] = res[i] + force;
                res[j] = res[j] - force;
            }
        }
    }

    /// Get the real-space force factor at distance `r` for charges `qi` and `qj`
    #[inline]
    fn real_space_force_factor(&self, r: f64, qi: f64, qj: f64) -> f64 {
        let mut factor = erfc(self.alpha * r) / r;
        factor += self.alpha * FRAC_2_SQRT_PI * f64::exp(-self.alpha * self.alpha * r * r);
        factor *= qi * qj / (r * r) / ELCC;
        return factor;
    }

    /// Real space contribution to the virial
    fn real_space_virial(&self, universe: &Universe) -> Matrix3 {
        let natoms = universe.size();
        let mut res = Matrix3::zero();
        for i in 0..natoms {
            for j in (i + 1)..natoms {
                if self.restriction.is_excluded_pair(universe, i, j) {continue}
                let s = self.restriction.scaling(universe, i, j);
                assert!(s == 1.0, "Scaling restriction scheme using Ewald are not implemented");

                let rij = universe.wraped_vector(i, j);
                let r = rij.norm();
                if r > self.rc {continue};

                let factor = s * self.real_space_force_factor(r, universe[i].charge, universe[j].charge);
                let force = -factor * rij;

                res = res + force.tensorial(&rij);
            }
        }
        return res;
    }
}

/// Self-interaction corection
impl Ewald {
    /// Self-interaction contribution to the energy
    fn self_energy(&self, universe: &Universe) -> f64 {
        let mut q2 = 0.0;
        for i in 0..universe.size() {
            q2 += universe[i].charge * universe[i].charge;
        }
        return -self.alpha / f64::sqrt(PI) * q2 / ELCC;
    }
}

/// k-space part of the summation
impl Ewald {
    /// Compute the Fourier transform of the electrostatic density
    fn density_fft(&self, universe: &Universe) {
        let natoms = universe.size();
        let mut fourier_phases = self.fourier_phases.borrow_mut();
        fourier_phases.resize((self.kmax, natoms, 3));

        // Do the k=0, 1 cases first
        for i in 0..natoms {
            let ri = universe.cell().fractional(&universe[i].position);
            for j in 0..3 {
                fourier_phases[(0, i, j)] = 0.0;
                fourier_phases[(1, i, j)] = -2.0 * PI * ri[j];
            }
        }

        // Use recursive definition for computing the factor for all the other values of k.
        for k in 2..self.kmax {
            for i in 0..natoms {
                for j in 0..3 {
                    fourier_phases[(k, i, j)] = fourier_phases[(k - 1, i, j)] + fourier_phases[(1, i, j)];
                }
            }
        }

        let mut rho = self.rho.borrow_mut();
        for ikx in 0..self.kmax {
            for iky in 0..self.kmax {
                for ikz in 0..self.kmax {
                    rho[(ikx, iky, ikz)] = Complex::polar(0.0, 0.0);
                    for j in 0..natoms {
                        let phi = fourier_phases[(ikx, j, 0)] + fourier_phases[(iky, j, 1)] + fourier_phases[(ikz, j, 2)];
                        rho[(ikx, iky, ikz)] = rho[(ikx, iky, ikz)] + Complex::polar(universe[j].charge, phi);
                    }
                }
            }
        }
    }

    /// k-space contribution to the energy
    fn kspace_energy(&self, universe: &Universe) -> f64 {
        self.density_fft(universe);
        let mut energy = 0.0;

        let expfactors = self.expfactors.borrow();
        let rho = self.rho.borrow();
        for ikx in 0..self.kmax {
            for iky in 0..self.kmax {
                for ikz in 0..self.kmax {
                    // The k = 0 case and the cutoff in k-space are already
                    // handled in expfactors
                    if expfactors[(ikx, iky, ikz)].abs() < f64::MIN {continue}
                    let density = rho[(ikx, iky, ikz)].norm();
                    energy += expfactors[(ikx, iky, ikz)] * density * density;
                }
            }
        }
        energy *= 2.0 * PI / (universe.cell().volume() * ELCC);
        return energy;
    }

    /// k-space contribution to the forces
    fn kspace_forces(&self, universe: &Universe, res: &mut Vec<Vector3D>) {
        assert!(res.len() == universe.size());
        self.density_fft(universe);

        let factor = 4.0 * PI / (universe.cell().volume() * ELCC);
        let (rec_kx, rec_ky, rec_kz) = universe.cell().reciprocal_vectors();

        let expfactors = self.expfactors.borrow();
        for ikx in 0..self.kmax {
            for iky in 0..self.kmax {
                for ikz in 0..self.kmax {
                    // The k = 0 and the cutoff in k-space are already handled in
                    // expfactors.
                    if expfactors[(ikx, iky, ikz)].abs() < f64::MIN {continue}
                    let k = (ikx as f64) * rec_kx + (iky as f64) * rec_ky + (ikz as f64) * rec_kz;
                    for i in 0..universe.size() {
                        let qi = universe[i].charge;
                        for j in (i + 1)..universe.size() {
                            let qj = universe[j].charge;
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
    #[inline]
    fn kspace_force_factor(&self, i: usize, j: usize, ikx: usize, iky: usize, ikz: usize, qi: f64, qj: f64) -> f64 {
        let fourier_phases = self.fourier_phases.borrow();
        let expfactors = self.expfactors.borrow();

        let fourier_i = fourier_phases[(ikx, i, 0)] + fourier_phases[(iky, i, 1)] + fourier_phases[(ikz, i, 2)];
        let fourier_j = fourier_phases[(ikx, j, 0)] + fourier_phases[(iky, j, 1)] + fourier_phases[(ikz, j, 2)];
        let sin_kr = fast_sin(fourier_i - fourier_j);

        return qi * qj * expfactors[(ikx, iky, ikz)] * sin_kr;
    }

    /// k-space contribution to the virial
    fn kspace_virial(&self, universe: &Universe) -> Matrix3 {
        self.density_fft(universe);
        let mut res = Matrix3::zero();

        let factor = 4.0 * PI / (universe.cell().volume() * ELCC);
        let (rec_kx, rec_ky, rec_kz) = universe.cell().reciprocal_vectors();

        let expfactors = self.expfactors.borrow();
        for ikx in 0..self.kmax {
            for iky in 0..self.kmax {
                for ikz in 0..self.kmax {
                    // The k = 0 and the cutoff in k-space are already handled in
                    // expfactors.
                    if expfactors[(ikx, iky, ikz)].abs() < f64::MIN {continue}
                    let k = (ikx as f64) * rec_kx + (iky as f64) * rec_ky + (ikz as f64) * rec_kz;
                    for i in 0..universe.size() {
                        let qi = universe[i].charge;
                        for j in (i + 1)..universe.size() {
                            let qj = universe[j].charge;
                            let force = factor * self.kspace_force_factor(i, j, ikx, iky, ikz, qi, qj) * k;
                            let rij = universe.wraped_vector(i, j);

                            res = res + force.tensorial(&rij);
                        }
                    }
                }
            }
        }
        return res;
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
    /// Molecular correction contribution to the energy
    fn molcorrect_energy(&self, universe: &Universe) -> f64 {
        let natoms = universe.size();
        let mut energy = 0.0;

        for i in 0..natoms {
            let qi = universe[i].charge;
            // I can not manage to get this work with a loop from (i+1) to N. The finite
            // difference test (testing that the force is the same that the finite difference
            // of the energy) always fail. So let's use it that way for now.
            for j in 0..natoms {
                if i == j {continue}
                // Only account for excluded pairs
                if !self.restriction.is_excluded_pair(universe, i, j) {continue}
                let s = self.restriction.scaling(universe, i, j);

                let qj = universe[j].charge;
                let r = universe.distance(i, j);
                assert!(r < self.rc, "Atoms in molecule are separated by more than the cutoff radius of Ewald sum.");

                energy += 0.5 * qi * qj * s / ELCC * erf(self.alpha * r)/r;
            }
        }
        return energy;
    }

    /// Molecular correction contribution to the forces
    fn molcorrect_forces(&self, universe: &Universe, res: &mut Vec<Vector3D>) {
        let natoms = universe.size();
        assert!(res.len() == natoms);

        for i in 0..natoms {
            let qi = universe[i].charge;
            for j in 0..natoms {
                if i == j {continue}
                // Only account for excluded pairs
                if !self.restriction.is_excluded_pair(universe, i, j) {continue}
                let s = self.restriction.scaling(universe, i, j);

                let qj = universe[j].charge;
                let rij = universe.wraped_vector(i, j);
                let r = rij.norm();
                assert!(r < self.rc, "Atoms in molecule are separated by more than the cutoff radius of Ewald sum.");

                let factor = s * self.molcorrect_force_factor(qi, qj, r);
                res[i] = res[i] - factor * rij;
            }
        }
    }

    /// Get the force factor for particles with charges `qi` and `qj`, at
    /// distance `r`.
    #[inline]
    fn molcorrect_force_factor(&self, qi: f64, qj: f64, r: f64) -> f64 {
        qi * qj / (ELCC * r * r) * (2.0 * self.alpha / f64::sqrt(PI) * f64::exp(-self.alpha * self.alpha * r * r) - erf(self.alpha * r) / r)
    }

    /// Molecular correction contribution to the virial
    fn molcorrect_virial(&self, universe: &Universe) -> Matrix3 {
        let natoms = universe.size();
        let mut res = Matrix3::zero();

        for i in 0..natoms {
            let qi = universe[i].charge;
            for j in 0..natoms {
                if i == j {continue}
                // Only account for excluded pairs
                if !self.restriction.is_excluded_pair(universe, i, j) {continue}
                let s = self.restriction.scaling(universe, i, j);

                let qj = universe[j].charge;
                let rij = universe.wraped_vector(i, j);
                let r = rij.norm();
                assert!(r < self.rc, "Atoms in molecule are separated by more than the cutoff radius of Ewald sum.");

                let force = s * self.molcorrect_force_factor(qi, qj, r) * rij;
                res = res + force.tensorial(&rij);
            }
        }
        return res;
    }
}

impl GlobalPotential for Ewald {
    fn energy(&self, universe: &Universe) -> f64 {
        self.precompute(universe.cell());
        let real = self.real_space_energy(universe);
        let self_e = self.self_energy(universe);
        let kspace = self.kspace_energy(universe);
        let molecular = self.molcorrect_energy(universe);
        return real + self_e + kspace + molecular;
    }

    fn forces(&self, universe: &Universe) -> Vec<Vector3D> {
        self.precompute(universe.cell());
        let mut res = vec![Vector3D::new(0.0, 0.0, 0.0); universe.size()];
        self.real_space_forces(universe, &mut res);
        /* No self force */
        self.kspace_forces(universe, &mut res);
        self.molcorrect_forces(universe, &mut res);
        return res;
    }

    fn virial(&self, universe: &Universe) -> Matrix3 {
        self.precompute(universe.cell());
        let real = self.real_space_virial(universe);
        /* No self virial */
        let kspace = self.kspace_virial(universe);
        let molecular = self.molcorrect_virial(universe);
        println!("{:?} {:?} {:?}", real, kspace, molecular);
        return real + kspace + molecular;
    }
}

impl CoulombicPotential for Ewald {
    fn set_restriction(&mut self, restriction: PairRestriction) {
        self.restriction = restriction;
    }
}


#[cfg(test)]
mod tests {
    pub use super::*;
    use universe::{Universe, UnitCell, Particle};
    use types::Vector3D;
    use potentials::GlobalPotential;

    const E_BRUTE_FORCE: f64 = -0.09262397663346732;

    pub fn testing_universe() -> Universe {
        let mut universe = Universe::from_cell(UnitCell::cubic(20.0));

        universe.add_particle(Particle::new("Cl"));
        universe[0].charge = -1.0;
        universe[0].position = Vector3D::new(0.0, 0.0, 0.0);

        universe.add_particle(Particle::new("Na"));
        universe[1].charge = 1.0;
        universe[1].position = Vector3D::new(1.5, 0.0, 0.0);

        return universe;
    }

    #[test]
    #[should_panic]
    fn infinite_cell() {
        let mut universe = testing_universe();
        universe.set_cell(UnitCell::new());
        let ewald = Ewald::new(8.0, 10);
        ewald.energy(&universe);
    }

    #[test]
    #[should_panic]
    fn triclinic_cell() {
        let mut universe = testing_universe();
        universe.set_cell(UnitCell::triclinic(0.0, 0.0, 0.0, 0.0, 0.0, 0.0));
        let ewald = Ewald::new(8.0, 10);
        ewald.energy(&universe);
    }

    #[test]
    fn energy() {
        let universe = testing_universe();
        let ewald = Ewald::new(8.0, 10);

        let e = ewald.energy(&universe);
        assert_approx_eq!(e, E_BRUTE_FORCE, 1e-4);
    }

    #[test]
    fn forces() {
        let mut universe = testing_universe();
        let ewald = Ewald::new(8.0, 10);

        let forces = ewald.forces(&universe);
        let norm = (forces[0] + forces[1]).norm();
        // Total force should be null
        assert_approx_eq!(norm, 0.0, 1e-9);

        // Finite difference computation of the force
        let e = ewald.energy(&universe);
        let eps = 1e-9;
        universe[0].position.x += eps;

        let e1 = ewald.energy(&universe);
        let force = ewald.forces(&universe)[0].x;
        assert_approx_eq!((e - e1)/eps, force, 1e-6);
    }

    #[test]
    fn virial() {
        let universe = testing_universe();
        let ewald = Ewald::new(8.0, 10);

        let virial = ewald.virial(&universe);

        let force = ewald.forces(&universe)[0];
        let w = force.tensorial(&Vector3D::new(1.5, 0.0, 0.0));

        for i in 0..3 {
            for j in 0..3 {
                assert_approx_eq!(virial[(i, j)], w[(i, j)]);
            }
        }
    }
}
