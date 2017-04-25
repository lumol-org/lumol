// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 Lumol's contributors — BSD license

use special::Error;
use std::f64::consts::PI;
use std::sync::Mutex;

use rayon::prelude::*;
use thread_local::CachedThreadLocal;

use sys::System;
use types::{Matrix3, Vector3D, Zero};
use consts::ELCC;
use energy::{PairRestriction, RestrictionInfo};

use super::{GlobalPotential, CoulombicPotential, GlobalCache};

/// Wolf summation for coulombic interactions.
///
/// This is a fast, direct, pairwise summation for coulombic potential
/// [Wolf1999].
///
/// # Examples
///
/// ```
/// use lumol::energy::Wolf;
/// use lumol::units;
///
/// // A relatively large cutoff is needed for Wolf summation
/// let wolf = Wolf::new(12.0);
///
/// use lumol::sys::System;
/// use lumol::sys::Particle;
/// use lumol::sys::UnitCell;
/// use lumol::types::Vector3D;
///
/// // Setup a system containing a NaCl pair
/// let mut system = System::new();
/// system.set_cell(UnitCell::cubic(10.0));
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
/// // Use Wolf summation for electrostatic interactions
/// system.interactions_mut().set_coulomb(Box::new(wolf));
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
        let alpha = PI/cutoff;
        let e_cst = f64::erfc(alpha*cutoff)/cutoff;
        let f_cst = f64::erfc(alpha*cutoff)/(cutoff*cutoff) + 2.0*alpha/f64::sqrt(PI) * f64::exp(-alpha*alpha*cutoff*cutoff)/cutoff;
        Wolf{
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
    fn energy_pair(&self, info: RestrictionInfo, qi: f64, qj: f64, rij: f64) -> f64 {
        if rij > self.cutoff || info.excluded {
            return 0.0;
        }
        info.scaling * qi * qj * (f64::erfc(self.alpha*rij)/rij - self.energy_cst) / ELCC
    }

    /// Compute the energy for self interaction of a particle with charge `qi`
    #[inline]
    fn energy_self(&self, qi: f64) -> f64 {
        qi * qi * (self.energy_cst/2.0 + self.alpha/f64::sqrt(PI)) / ELCC
    }

    /// Compute the force for self the pair of particles with charge `qi` and
    /// `qj`, at the distance of `rij`. The `scaling` parameter comes from the
    /// restriction associated with this potential.
    #[inline]
    fn force_pair(&self, info: RestrictionInfo, qi: f64, qj: f64, rij: Vector3D) -> Vector3D {
        let d = rij.norm();
        if d > self.cutoff || info.excluded {
            return Vector3D::zero();
        }
        let factor = f64::erfc(self.alpha*d)/(d*d) + 2.0*self.alpha/f64::sqrt(PI) * f64::exp(-self.alpha*self.alpha*d*d)/d;
        return info.scaling * qi * qj * (factor - self.force_cst) * rij.normalized() / ELCC;
    }
}

impl GlobalCache for Wolf {
    fn move_particles_cost(&self, system: &System, idxes: &[usize], newpos: &[Vector3D]) -> f64 {
        let mut e_old = 0.0;
        let mut e_new = 0.0;

        // Iterate over all interactions between a moved particle and a
        // particle not moved
        for (idx, &i) in idxes.iter().enumerate() {
            let qi = system[i].charge;
            if qi == 0.0 {continue;}
            for j in (0..system.size()).filter(|x| !idxes.contains(x)) {
                let qj = system[j].charge;
                if qj == 0.0 {continue;}

                let r_old = system.cell().distance(&system[i].position, &system[j].position);
                let r_new = system.cell().distance(&newpos[idx], &system[j].position);

                let distance = system.bond_distance(i, j);
                let info = self.restriction.information(distance);

                e_old += self.energy_pair(info, qi, qj, r_old);
                e_new += self.energy_pair(info, qi, qj, r_new);
            }
        }

        // Iterate over all interactions between two moved particles
        for (idx, &i) in idxes.iter().enumerate() {
            let qi = system[i].charge;
            if qi == 0.0 {continue;}
            for (jdx, &j) in idxes.iter().enumerate().skip(idx + 1) {
                let qj = system[j].charge;
                if qj == 0.0 {continue;}

                let r_old = system.distance(i, j);
                let r_new = system.cell().distance(&newpos[idx], &newpos[jdx]);

                let distance = system.bond_distance(i, j);
                let info = self.restriction.information(distance);

                e_old += self.energy_pair(info, qi, qj, r_old);
                e_new += self.energy_pair(info, qi, qj, r_new);
            }
        }

        return e_new - e_old;
    }

    fn update(&self) {
        // Nothing to do
    }
}

impl GlobalPotential for Wolf {
    fn cutoff(&self) -> Option<f64> {
        Some(self.cutoff)
    }

    fn energy(&self, system: &System) -> f64 {
        let natoms = system.size();
        let mut res = 0.0;
        for i in 0..natoms {
            let qi = system[i].charge;
            if qi == 0.0 {continue}
            for j in i+1..natoms {
                let qj = system[j].charge;
                if qj == 0.0 {continue}

                let distance = system.bond_distance(i, j);
                let info = self.restriction.information(distance);

                let rij = system.distance(i, j);
                res += self.energy_pair(info, qi, qj, rij);
            }
            // Remove self term
            res -= self.energy_self(qi);
        }
        return res;
    }

    fn forces(&self, system: &System) -> Vec<Vector3D> {
        // To avoid race conditions, each thread needs its
        // own local forces Vec
        let thread_forces_mutex : CachedThreadLocal<Mutex<Vec<Vector3D>>> = CachedThreadLocal::new();
        let natoms = system.size();

        (0..natoms).into_par_iter().for_each(|i| {

            /// Get the thread local forces Vec
            let mut thread_forces = thread_forces_mutex
                .get_or(|| Box::new(Mutex::new(vec![Vector3D::zero() ; natoms])))
                .lock().unwrap();

            let qi = system[i].charge;
            if qi == 0.0 { return; }
            for j in i+1..natoms {
                let qj = system[j].charge;
                if qj == 0.0 {continue}

                let distance = system.bond_distance(i, j);
                let info = self.restriction.information(distance);

                let rij = system.nearest_image(i, j);
                let force = self.force_pair(info, qi, qj, rij);
                thread_forces[i] += force;
                thread_forces[j] -= force;
            }
        });

        // At this point all the forces are computed, but the
        // results are scattered across all thread local Vecs,
        // here we gather them.
        let mut forces = vec![Vector3D::zero(); natoms];
        for mutex in thread_forces_mutex {
            let thread_forces = mutex.into_inner().unwrap();
            for (f, thread_f) in forces.iter_mut().zip(thread_forces) {
                *f += thread_f
            }
        }

        return forces;
    }

    fn virial(&self, system: &System) -> Matrix3 {
        let natoms = system.size();
        let mut res = Matrix3::zero();
        for i in 0..natoms {
            let qi = system[i].charge;
            if qi == 0.0 {continue}
            for j in i+1..natoms {
                let qj = system[j].charge;
                if qj == 0.0 {continue}

                let distance = system.bond_distance(i, j);
                let info = self.restriction.information(distance);

                let rij = system.nearest_image(i, j);
                let force = self.force_pair(info, qi, qj, rij);
                res += force.tensorial(&rij);
            }
        }
        return res;
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
    use sys::{System, UnitCell, Particle};
    use types::{Vector3D, Zero};
    use energy::GlobalPotential;

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
    fn energy() {
        let system = testing_system();
        let wolf = Wolf::new(8.0);

        let e = wolf.energy(&system);
        // Wolf is not very good for heterogeneous systems
        assert_ulps_eq!(e, E_BRUTE_FORCE, epsilon=1e-2);
    }

    #[test]
    fn forces() {
        let mut system = testing_system();
        let wolf = Wolf::new(8.0);

        let forces = wolf.forces(&system);
        let norm = (forces[0] + forces[1]).norm();
        // Total force should be null
        assert_ulps_eq!(norm, 0.0);

        // Finite difference computation of the force
        let e = wolf.energy(&system);
        let eps = 1e-9;
        system[0].position[0] += eps;

        let e1 = wolf.energy(&system);
        let force = wolf.forces(&system)[0][0];
        assert_relative_eq!((e - e1) / eps, force, epsilon=1e-6);
    }
}
