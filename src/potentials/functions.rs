/*
 * Cymbalum, Molecular Simulation in Rust
 * Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
*/

//! Pair potentials traits and implementations.

use ::types::Vector3D;
use ::types::Matrix3;

/// A `PotentialFunction` is a set of two parametric functions which takes a
/// single scalar variable and return the energy or the norm of the force
/// corresponding to the value of that variable.
///
/// The scalar variable will be the distance for pair potentials, the angle for
/// angles or dihedral angles potentials, *etc.*
pub trait PotentialFunction {
    /// Get the energy corresponding to the variable `x`
    fn energy(&self, x: f64) -> f64;
    /// Get the force norm corresponding to the variable `x`
    fn force(&self, x: f64) -> f64;
}

/// Pair potential can be expressed by implenting the `PairPotential` trait.
pub trait PairPotential : PotentialFunction {
    /// Compute the virial contribution corresponding to the distance `r` between the particles
    fn virial(&self, r: &Vector3D) -> Matrix3 {
        let fact = self.force(r.norm());
        let rn = r.normalized();
        let force = - fact * rn;
        force.tensorial(r)
    }
}

/// Angle potential can be expressed by implenting the `AnglePotential` trait.
pub trait AnglePotential : PotentialFunction {}

/// Dihedral angles potential can be expressed by implenting the
/// `DihedralPotential` trait.
pub trait DihedralPotential : PotentialFunction {}

/******************************************************************************/
/// The `NullPotential` always returns 0.0 as energy and force.
///
/// It is to be used when there is no potential interaction between two
/// particles kinds.
#[derive(Clone, Copy)]
pub struct NullPotential;
impl PotentialFunction for NullPotential {
    #[inline] fn energy(&self, _: f64) -> f64 {0.0}
    #[inline] fn force(&self, _: f64) -> f64 {0.0}
}

impl PairPotential for NullPotential {}
impl AnglePotential for NullPotential {}
impl DihedralPotential for NullPotential {}

/******************************************************************************/
/// Lennard-Jones potential, using the following form:
/// $$ V(r) = 4 \epsilon \left[ (\sigma/r)^{12} - (\sigma/r)^{6} \right] $$
/// where $\sigma$ is a distance constant, and $\epsilon$ an energetic constant.
#[derive(Clone, Copy)]
pub struct LennardJones {
    /// Distance constant of the Lennard-Jones potential
    pub sigma: f64,
    /// Energy constant of the Lennard-Jones potential
    pub epsilon: f64,
}

impl PotentialFunction for LennardJones {
    #[inline]
    fn energy(&self, r: f64) -> f64 {
        let s6 = f64::powi(self.sigma/r, 6);
        4.0 * self.epsilon * (f64::powi(s6, 2) - s6)
    }

    #[inline]
    fn force(&self, r: f64) -> f64 {
        let s6 = f64::powi(self.sigma/r, 6);
        -24.0 * self.epsilon * (s6 - 2.0 * f64::powi(s6, 2)) / r
    }
}

impl PairPotential for LennardJones {}

/******************************************************************************/
/// Harmonic potential, using the following form:
/// $$ V(x) = \frac 12 k (x - x_0)^2 $$
/// where $x_0$ is the distance equilibrium, and $k$ the elastic constant.
#[derive(Clone, Copy)]
pub struct Harmonic {
    /// Spring constant
    pub k: f64,
    /// Equilibrium value
    pub x0: f64,
}

impl PotentialFunction for Harmonic {
    #[inline]
    fn energy(&self, x: f64) -> f64 {
        let dx = x - self.x0;
        0.5 * self.k * dx * dx
    }

    #[inline]
    fn force(&self, x: f64) -> f64 {
        self.k * (self.x0 - x)
    }
}

impl PairPotential for Harmonic {}
impl AnglePotential for Harmonic {}
impl DihedralPotential for Harmonic {}

/******************************************************************************/
/// Cosine harmonic potential, using the following form:
/// $$ V(r) = \frac 12 k (\cos r - \cos x_0)^2 $$
/// where $x_0$ is the distance equilibrium, and $k$ the elastic constant.
#[derive(Clone, Copy)]
pub struct CosineHarmonic {
    /// Spring constant
    k: f64,
    /// Cosine of the equilibrium value
    cos_x0: f64,
}

impl CosineHarmonic {
    /// Create a new `CosineHarmonic` potentials, with elastic constant of `k`
    /// and equilibrium value of `x0`
    pub fn new(k: f64, x0:f64) -> CosineHarmonic {
        CosineHarmonic{k: k, cos_x0: f64::cos(x0)}
    }
}

impl PotentialFunction for CosineHarmonic {
    #[inline]
    fn energy(&self, x: f64) -> f64 {
        let dr = f64::cos(x) - self.cos_x0;
        0.5 * self.k * dr * dr
    }

    #[inline]
    fn force(&self, x: f64) -> f64 {
        self.k * (f64::cos(x) - self.cos_x0) * f64::sin(x)
    }
}

impl AnglePotential for CosineHarmonic {}
impl DihedralPotential for CosineHarmonic {}

/******************************************************************************/

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn energy_null() {
        let null = NullPotential;
        assert_eq!(null.energy(2.0), 0.0);
        assert_eq!(null.energy(2.5), 0.0);

        assert_eq!(null.force(2.0), 0.0);
        assert_eq!(null.force(2.5), 0.0);
    }

    #[test]
    fn energy_lj() {
        let lj = LennardJones{epsilon: 0.8, sigma: 2.0};
        assert_eq!(lj.energy(2.0), 0.0);
        assert_eq!(lj.energy(2.5), -0.6189584744448002);
    }

    #[test]
    fn force_lj() {
        let lj = LennardJones{epsilon: 0.8, sigma: 2.0};
        assert_approx_eq!(lj.force(f64::powf(2.0, 1.0/6.0) * 2.0), 0.0, 1e-15);
        assert_approx_eq!(lj.force(2.5), -0.95773475733504, 1e-15);
    }

    #[test]
    fn energy_harmonic() {
        let harm = Harmonic{k: 50.0, x0: 2.0};
        assert_eq!(harm.energy(2.0), 0.0);
        assert_eq!(harm.energy(2.5), 6.25);
    }

    #[test]
    fn force_harmonic() {
        let harm = Harmonic{k: 50.0, x0: 2.0};
        assert_eq!(harm.force(2.0), 0.0);
        assert_eq!(harm.force(2.5), -25.0);
    }

    #[test]
    fn energy_cosine_harmonic() {
        let harm = CosineHarmonic::new(50.0, 2.0);
        assert_eq!(harm.energy(2.0), 0.0);
        let dcos = f64::cos(2.5) - f64::cos(2.0);
        assert_eq!(harm.energy(2.5), 0.5 * 50.0 * dcos * dcos);
    }

    #[test]
    fn force_cosine_harmonic() {
        let harm = CosineHarmonic::new(50.0, 2.0);
        assert_eq!(harm.force(2.0), 0.0);
        let dcos = f64::cos(2.5) - f64::cos(2.0);
        assert_eq!(harm.force(2.5), 50.0 * dcos * f64::sin(2.5));
    }
}
