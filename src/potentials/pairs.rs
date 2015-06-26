/*
 * Cymbalum, Molecular Simulation in Rust
 * Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
*/

use ::types::Vector3D;
use ::types::Matrix3;

pub trait PairPotential {
    fn energy(&self, r: f64) -> f64;
    fn force(&self, r: f64) -> f64;

    fn virial(&self, r: &Vector3D) -> Matrix3 {
        let fact = self.force(r.norm());
        let rn = r.normalized();
        let force = fact * rn;
        force.tensorial(r)
    }
}

/******************************************************************************/

#[derive(Clone, Copy)]
pub struct LennardJones {
    pub sigma: f64,
    pub epsilon: f64,
}

impl PairPotential for LennardJones {
    fn energy(&self, r: f64) -> f64 {
        let s6 = f64::powi(self.sigma/r, 6);
        4.0 * self.epsilon * (f64::powi(s6, 2) - s6)
    }

    fn force(&self, r: f64) -> f64 {
        let s6 = f64::powi(self.sigma/r, 6);
        -24.0 * self.epsilon * (s6 - 2.0 * f64::powi(s6, 2)) / r
    }
}

/******************************************************************************/

#[derive(Clone, Copy)]
pub struct Harmonic {
    pub k: f64,
    pub r0: f64,
}

impl PairPotential for Harmonic {
    fn energy(&self, r: f64) -> f64 {
        let dr = r - self.r0;
        0.5 * self.k * dr * dr
    }

    fn force(&self, r: f64) -> f64 {
        self.k * (self.r0 - r)
    }
}

/******************************************************************************/

#[cfg(test)]
mod tests {
    use super::*;

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
        let harm = Harmonic{k: 50.0, r0: 2.0};
        assert_eq!(harm.energy(2.0), 0.0);
        assert_eq!(harm.energy(2.5), 6.25);
    }

    #[test]
    fn force_harmonic() {
        let harm = Harmonic{k: 50.0, r0: 2.0};
        assert_eq!(harm.force(2.0), 0.0);
        assert_eq!(harm.force(2.5), -25.0);
    }
}
