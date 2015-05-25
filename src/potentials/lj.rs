/*
 * Cymbalum, Molecular Simulation in Rust
 * Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
*/

use super::traits::PairPotential;

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

#[cfg(test)]
mod tests {
    use super::*;
    use ::potentials::traits::*;

    #[test]
    fn energy() {
        let lj = LennardJones{epsilon: 0.8, sigma: 2.0};
        assert_eq!(lj.energy(2.0), 0.0);
        assert_eq!(lj.energy(2.5), -0.6189584744448002);
    }

    #[test]
    fn force() {
        let lj = LennardJones{epsilon: 0.8, sigma: 2.0};
        assert_approx_eq!(lj.force(f64::powf(2.0, 1.0/6.0) * 2.0), 0.0, 1e-15);
        assert_approx_eq!(lj.force(2.5), -0.95773475733504, 1e-15);
    }
}
