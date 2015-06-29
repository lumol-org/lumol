/*
 * Cymbalum, Molecular Simulation in Rust
 * Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
*/
extern crate rand;
use self::rand::distributions::{Range, Normal};
use self::rand::distributions::Sample;
use self::rand::Isaac64Rng;
use self::rand::SeedableRng;

use ::constants::K_BOLTZMANN;
use ::types::Vector3D;
use super::Universe;

/// Random initializer for the velocities of an universe.
pub trait InitVelocities {
    /// Initialize the velocities of the universe.
    fn init(&mut self, universe: &mut Universe);
    /// Set the seed of the random number generator. The default seed is 42.
    fn seed(&mut self, seed: u64);
}

/// Initialize the velocities from a Boltzman distribution.
pub struct BoltzmanVelocities {
    dist: Normal,
    rng: Isaac64Rng,
}

impl BoltzmanVelocities {
    /// This `BoltzmanVelocities` initializer will initialize velocities at
    /// temperature `T`.
    pub fn new(T: f64) -> BoltzmanVelocities {
        BoltzmanVelocities{
            dist: Normal::new(0.0, f64::sqrt(K_BOLTZMANN * T)),
            rng: Isaac64Rng::from_seed(&[42]),
        }
    }
}

impl InitVelocities for BoltzmanVelocities {
    fn init(&mut self, universe: &mut Universe) {
        for i in 0..universe.size() {
            let part = &mut universe[i];
            let m_inv = 1.0 / part.mass();
            let x = f64::sqrt(m_inv) * self.dist.sample(&mut self.rng);
            let y = f64::sqrt(m_inv) * self.dist.sample(&mut self.rng);
            let z = f64::sqrt(m_inv) * self.dist.sample(&mut self.rng);
            part.set_velocity(Vector3D::new(x, y, z));
        }
    }

    fn seed(&mut self, seed: u64) {
        self.rng.reseed(&[seed]);
    }
}

/// Initialize the velocities from an uniform distribution.
pub struct UniformVelocities {
    dist: Range<f64>,
    rng: Isaac64Rng,
}

impl UniformVelocities {
    /// This `UniformVelocities` initializer will initialize velocities at
    /// temperature `T`.
    pub fn new(T: f64) -> UniformVelocities {
        let factor = f64::sqrt(3.0*K_BOLTZMANN * T);
        UniformVelocities{
            dist: Range::new(-factor, factor),
            rng: Isaac64Rng::from_seed(&[42]),
        }
    }
}

impl InitVelocities for UniformVelocities {
    fn init(&mut self, universe: &mut Universe) {
        for i in 0..universe.size() {
            let part = &mut universe[i];
            let m_inv = 1.0 / part.mass();
            let x = f64::sqrt(m_inv) * self.dist.sample(&mut self.rng);
            let y = f64::sqrt(m_inv) * self.dist.sample(&mut self.rng);
            let z = f64::sqrt(m_inv) * self.dist.sample(&mut self.rng);
            part.set_velocity(Vector3D::new(x, y, z));
        }
    }

    fn seed(&mut self, seed: u64) {
        self.rng.reseed(&[seed]);
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use ::simulation::{Compute, Temperature};

    use ::universe::{Universe, Particle};

    fn testing_universe() -> Universe {
        let mut universe = Universe::new();
        for _ in 0..10000 {
            universe.add_particle(Particle::new("F"));
        }
        return universe;
    }

    #[test]
    fn init_boltzmann() {
        let mut universe = testing_universe();
        let mut velocities = BoltzmanVelocities::new(300.0);
        velocities.seed(1234);
        velocities.init(&mut universe);
        let T = Temperature.compute(&universe);
        assert!(f64::abs(T - 300.0) < 2.0);
    }

    #[test]
    fn init_uniform() {
        let mut universe = testing_universe();
        let mut velocities = UniformVelocities::new(300.0);
        velocities.seed(1234);
        velocities.init(&mut universe);
        let T = Temperature.compute(&universe);
        assert!(f64::abs(T - 300.0) < 2.0);
    }
}
