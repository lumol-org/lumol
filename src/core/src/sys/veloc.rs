// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

//! This module provides some ways to initialize the velocities in a `System`
use rand::distributions::{Range, Normal, Sample};
use rand::Isaac64Rng;
use rand::SeedableRng;

use consts::K_BOLTZMANN;
use types::Vector3D;
use sys::System;
use sim::md::{RemoveRotation, RemoveTranslation, Control};

/// Scale all velocities in the `System` such that the `system` temperature
/// is `temperature`.
pub fn scale(system: &mut System, temperature: f64) {
    let instant_temperature = system.temperature();
    let factor = f64::sqrt(temperature / instant_temperature);
    for velocity in system.particles_mut().velocity {
        *velocity *= factor;
    }
}

/// A method to initialize the velocities of a system.
pub trait InitVelocities {
    /// Initialize the velocities of the system.
    fn init(&mut self, system: &mut System);
    /// Set the seed of the random number generator. The default seed is 42.
    fn seed(&mut self, seed: u64);
}

/// Initialize the velocities from a Boltzmann distribution.
pub struct BoltzmannVelocities {
    temperature: f64,
    dist: Normal,
    rng: Isaac64Rng,
}

impl BoltzmannVelocities {
    /// Create a new `BoltzmannVelocities` at the given `temperature`.
    pub fn new(temperature: f64) -> BoltzmannVelocities {
        BoltzmannVelocities{
            temperature: temperature,
            dist: Normal::new(0.0, f64::sqrt(K_BOLTZMANN * temperature)),
            rng: Isaac64Rng::from_seed(&[42]),
        }
    }
}

impl InitVelocities for BoltzmannVelocities {
    fn init(&mut self, system: &mut System) {
        for particle in system.particles_mut() {
            let m_inv = 1.0 / (*particle.mass);
            let x = f64::sqrt(m_inv) * self.dist.sample(&mut self.rng);
            let y = f64::sqrt(m_inv) * self.dist.sample(&mut self.rng);
            let z = f64::sqrt(m_inv) * self.dist.sample(&mut self.rng);
            *particle.velocity = Vector3D::new(x, y, z);
        }
        RemoveTranslation.control(system);
        RemoveRotation.control(system);
        scale(system, self.temperature);
    }

    fn seed(&mut self, seed: u64) {
        self.rng.reseed(&[seed]);
    }
}

/// Initialize the velocities from an uniform distribution.
pub struct UniformVelocities {
    temperature: f64,
    dist: Range<f64>,
    rng: Isaac64Rng,
}

impl UniformVelocities {
    /// Create a new `UniformVelocities` at the given `temperature`.
    pub fn new(temperature: f64) -> UniformVelocities {
        let factor = f64::sqrt(3.0 * K_BOLTZMANN * temperature);
        UniformVelocities{
            temperature: temperature,
            dist: Range::new(-factor, factor),
            rng: Isaac64Rng::from_seed(&[42]),
        }
    }
}

impl InitVelocities for UniformVelocities {
    fn init(&mut self, system: &mut System) {
        for particle in system.particles_mut() {
            let m_inv = 1.0 / (*particle.mass);
            let x = f64::sqrt(m_inv) * self.dist.sample(&mut self.rng);
            let y = f64::sqrt(m_inv) * self.dist.sample(&mut self.rng);
            let z = f64::sqrt(m_inv) * self.dist.sample(&mut self.rng);
            *particle.velocity = Vector3D::new(x, y, z);
        }
        RemoveTranslation.control(system);
        RemoveRotation.control(system);
        scale(system, self.temperature);
    }

    fn seed(&mut self, seed: u64) {
        self.rng.reseed(&[seed]);
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use rand::random;
    use sys::{System, Particle};
    use utils::system_from_xyz;

    fn testing_system() -> System {
        let mut system = System::new();
        for _ in 0..10000 {
            let mut particle = Particle::new("F");
            particle.position = Vector3D::new(
                random::<f64>() * 100.0, random::<f64>() * 100.0, random::<f64>() * 100.0
            );
            system.add_particle(particle);
        }
        return system;
    }

    fn global_translation(system: &System) -> f64 {
        use sys::zip_particle::{Mass, Velocity};
        use types::{Vector3D, Zero};

        let total_mass = system.particles().mass.iter().sum();
        let mut com_velocity = Vector3D::zero();
        for (&mass, velocity) in system.particles().zip((&Mass, &Velocity)) {
            com_velocity += velocity * mass / total_mass;
        }
        return com_velocity.norm();
    }

    #[test]
    fn init_boltzmann() {
        let mut system = testing_system();
        let mut velocities = BoltzmannVelocities::new(300.0);
        velocities.seed(1234);
        velocities.init(&mut system);
        let temperature = system.temperature();
        assert_ulps_eq!(temperature, 300.0, epsilon=1e-9);
        assert_ulps_eq!(global_translation(&system), 0.0);
    }

    #[test]
    fn init_uniform() {
        let mut system = testing_system();
        let mut velocities = UniformVelocities::new(300.0);
        velocities.seed(1234);
        velocities.init(&mut system);
        let temperature = system.temperature();
        assert_ulps_eq!(temperature, 300.0, epsilon=1e-9);
        assert_ulps_eq!(global_translation(&system), 0.0);
    }

    #[test]
    fn scaling_keeps_global_velocity() {
        let mut system = system_from_xyz("2
        cell: 20.0
        Ag 0 0 0 1 0 0
        Ag 1 1 1 2 0 0
        ");
        assert!(global_translation(&system) > 1.0);

        scale(&mut system, 452.0);
        let temperature = system.temperature();
        assert_ulps_eq!(temperature, 452.0, epsilon=1e-9);

        println!("{:?}", global_translation(&system));
        assert!(global_translation(&system) > 1e-5);
    }
}
