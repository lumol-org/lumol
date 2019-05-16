// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

//! This module provides some ways to initialize the velocities in a `System`
use rand_xorshift::XorShiftRng;
use rand::SeedableRng;
use rand_distr::{Normal, Uniform, Distribution};

use lumol_core::consts::K_BOLTZMANN;
use lumol_core::{System, Vector3D};

use crate::md::{Control, RemoveRotation, RemoveTranslation};

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
    dist: Normal<f64>,
    rng: XorShiftRng,
}

impl BoltzmannVelocities {
    /// Create a new `BoltzmannVelocities` at the given `temperature`.
    pub fn new(temperature: f64) -> BoltzmannVelocities {
        let dist = Normal::new(0.0, f64::sqrt(K_BOLTZMANN * temperature))
                          .expect("bad normal distribution");
        BoltzmannVelocities {
            temperature: temperature,
            dist: dist,
            rng: XorShiftRng::from_seed([
                0xeb, 0xa8, 0xe4, 0x29, 0xca, 0x60, 0x44, 0xb0,
                0xd3, 0x77, 0xc6, 0xa0, 0x21, 0x71, 0x37, 0xf7,
            ]),
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
        let b1 = ((seed >> 56) & 0xff) as u8;
        let b2 = ((seed >> 48) & 0xff) as u8;
        let b3 = ((seed >> 40) & 0xff) as u8;
        let b4 = ((seed >> 32) & 0xff) as u8;
        let b5 = ((seed >> 24) & 0xff) as u8;
        let b6 = ((seed >> 16) & 0xff) as u8;
        let b7 = ((seed >> 8) & 0xff) as u8;
        let b8 = (seed & 0xff) as u8;
        let seed = [
            b1, 0xa8, b2, 0x29, b3, 0x60, b4, 0xb0, b5, 0x77, b6, 0xa0, b7, 0x71, b8, 0xf7,
        ];
        self.rng = XorShiftRng::from_seed(seed);
    }
}

/// Initialize the velocities from an uniform distribution.
pub struct UniformVelocities {
    temperature: f64,
    dist: Uniform<f64>,
    rng: XorShiftRng,
}

impl UniformVelocities {
    /// Create a new `UniformVelocities` at the given `temperature`.
    pub fn new(temperature: f64) -> UniformVelocities {
        let factor = f64::sqrt(3.0 * K_BOLTZMANN * temperature);
        UniformVelocities {
            temperature: temperature,
            dist: Uniform::new(-factor, factor),
            rng: XorShiftRng::from_seed([
                0xeb, 0xa8, 0xe4, 0x29, 0xca, 0x60, 0x44, 0xb0,
                0xd3, 0x77, 0xc6, 0xa0, 0x21, 0x71, 0x37, 0xf7,
            ]),
        }
    }
}

impl InitVelocities for UniformVelocities {
    fn init(&mut self, system: &mut System) {
        for particle in system.particles_mut() {
            let m_inv = 1.0 / (*particle.mass);
            *particle.velocity = f64::sqrt(m_inv) * Vector3D::new(
                self.dist.sample(&mut self.rng),
                self.dist.sample(&mut self.rng),
                self.dist.sample(&mut self.rng),
            );
        }
        RemoveTranslation.control(system);
        RemoveRotation.control(system);
        scale(system, self.temperature);
    }

    fn seed(&mut self, seed: u64) {
        let b1 = ((seed >> 56) & 0xff) as u8;
        let b2 = ((seed >> 48) & 0xff) as u8;
        let b3 = ((seed >> 40) & 0xff) as u8;
        let b4 = ((seed >> 32) & 0xff) as u8;
        let b5 = ((seed >> 24) & 0xff) as u8;
        let b6 = ((seed >> 16) & 0xff) as u8;
        let b7 = ((seed >> 8) & 0xff) as u8;
        let b8 = (seed & 0xff) as u8;
        let seed = [
            b1, 0xa8, b2, 0x29, b3, 0x60, b4, 0xb0, b5, 0x77, b6, 0xa0, b7, 0x71, b8, 0xf7,
        ];
        self.rng = XorShiftRng::from_seed(seed);
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use rand::random;
    use lumol_core::{Molecule, Particle, System, Vector3D, UnitCell};

    use approx::assert_ulps_eq;
    use soa_derive::soa_zip;

    fn testing_system() -> System {
        let mut system = System::new();
        for _ in 0..10000 {
            let mut particle = Particle::new("F");
            particle.position = Vector3D::new(
                random::<f64>() * 100.0,
                random::<f64>() * 100.0,
                random::<f64>() * 100.0,
            );
            system.add_molecule(Molecule::new(particle));
        }
        return system;
    }

    fn global_translation(system: &System) -> f64 {
        let mut total_mass = 0.0;
        let mut com_velocity = Vector3D::zero();
        for (&mass, velocity) in soa_zip!(system.particles(), [mass, velocity]) {
            total_mass += mass;
            com_velocity += velocity * mass;
        }
        com_velocity /= total_mass;
        return com_velocity.norm();
    }

    #[test]
    fn init_boltzmann() {
        let mut system = testing_system();
        let mut velocities = BoltzmannVelocities::new(300.0);
        velocities.seed(1234);
        velocities.init(&mut system);
        let temperature = system.temperature();
        assert_ulps_eq!(temperature, 300.0, epsilon = 1e-9);
        assert_ulps_eq!(global_translation(&system), 0.0);
    }

    #[test]
    fn init_uniform() {
        let mut system = testing_system();
        let mut velocities = UniformVelocities::new(300.0);
        velocities.seed(1234);
        velocities.init(&mut system);
        let temperature = system.temperature();
        assert_ulps_eq!(temperature, 300.0, epsilon = 1e-9);
        assert_ulps_eq!(global_translation(&system), 0.0);
    }

    #[test]
    fn scaling_keeps_global_velocity() {
        let mut system = System::with_cell(UnitCell::cubic(10.0));
        system.add_molecule(Molecule::new(Particle::with_position("Ag", [0.0, 0.0, 0.0].into())));
        system.add_molecule(Molecule::new(Particle::with_position("Ag", [1.0, 1.0, 1.0].into())));

        system.particles_mut().velocity[0] = [1.0, 0.0, 0.0].into();
        system.particles_mut().velocity[1] = [2.0, 0.0, 0.0].into();
        assert!(global_translation(&system) > 1.0);

        scale(&mut system, 452.0);
        let temperature = system.temperature();
        assert_ulps_eq!(temperature, 452.0, epsilon = 1e-9);

        assert!(global_translation(&system) > 1e-5);
    }
}
