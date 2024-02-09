// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license
#![allow(clippy::iter_without_into_iter)]

use std::fmt;
use soa_derive::StructOfArray;

use crate::sys::get_atomic_mass;
use crate::Vector3D;

/// A particle kind. Particles with the same name will have the same kind. This
/// is used for faster potential lookup.
#[derive(Clone, Copy, Hash, PartialOrd, Ord, PartialEq, Eq, Debug)]
pub struct ParticleKind(pub u32);

impl ParticleKind {
    /// Get an invalid value (`u32::max_value()`) to use as a marker
    pub fn invalid() -> ParticleKind {
        ParticleKind(u32::max_value())
    }
}

impl fmt::Display for ParticleKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

/// The Particle type hold basic data about a particle in the system. It is self
/// contained, so that it will be easy to send data between parallels processes.
#[derive(Clone, Debug, StructOfArray)]
#[soa_derive(Clone, Debug)]
pub struct Particle {
    /// Particle name.
    pub name: String,
    /// Particle kind, an index for potentials lookup
    pub kind: ParticleKind,
    /// Particle charge
    pub charge: f64,
    /// Particle mass
    pub mass: f64,
    /// Particle positions
    pub position: Vector3D,
    /// Particle velocity, if needed
    pub velocity: Vector3D,
}

impl Particle {
    /// Create a new `Particle` from a `name`, setting the mass to the atomic
    /// mass if the `name` can be found in the periodic table. The charge,
    /// position, and velocity are set to 0.
    pub fn new<S: Into<String>>(name: S) -> Particle {
        Particle::with_position(name, Vector3D::zero())
    }

    /// Create a new `Particle` from a `name` and a `position`, setting the
    /// mass to the atomic mass if the `name` can be found in the periodic
    /// table. The charge and velocity are set to 0.
    pub fn with_position<S: Into<String>>(name: S, position: Vector3D) -> Particle {
        let name = name.into();
        let mass = get_atomic_mass(&name).unwrap_or(0.0);
        Particle {
            name: name,
            mass: mass,
            charge: 0.0,
            kind: ParticleKind::invalid(),
            position: position,
            velocity: Vector3D::zero(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Vector3D;

    #[test]
    fn mass_initialization() {
        let particle = Particle::new("O");
        assert_eq!(particle.mass, 15.999);
    }

    #[test]
    fn name() {
        let particle = Particle::new("");
        assert_eq!(particle.name, "");

        assert_eq!(particle.mass, 0.0);
        assert_eq!(particle.charge, 0.0);
        assert_eq!(particle.kind, ParticleKind::invalid());
        assert_eq!(particle.position, Vector3D::new(0.0, 0.0, 0.0));
        assert_eq!(particle.velocity, Vector3D::new(0.0, 0.0, 0.0));
    }

    #[test]
    fn with_position() {
        let particle = Particle::with_position("", Vector3D::new(1.0, 2.0, 3.0));
        assert_eq!(particle.name, "");
        assert_eq!(particle.position, Vector3D::new(1.0, 2.0, 3.0));

        assert_eq!(particle.mass, 0.0);
        assert_eq!(particle.charge, 0.0);
        assert_eq!(particle.kind, ParticleKind::invalid());
        assert_eq!(particle.velocity, Vector3D::new(0.0, 0.0, 0.0));
    }
}
