/* Cymbalum, Molecular Simulation in Rust - Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
 */

//! `Particle` type and manipulation.

use ::types::Vector3D;
use super::PeriodicTable;

/// The Particle type hold basic data about a particle in the system. It is self
/// contained, so that it will be easy to send data between parrallels
/// processes.
#[derive(Clone, Debug)]
pub struct Particle {
    /// Particle name
    name: String,
    /// Particle kind, an index for potentials lookup
    kind: u16,
    /// Particle mass
    mass: f64,
    /// Particle charge
    charge: f64,
    /// Particle positions
    position: Vector3D,
    /// Particle velocity, if needed
    velocity: Vector3D,
}


impl Particle {
    /// Create a new `Particle` from a name
    pub fn new<S: Into<String>>(n: S) -> Particle {
        let name = n.into();
        let mass = match PeriodicTable::mass(&*name) {
            Some(val) => val,
            None => panic!("Could not find the mass for the {} particle", name)
        };
        Particle{name: name,
                 mass: mass as f64,
                 charge: 0.0,
                 kind: u16::max_value(),
                 position: Vector3D::new(0.0, 0.0, 0.0),
                 velocity: Vector3D::new(0.0, 0.0, 0.0)}
    }

    /// Get the particle name
    #[inline] pub fn name<'a>(&'a self) -> &'a str {
        &self.name
    }

    /// Set the particle name
    #[inline] pub fn set_name<S>(&mut self, name: S) where S: Into<String> {
        self.name = name.into();
    }

    /// Get the particle kind. This kind is used for fast lookup in interactions
    /// tables
    #[inline] pub fn kind(&self) -> u16 {
        self.kind
    }

    /// Set the particle kind.
    #[inline] pub fn set_kind(&mut self, kind: u16) {
        self.kind = kind;
    }

    /// Get the particle mass
    #[inline] pub fn mass(&self) -> f64 {
        self.mass
    }

    /// Set the particle mass
    #[inline] pub fn set_mass(&mut self, mass: f64) {
        self.mass = mass;
    }

    /// Get the particle charge
    #[inline] pub fn charge(&self) -> f64 {
        self.charge
    }

    /// Set the particle charge
    #[inline] pub fn set_charge(&mut self, charge: f64) {
        self.charge = charge;
    }

    /// Get the particle position
    #[inline] pub fn position<'a>(&'a self) -> &'a Vector3D {
        &self.position
    }

    /// Set the particle position
    #[inline] pub fn set_position(&mut self, pos: Vector3D) {
        self.position = pos;
    }

    /// Add `dr` to the particle position
    #[inline] pub fn add_position(&mut self, dr: Vector3D) {
        self.position = self.position + dr;
    }

    /// Get the particle velocity
    #[inline] pub fn velocity<'a>(&'a self) -> &'a Vector3D {
        &self.velocity
    }

    /// Set the particle velocity
    #[inline] pub fn set_velocity(&mut self, vel: Vector3D) {
        self.velocity = vel;
    }

    /// Add `dv` to the particle velocity
    #[inline] pub fn add_velocity(&mut self, dv: Vector3D) {
        self.velocity = self.velocity + dv;
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use ::types::*;

    #[test]
    fn name() {
        let mut part = Particle::new("O");
        assert_eq!(part.name(), "O");

        part.set_name("H");
        assert_eq!(part.name(), "H");

        let name = "Zn";
        part.set_name(name);
        assert_eq!(part.name(), "Zn");
    }

    #[test]
    fn mass() {
        let mut part = Particle::new("O");
        assert_eq!(part.mass(), 15.999f32 as f64);
        part.set_mass(10.0);
        assert_eq!(part.mass(), 10.0);
    }

    #[test]
    fn charge() {
        let mut part = Particle::new("O");
        assert_eq!(part.charge(), 0.0);
        part.set_charge(-1.0);
        assert_eq!(part.charge(), -1.0);
    }

    #[test]
    fn kind() {
        let mut part = Particle::new("O");
        part.set_kind(42);
        assert_eq!(part.kind(), 42);
    }

    #[test]
    fn coordinates() {
        let mut part = Particle::new("O");
        assert_eq!(part.position(), &Vector3D::new(0.0, 0.0, 0.0));

        part.set_position(Vector3D::new(1.0, 2.0, 3.0));
        assert_eq!(part.position(), &Vector3D::new(1.0, 2.0, 3.0));

        part.set_velocity(Vector3D::new(1.0, 2.0, 3.0));
        assert_eq!(part.velocity(), &Vector3D::new(1.0, 2.0, 3.0));
    }
}
