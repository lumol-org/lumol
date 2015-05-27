/*
 * Cymbalum, Molecular Simulation in Rust
 * Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
*/

use ::types::*;

/// The Particle type hold basic data about a particle in the system. It is self
/// contained, so that it will be easy to send data between parrallels
/// processes.
#[derive(Clone, Debug)]
pub struct Particle {
    /// Particle name
    pub name: String,
    /// Particle mass
    pub mass: f64,
    /// Particle positions
    pub position: Vector3D,
    /// Particle velocity, if needed
    pub velocity: Vector3D,
}


impl Particle {
    pub fn new(name: &str) -> Particle {
        // TODO: get the mass here
        Particle{name: name.to_string(),
                 mass: 1.0,
                 position: Vector3D::new(0.0, 0.0, 0.0),
                 velocity: Vector3D::new(0.0, 0.0, 0.0)}
    }

    pub fn mass(&mut self, mass: f64) {
        (*self).mass = mass;
    }

    pub fn position(&mut self, pos: Vector3D) {
        (*self).position = pos;
    }

    pub fn velocity(&mut self, vel: Vector3D) {
        (*self).velocity = vel;
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use ::types::*;

    #[test]
    fn set_mass() {
        let mut part = Particle::new("O");
        assert_eq!(part.name, "O");

        part.mass(16.0);
        assert_eq!(part.mass, 16.0);
    }

    #[test]
    fn set_coords() {
        let mut part = Particle::new("O");
        assert_eq!(part.position, Vector3D::new(0.0, 0.0, 0.0));

        part.position(Vector3D::new(1.0, 2.0, 3.0));
        assert_eq!(part.position, Vector3D::new(1.0, 2.0, 3.0));

        part.velocity(Vector3D::new(1.0, 2.0, 3.0));
        assert_eq!(part.velocity, Vector3D::new(1.0, 2.0, 3.0));
    }
}
