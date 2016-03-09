// Cymbalum, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

//! `Particle` type and manipulation.
use types::Vector3D;
use super::PeriodicTable;

/// The Particle type hold basic data about a particle in the system. It is self
/// contained, so that it will be easy to send data between parrallels
/// processes.
#[derive(Clone, Debug)]
pub struct Particle {
    /// Particle name. This one is not public, as we always want to get &str,
    /// and to use either `String` of `&str` to set it.
    name: String,
    /// Particle kind, an index for potentials lookup
    pub kind: u16,
    /// Particle mass
    pub mass: f64,
    /// Particle charge
    pub charge: f64,
    /// Particle positions
    pub position: Vector3D,
    /// Particle velocity, if needed
    pub velocity: Vector3D,
}


impl Particle {
    /// Create a new `Particle` from a name
    pub fn new<S: Into<String>>(name: S) -> Particle {
        let name = name.into();
        let mass = PeriodicTable::mass(&name).unwrap_or_default();
        if mass == 0.0f32 {
            warn!("Could not find the mass for the {} particle", name);
        }
        Particle {
            name: name,
            mass: mass as f64,
            charge: 0.0,
            kind: u16::max_value(),
            position: Vector3D::new(0.0, 0.0, 0.0),
            velocity: Vector3D::new(0.0, 0.0, 0.0)
        }
    }

    /// Get the particle name
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Set the particle name to `name`
    pub fn set_name<'a, S>(&mut self, name: S) where S: Into<&'a str> {
        self.name = String::from(name.into());
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mass_initialization() {
        let part = Particle::new("O");
        assert_eq!(part.mass, 15.999f32 as f64);
    }

    #[test]
    fn name() {
        let mut part = Particle::new("O");
        assert_eq!(part.name(), "O");

        part.set_name("H");
        assert_eq!(part.name(), "H");

    }
}
