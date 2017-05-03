// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 Lumol's contributors — BSD license

use std::ops::{Deref, DerefMut};
use std::slice;
use std::collections::BTreeMap;

use sys::{Configuration, Particle, ParticleKind, UnitCell};
use sys::Composition;

#[derive(Clone)]
pub struct System {
    configuration: Configuration,
    kinds: BTreeMap<String, ParticleKind>,
}

impl System {
    /// Create a new empty `System`
    pub fn new() -> System {
        System {
            configuration: Configuration::new(),
            kinds: BTreeMap::new(),
        }
    }

    /// Create an empty system with a specific unit cell
    pub fn with_cell(cell: UnitCell) -> System {
        let mut configuration = Configuration::new();
        configuration.cell = cell;
        System {
            configuration: configuration,
            kinds: BTreeMap::new(),
        }
    }

    fn get_kind(&mut self, name: &str) -> ParticleKind {
        match self.kinds.get(name).cloned() {
            Some(kind) => kind,
            None => {
                let kind = ParticleKind(self.kinds.len() as u32);
                let _ = self.kinds.insert(String::from(name), kind);
                kind
            }
        }
    }

    /// Insert a particle at the end of the internal list.
    pub fn add_particle(&mut self, mut particle: Particle) {
        if particle.kind == ParticleKind::invalid() {
            particle.kind = self.get_kind(particle.name());
        }
        self.configuration.add_particle(particle);
    }

    /// Get the number of particles of each kind in the configuration
    pub fn composition(&self) -> Composition {
        let mut composition = Composition::new();
        composition.resize(self.kinds.len());
        for particle in self {
            composition[particle.kind] += 1;
        }
        return composition;
    }

    /// Get a list of all the particles kinds in the system.
    pub fn particle_kinds(&self) -> Vec<ParticleKind> {
        self.kinds.values().cloned().collect()
    }
}

impl Deref for System {
    type Target = Configuration;

    fn deref(&self) -> &Configuration {
        &self.configuration
    }
}

impl DerefMut for System {
    fn deref_mut(&mut self) -> &mut Configuration {
        &mut self.configuration
    }
}

impl<'a> IntoIterator for &'a System {
    type Item = &'a Particle;
    type IntoIter = slice::Iter<'a, Particle>;

    #[inline]
    fn into_iter(self) -> slice::Iter<'a, Particle> {
        self.configuration.iter()
    }
}

impl<'a> IntoIterator for &'a mut System {
    type Item = &'a mut Particle;
    type IntoIter = slice::IterMut<'a, Particle>;

    #[inline]
    fn into_iter(self) -> slice::IterMut<'a, Particle> {
        self.configuration.iter_mut()
    }
}

#[cfg(test)]
mod tests {
    use super::System;
    use sys::*;

    #[test]
    fn deref() {
        let mut system = System::new();
        system.add_particle(Particle::new("H"));
        system.add_particle(Particle::new("O"));
        system.add_particle(Particle::new("H"));
        assert_eq!(system.molecules().len(), 3);

        // This uses deref_mut
        let _ = system.add_bond(0, 1);
        let _ = system.add_bond(2, 1);

        // This uses deref
        assert_eq!(system.molecules().len(), 1);
    }

    #[test]
    fn add_particle() {
        let mut system = System::new();
        system.add_particle(Particle::new("H"));
        system.add_particle(Particle::new("O"));
        system.add_particle(Particle::new("H"));

        assert_eq!(system[0].kind, ParticleKind(0));
        assert_eq!(system[1].kind, ParticleKind(1));
        assert_eq!(system[2].kind, ParticleKind(0));
    }

    #[test]
    fn particle_kinds() {
        let mut system = System::new();
        system.add_particle(Particle::new("H"));
        system.add_particle(Particle::new("O"));
        system.add_particle(Particle::new("O"));
        system.add_particle(Particle::new("H"));
        system.add_particle(Particle::new("C"));
        system.add_particle(Particle::new("U"));

        let kinds = system.particle_kinds();
        assert_eq!(kinds.len(), 4);
        assert!(kinds.contains(&ParticleKind(0)));
        assert!(kinds.contains(&ParticleKind(1)));
        assert!(kinds.contains(&ParticleKind(2)));
        assert!(kinds.contains(&ParticleKind(3)));
    }

    #[test]
    fn composition() {
        let mut system = System::new();
        system.add_particle(Particle::new("H"));
        system.add_particle(Particle::new("O"));
        system.add_particle(Particle::new("O"));
        system.add_particle(Particle::new("H"));
        system.add_particle(Particle::new("C"));
        system.add_particle(Particle::new("U"));
        system.add_particle(Particle::new("H"));

        let composition = system.composition();
        assert_eq!(composition.len(), 4);
        assert_eq!(composition[ParticleKind(0)], 3);
        assert_eq!(composition[ParticleKind(1)], 2);
        assert_eq!(composition[ParticleKind(2)], 1);
        assert_eq!(composition[ParticleKind(3)], 1);
    }
}
