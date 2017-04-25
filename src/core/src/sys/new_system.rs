// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 Lumol's contributors — BSD license

//! The base type for simulation data in `lumol` is the `System` type.
//!
//! An `System` consists of a list of `Particle`; a list of `Molecule`
//! specifying how the particles are bonded together; an unit cell for boundary
//! conditions; and the interactions between these particles.
use std::ops::{Index, IndexMut};
use std::slice;
use std::cmp::{min, max};
use std::iter::IntoIterator;
use std::i8;
use std::collections::BTreeMap;

use energy::PairInteraction;
use energy::{BondPotential, AnglePotential, DihedralPotential};
use types::{Vector3D, Matrix3, Zero};

use super::{Particle, ParticleKind};
use super::Molecule;
use super::{CONNECT_12, CONNECT_13, CONNECT_14, CONNECT_FAR};
use super::UnitCell;
use super::interactions::Interactions;
use super::EnergyEvaluator;
use super::molecules::molecule_type;


/// A particle kind. Particles with the same name will have the same kind. This
/// is used for faster potential lookup.
// #[derive(Clone, Copy, Hash, PartialOrd, Ord, PartialEq, Eq, Debug)]
// pub struct ParticleKind(pub u32);

#[derive(Clone)]
pub struct Particles {
    names: Vec<String>,
    kinds: Vec<ParticleKind>,
    masses: Vec<f64>,
    charges: Vec<f64>,
    positions: Vec<Vector3D>,
    velocities: Vec<Vector3D>,
    molecule_ids: Vec<usize>
}

pub struct ParticlesRef<'a> {
    pub names: &'a [String],
    pub kinds: &'a [ParticleKind],
    pub masses: &'a [f64],
    pub charges: &'a [f64],
    pub positions: &'a [Vector3D],
    pub velocities: &'a [Vector3D],
    pub molecule_ids: &'a [usize]
}

pub struct ParticlesRefMut<'a> {
    pub names: &'a mut [String],
    pub kinds: &'a mut [ParticleKind],
    pub masses: &'a mut [f64],
    pub charges: &'a mut [f64],
    pub positions: &'a mut [Vector3D],
    pub velocities: &'a mut [Vector3D],
    pub molecule_ids: &'a mut [usize]
}


#[derive(Clone)]
pub struct SystemGeometry {
    particles: Particles,
    cell: UnitCell,
    molecules: Vec<Molecule>,
}

#[derive(Clone)]
pub struct System {
    geometry: SystemGeometry,
    interactions: Interactions,
    step: u64,
    external_temperature: Option<f64>
}

impl Default for Particles {
    fn default() -> Particles {
        Particles {
            names: Vec::new(),
            kinds: Vec::new(),
            masses: Vec::new(),
            charges: Vec::new(),
            positions: Vec::new(),
            velocities: Vec::new(),
            molecule_ids: Vec::new(),
        }
    }
}

impl Particles {
    #[inline] pub fn len(&self) -> usize {
        self.kinds.len()
    }

    pub fn as_slices(&self) -> ParticlesRef {
        ParticlesRef {
            names: &self.names,
            kinds: &self.kinds,
            masses: &self.masses,
            charges: &self.charges,
            positions: &self.positions,
            velocities: &self.velocities,
            molecule_ids: &self.molecule_ids
        }
    }

    pub fn as_slices_mut(&mut self) -> ParticlesRefMut {
        ParticlesRefMut {
            names: &mut self.names,
            kinds: &mut self.kinds,
            masses: &mut self.masses,
            charges: &mut self.charges,
            positions: &mut self.positions,
            velocities: &mut self.velocities,
            molecule_ids: &mut self.molecule_ids
        }
    }


    pub fn push(&mut self, particle: Particle, molecule_id: usize) {
        self.names.push(particle.name);
        self.kinds.push(particle.kind);
        self.masses.push(particle.mass);
        self.charges.push(particle.charge);
        self.positions.push(particle.position);
        self.velocities.push(particle.velocity);
        self.molecule_ids.push(molecule_id);
    }
}

impl<'a> ParticlesRef<'a> {
    #[inline] pub fn len(&self) -> usize {
        self.kinds.len()
    }
}

impl<'a> ParticlesRefMut<'a> {
    #[inline] pub fn len(&self) -> usize {
        self.kinds.len()
    }
}

impl Default for SystemGeometry {
    fn default() -> SystemGeometry {
        SystemGeometry {
            particles: Particles::default(),
            cell: UnitCell::new(),
            molecules: Vec::new()
        }
    }
}

impl SystemGeometry {
    pub fn particles(&self) -> ParticlesRef {
        self.particles.as_slices()
    }

    pub fn particles_mut(&mut self) -> ParticlesRefMut {
        self.particles.as_slices_mut()
    }
}


impl Default for System {
    fn default() -> System {
        System {
            geometry: SystemGeometry::default(),
            interactions: Interactions::new(),
            step: 0,
            external_temperature: None
        }
    }
}

impl System {
    pub fn from_cell(cell: UnitCell) -> System {
        let mut system = System::default();
        system.geometry.cell = cell;
        system
    }

    /// Get the current step of the system
    #[inline] pub fn step(&self) -> u64 {
        self.step
    }

    /// Increment the system step
    #[inline] pub fn increment_step(&mut self) {
        self.step += 1;
    }

    pub fn particles(&self) -> ParticlesRef {
        self.geometry.particles()
    }

    pub fn particles_mut(&mut self) -> ParticlesRefMut {
        self.geometry.particles_mut()
    }
}



/// Topology and particles related functions
impl System {
    /// Insert a particle at the end of the internal list
    pub fn add_particle(&mut self, p: Particle) {
        let mut part = p;
        if part.kind == ParticleKind::default() {
            // If no value have been precised, set one from the internal list
            // of particles kinds.
            part.kind = self.interactions.get_kind(part.name());
        }

        let i = self.particles().len();
        self.geometry.molecules.push(Molecule::new(i));
        let molecule_id = self.geometry.molecules.len() - 1;
        self.geometry.particles.push(part, molecule_id);
    }
}

#[cfg(test)]
mod tests {
    use super::System;
    use sys::Particle;

    #[test]
    fn step() {
        let mut system = System::default();

        assert_eq!(system.step(), 0);

        system.increment_step();
        system.increment_step();
        system.increment_step();

        assert_eq!(system.step(), 3);
    }

    #[test]
    fn particles() {
        let mut system = System::default();
        system.add_particle(Particle::new("O"));
        system.add_particle(Particle::new("H"));
        system.add_particle(Particle::new("H"));

        let particles = system.particles();
        assert_eq!(particles.len(), 3);
        assert_eq!(particles.names[0], "O");
        assert_eq!(particles.names[1], "H");
        assert_eq!(particles.names[2], "H");
    }
}
