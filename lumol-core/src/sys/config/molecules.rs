// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

use std::ops::Deref;
use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};

use soa_derive::soa_zip;

use crate::{Particle, ParticleVec, ParticleSlice, ParticleSliceMut};
use crate::{Bonding, UnitCell};
use crate::Vector3D;

/// A molecule hash allow to identify a molecule from its atoms and bonds, and
/// to know wether two molecules are the same without checking each atom and
/// bond.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct MoleculeHash(u64);

#[cfg(test)]
impl MoleculeHash {
    pub(crate) fn new(value: u64) -> MoleculeHash {
        MoleculeHash(value)
    }
}

/// A Molecule associate some particles bonded together.
///
/// [`Molecule`] implement `Deref` to a [`Bonding`] struct, to give read access
/// to all the bonds. It does not implement `DerefMut`, adding new bonds should
/// be done through [`Molecule::add_bond()`].
///
/// [`Molecule`]: struct.Molecule.html
/// [`Bonding`]: struct.Bonding.html
/// [`Molecule::add_bond()`]: struct.Molecule.html#method.add_bond
#[derive(Debug, Clone)]
pub struct Molecule {
    pub(crate) bonding: Bonding,
    pub(crate) particles: ParticleVec
}

/// An analog to [`&Molecule`] using particles stored elsewhere (in a system or
/// an [`Molecule`]).
///
/// [`MoleculeRef`] implement `Deref` to a [`Bonding`] struct, to give read
/// access to all the bonds.
///
/// [`Molecule`]: struct.Molecule.html
/// [`MoleculeRef`]: struct.MoleculeRef.html
#[derive(Debug)]
pub struct MoleculeRef<'a> {
    bonding: &'a Bonding,
    particles: ParticleSlice<'a>
}

/// An analog to [`&mut Molecule`] using particles stored elsewhere (in a
/// system or an [`Molecule`]).
///
/// [`MoleculeRefMut`] implement `Deref` to a [`Bonding`] struct, to give read
/// access to all the bonds. It does not implement `DerefMut`, adding new bonds
/// should be done through [`Molecule::add_bond()`] or
/// [`Configuration::add_bond()`].
///
/// [`Molecule`]: struct.Molecule.html
/// [`MoleculeRefMut`]: struct.MoleculeRefMut.html
/// [`Molecule::add_bond()`]: struct.Molecule.html#method.add_bond
/// [`Configuration::add_bond()`]: struct.Configuration.html#method.add_bond
#[derive(Debug)]
pub struct MoleculeRefMut<'a> {
    bonding: &'a Bonding,
    particles: ParticleSliceMut<'a>
}

impl Molecule {
    /// Create a new `Molecule` containing a single `particle`
    pub fn new(particle: Particle) -> Molecule {
        let mut particles = ParticleVec::new();
        particles.push(particle);
        Molecule {
            bonding: Bonding::new(0),
            particles: particles,
        }
    }

    /// Borrow `self` as a `MoleculeRef`.
    pub fn as_ref(&self) -> MoleculeRef<'_> {
        MoleculeRef {
            bonding: &self.bonding,
            particles: self.particles.as_slice(),
        }
    }

    /// Mutablely borrow `self` as a `MoleculeRefMut`.
    pub fn as_mut(&mut self) -> MoleculeRefMut<'_> {
        MoleculeRefMut {
            bonding: &self.bonding,
            particles: self.particles.as_mut_slice(),
        }
    }

    /// Get access to the particles in this molecule
    pub fn particles(&self) -> ParticleSlice<'_> {
        self.particles.as_slice()
    }

    /// Get mutable access to the particles in this molecule
    pub fn particles_mut(&mut self) -> ParticleSliceMut<'_> {
        self.particles.as_mut_slice()
    }

    /// Add a new `particle` in this molecule, bonded to an `other` particle
    /// in the molecule.
    pub fn add_particle_bonded_to(&mut self, other: usize, particle: Particle) {
        assert!(self.contains(other));
        self.particles.push(particle);
        let i = self.particles.len() - 1;
        self.bonding.merge_with(Bonding::new(i));
        self.bonding.add_bond(i, other);
    }

    /// Add bond between particles at indexes `i` and `j` in this molecule.
    ///
    /// # Panics
    ///
    /// If `i` or `j` are not in this molecule.
    pub fn add_bond(&mut self, i: usize, j: usize) {
        self.bonding.add_bond(i, j);
    }
}

impl Deref for Molecule {
    type Target = Bonding;

    fn deref(&self) -> &Self::Target {
        &self.bonding
    }
}

impl<'a> MoleculeRef<'a> {
    /// Create a new `MoleculeRef` associating the given `bonding` and
    /// `particles`.
    ///
    /// # Panics
    ///
    /// If the `bonding` and the `particles` do not containe the same number
    /// of particles.
    pub fn new(bonding: &'a Bonding, particles: ParticleSlice<'a>) -> MoleculeRef<'a> {
        assert_eq!(bonding.size(), particles.len());
        MoleculeRef {
            bonding: bonding,
            particles: particles,
        }
    }

    /// Get access to the particles in this molecule
    pub fn particles(&self) -> ParticleSlice<'_> {
        self.particles
    }

    /// Copies `self` into a new `Molecule`
    pub fn to_owned(&self) -> Molecule {
        Molecule {
            bonding: self.bonding.clone(),
            particles: self.particles.to_vec(),
        }
    }
}

impl<'a> Deref for MoleculeRef<'a> {
    type Target = Bonding;

    fn deref(&self) -> &Self::Target {
        self.bonding
    }
}

impl<'a> MoleculeRefMut<'a> {
    /// Create a new `MoleculeRefMut` associating the given `bonding` and
    /// `particles`.
    ///
    /// # Panics
    ///
    /// If the `bonding` and the `particles` do not containe the same number
    /// of particles.
    pub fn new(bonding: &'a Bonding, particles: ParticleSliceMut<'a>) -> MoleculeRefMut<'a> {
        assert_eq!(bonding.size(), particles.len());
        MoleculeRefMut {
            bonding: bonding,
            particles: particles,
        }
    }

    /// Borrow `self` as a `MoleculeRef`.
    pub fn as_ref(&self) -> MoleculeRef<'_> {
        MoleculeRef {
            bonding: self.bonding,
            particles: self.particles.as_ref(),
        }
    }

    /// Get access to the particles in this molecule
    pub fn particles(&self) -> ParticleSlice<'_> {
        self.particles.as_ref()
    }

    /// Get mutable access to the particles in this molecule
    pub fn particles_mut(&mut self) -> ParticleSliceMut<'_> {
        // Explicity re-borrow all the fiels, as ParticleSliceMut can not be
        // copied
        ParticleSliceMut {
            name: self.particles.name,
            mass: self.particles.mass,
            kind: self.particles.kind,
            charge: self.particles.charge,
            position: self.particles.position,
            velocity: self.particles.velocity,
        }
    }

    /// Copies `self` into a new `Molecule`
    pub fn to_owned(&self) -> Molecule {
        // This can not be a `ToOwned` implementation, as ToOwned requires
        // `Borrow`, and `Borrow` requires a reference, not a reference
        // wrapper.
        Molecule {
            bonding: self.bonding.clone(),
            particles: self.particles.to_vec(),
        }
    }
}

impl<'a> Deref for MoleculeRefMut<'a> {
    type Target = Bonding;

    fn deref(&self) -> &Self::Target {
        self.bonding
    }
}


// Add inherent functions in $body to all types in $Type
macro_rules! impl_on {
    ($($Type:ty,)+ => $body: tt) => (
        $(impl<'a> $Type $body)*
    );
}

impl_on!(Molecule, MoleculeRef<'a>, MoleculeRefMut<'a>, => {
    /// Return the center-of-mass of a molecule
    ///
    /// # Warning
    ///
    /// This function does not check for the particles' positions' nearest
    /// images. To use this function properly, make sure that all particles of
    /// the molecule are adjacent.
    pub fn center_of_mass(&self) -> Vector3D {
        let mut total_mass = 0.0;
        let mut com = Vector3D::zero();
        for (&mass, position) in soa_zip!(&self.particles, [mass, position]) {
            total_mass += mass;
            com += mass * position;
        }
        com / total_mass
    }

    /// Get a hash of this molecule. This is a hash of the particles names (in
    /// order), and the set of bonds in the molecule. This means that two
    /// molecules will have the same type if and only if they contains the same
    /// atoms and the same bonds, **in the same order**.
    pub fn hash(&self) -> MoleculeHash {
        let mut hasher = DefaultHasher::new();
        self.bonding.hash(&mut hasher);
        for name in self.particles().name {
            name.hash(&mut hasher);
        }
        MoleculeHash(hasher.finish())
    }
});

impl_on!(Molecule, MoleculeRefMut<'a>, => {
    /// Move all particles of a molecule such that the molecules center-of-mass
    /// position resides inside the simulation cell.
    ///
    /// # Note
    ///
    /// If the `CellShape` is `Infinite` there are no changes to the positions.
    pub fn wrap(&mut self, cell: &UnitCell) {
        let com = self.as_ref().center_of_mass();
        let mut com_wrapped = com;
        cell.wrap_vector(&mut com_wrapped);
        let delta = com_wrapped - com;
        // iterate over all positions and move them accordingly
        for position in self.particles_mut().position.iter_mut() {
            *position += delta;
        }
    }
});

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ParticleKind;

    use lazy_static::lazy_static;

    /// Create particles with intialized kind for the tests
    fn particle(name: &str) -> Particle {
        use std::collections::HashMap;
        use std::sync::Mutex;
        lazy_static! {
            static ref KINDS: Mutex<HashMap<String, u32>> = Mutex::new(HashMap::new());
        }
        let mut kind_map = KINDS.lock().unwrap();
        let nkinds = kind_map.len() as u32;
        let &mut kind = kind_map.entry(name.into()).or_insert(nkinds);

        let mut particle = Particle::new(name);
        particle.kind = ParticleKind(kind);
        return particle;
    }

    #[test]
    fn center_of_mass() {
        let mut molecule = Molecule::new(particle("O"));
        molecule.add_particle_bonded_to(0, particle("O"));

        molecule.particles_mut().position[0] = Vector3D::new(1.0, 0.0, 0.0);
        molecule.particles_mut().position[1] = Vector3D::zero();

        assert_eq!(molecule.center_of_mass(), Vector3D::new(0.5, 0.0, 0.0));
    }

    #[test]
    fn test_wrap_molecule() {
        let mut molecule = Molecule::new(particle("O"));
        molecule.add_particle_bonded_to(0, particle("O"));

        molecule.particles_mut().position[0] = Vector3D::new(-2.0, 0.0, 0.0);
        molecule.particles_mut().position[1] = Vector3D::zero();
        molecule.wrap(&UnitCell::cubic(5.0));

        assert_eq!(molecule.particles().position[0], Vector3D::new(3.0, 0.0, 0.0));
        assert_eq!(molecule.particles().position[1], Vector3D::new(5.0, 0.0, 0.0));
        assert_eq!(molecule.center_of_mass(), Vector3D::new(4.0, 0.0, 0.0));
    }
}
