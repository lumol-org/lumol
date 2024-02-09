// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

use std::collections::BTreeMap;
use std::collections::btree_map::Entry;

use crate::{ParticleKind, MoleculeHash};

/// The system composition contains the number of particles of each kind
/// in the system, as well as the number of molecules of each molecule type.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Composition {
    /// The particles composition, indexes by particles kind
    particles: Vec<usize>,
    /// The molecules compostion, indexes by molecule type
    molecules: BTreeMap<MoleculeHash, usize>,
}

impl Composition {
    /// Create a new empty composition
    ///
    /// # Examples
    ///
    /// ```
    /// # use lumol_core::sys::Composition;
    /// let composition = Composition::new();
    ///
    /// // no particles
    /// assert_eq!(composition.all_particles().count(), 0);
    ///
    /// // no molecules
    /// assert_eq!(composition.all_molecules().count(), 0);
    /// ```
    pub fn new() -> Composition {
        Composition {
            particles: Vec::new(),
            molecules: BTreeMap::new(),
        }
    }

    /// Add a particle with the given `kind` to the internal counter
    ///
    /// # Examples
    ///
    /// ```
    /// # use lumol_core::sys::{Composition, ParticleKind};
    /// let mut composition = Composition::new();
    ///
    /// composition.add_particle(ParticleKind(3));
    /// composition.add_particle(ParticleKind(3));
    ///
    /// assert_eq!(composition.particles(ParticleKind(3)), 2);
    /// ```
    pub fn add_particle(&mut self, kind: ParticleKind) {
        let i = kind.0 as usize;
        if i >= self.particles.len() {
            self.particles.resize(i + 1, 0);
        }
        self.particles[i] += 1;
    }

    /// Remove a particle with the given `kind` to the internal counter
    ///
    /// # Examples
    ///
    /// ```
    /// # use lumol_core::sys::{Composition, ParticleKind};
    /// let mut composition = Composition::new();
    ///
    /// composition.add_particle(ParticleKind(3));
    /// composition.add_particle(ParticleKind(3));
    ///
    /// assert_eq!(composition.particles(ParticleKind(3)), 2);
    ///
    /// composition.remove_particle(ParticleKind(3));
    ///
    /// assert_eq!(composition.particles(ParticleKind(3)), 1);
    /// ```
    pub fn remove_particle(&mut self, kind: ParticleKind) {
        let i = kind.0 as usize;
        assert!(i < self.particles.len());
        assert!(self.particles[i] > 0);
        self.particles[i] -= 1;
    }

    /// Get the number of particles with a given kind
    ///
    /// # Examples
    ///
    /// ```
    /// # use lumol_core::sys::{Composition, ParticleKind};
    /// let mut composition = Composition::new();
    ///
    /// assert_eq!(composition.particles(ParticleKind(3)), 0);
    ///
    /// composition.add_particle(ParticleKind(3));
    /// composition.add_particle(ParticleKind(3));
    ///
    /// assert_eq!(composition.particles(ParticleKind(3)), 2);
    /// ```
    pub fn particles(&self, kind: ParticleKind) -> usize {
        let i = kind.0 as usize;
        return self.particles.get(i).copied().unwrap_or(0);
    }

    /// Get an iterator over the particles kind and count
    ///
    /// # Examples
    ///
    /// ```
    /// # use lumol_core::sys::{Composition, ParticleKind};
    /// let mut composition = Composition::new();
    ///
    /// composition.add_particle(ParticleKind(10));
    /// composition.add_particle(ParticleKind(4));
    /// composition.add_particle(ParticleKind(4));
    ///
    /// let mut iter = composition.all_particles();
    /// assert_eq!(iter.next(), Some((ParticleKind(4), 2)));
    /// assert_eq!(iter.next(), Some((ParticleKind(10), 1)));
    /// assert_eq!(iter.next(), None);
    /// ```
    pub fn all_particles(&self) -> impl Iterator<Item = (ParticleKind, usize)> + '_ {
        self.particles
            .iter()
            .enumerate()
            .filter(|&(_, &n)| n != 0)
            .map(|(i, &n)| (ParticleKind(i as u32), n))
    }

    /// Add a molecule with the given `moltype` to the internal counter
    ///
    /// # Examples
    ///
    /// ```
    /// # use lumol_core::sys::{Composition, Molecule, Particle};
    /// // Getting a molecule hash
    /// let he = Molecule::new(Particle::new("He")).hash();
    ///
    /// let mut composition = Composition::new();
    /// assert_eq!(composition.molecules(he), 0);
    ///
    /// composition.add_molecule(he);
    /// composition.add_molecule(he);
    ///
    /// assert_eq!(composition.molecules(he), 2);
    /// ```
    pub fn add_molecule(&mut self, hash: MoleculeHash) {
        *self.molecules.entry(hash).or_insert(0) += 1;
    }

    /// Add a molecule with the given `moltype` to the internal counter
    ///
    /// # Examples
    ///
    /// ```
    /// # use lumol_core::sys::{Composition, Molecule, Particle};
    /// // Getting a molecule hash
    /// let he = Molecule::new(Particle::new("He")).hash();
    ///
    /// let mut composition = Composition::new();
    /// assert_eq!(composition.molecules(he), 0);
    ///
    /// composition.add_molecule(he);
    /// composition.add_molecule(he);
    ///
    /// assert_eq!(composition.molecules(he), 2);
    ///
    /// composition.remove_molecule(he);
    ///
    /// assert_eq!(composition.molecules(he), 1);
    /// ```
    pub fn remove_molecule(&mut self, hash: MoleculeHash) {
        if let Entry::Occupied(mut count) = self.molecules.entry(hash) {
            let count = count.get_mut();
            if *count > 0 {
                *count -= 1;
            }
        }
    }

    /// Get the number of particles with the given `hash`
    ///
    /// # Examples
    ///
    /// ```
    /// # use lumol_core::sys::{Composition, Molecule, Particle};
    /// // Getting some molecules hashes
    /// let he = Molecule::new(Particle::new("He")).hash();
    ///
    /// let ar = Molecule::new(Particle::new("Ar")).hash();
    ///
    /// let mut n2 = Molecule::new(Particle::new("N"));
    /// n2.add_particle_bonded_to(0, Particle::new("N"));
    /// let n2 = n2.hash();
    ///
    /// let mut composition = Composition::new();
    ///
    /// composition.add_molecule(he);
    /// composition.add_molecule(he);
    /// composition.add_molecule(n2);
    ///
    /// assert_eq!(composition.molecules(he), 2);
    /// assert_eq!(composition.molecules(n2), 1);
    /// assert_eq!(composition.molecules(ar), 0);
    /// ```
    pub fn molecules(&self, hash: MoleculeHash) -> usize {
        self.molecules.get(&hash).copied().unwrap_or(0)
    }

    /// Get an iterator over the molecules hashes and count
    ///
    /// # Examples
    ///
    /// ```
    /// # use lumol_core::sys::{Composition, Particle, Molecule};
    /// let mut composition = Composition::new();
    /// let he = Molecule::new(Particle::new("He")).hash();
    /// let ar = Molecule::new(Particle::new("Ar")).hash();
    ///
    /// composition.add_molecule(he);
    /// composition.add_molecule(he);
    /// composition.add_molecule(ar);
    ///
    /// let mut iter = composition.all_molecules();
    /// assert_eq!(iter.next(), Some((he, 2)));
    /// assert_eq!(iter.next(), Some((ar, 1)));
    /// assert_eq!(iter.next(), None);
    /// ```
    pub fn all_molecules(&self) -> impl Iterator<Item = (MoleculeHash, usize)> + '_ {
        self.molecules
            .iter()
            .filter(|(_, n)| **n != 0)
            .map(|(&m, &n)| (m, n))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn particles() {
        let mut composition = Composition::new();
        for _ in 0..2 {
            composition.add_particle(ParticleKind(0));
            for _ in 0..2 {
                composition.add_particle(ParticleKind(1));
                for _ in 0..2 {
                    composition.add_particle(ParticleKind(2));
                }
            }
        }

        let mut expected = 2;
        for (_, count) in composition.all_particles() {
            assert_eq!(count, expected);
            expected *= 2;
        }
    }

    #[test]
    fn molecules() {
        let mut composition = Composition::new();

        composition.add_molecule(MoleculeHash::new(22));
        composition.add_molecule(MoleculeHash::new(10));
        composition.add_molecule(MoleculeHash::new(22));

        assert_eq!(composition.molecules(MoleculeHash::new(22)), 2);
        assert_eq!(composition.molecules(MoleculeHash::new(10)), 1);
        assert_eq!(composition.molecules(MoleculeHash::new(124)), 0);
    }
}
