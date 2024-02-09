// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 Lumol's contributors — BSD license

use std::cmp::{max, min};
use std::marker::PhantomData;

use log::trace;
use log_once::warn_once;

use crate::Vector3D;
use crate::{BondDistances, Bonding, ParticleKind, UnitCell};
use crate::{ParticleSlice, ParticleSliceMut, ParticleVec, ParticlePtr, ParticlePtrMut};
use crate::{Molecule, MoleculeRef, MoleculeRefMut};
use crate::BondPath;

/// The `Permutation` struct contains the old and new particle index in a
/// `Configuration` after the particles where moved due to a new bond being
/// added
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct Permutation {
    /// The old particle index
    pub old: usize,
    /// The new particle index
    pub new: usize,
}

impl Permutation {
    fn new(old: usize, new: usize) -> Permutation {
        Permutation {
            old: old,
            new: new,
        }
    }
}

/// The `Configuration` contains the physical data of the system:
///
/// - The unit cell;
/// - The list of particles in the system;
/// - The list of molecules in the system.
#[derive(Clone)]
pub struct Configuration {
    /// Unit cell of the system
    pub cell: UnitCell,
    /// List of particles in the system
    particles: ParticleVec,
    /// Bonding information in the system
    bondings: Vec<Bonding>,
    /// Molecules indexes for all the particles
    molecule_ids: Vec<usize>,
}

impl Configuration {
    /// Create a new empty `Configuration`
    pub fn new() -> Configuration {
        Configuration {
            particles: ParticleVec::new(),
            bondings: Vec::new(),
            molecule_ids: Vec::new(),
            cell: UnitCell::infinite(),
        }
    }
}

/// Topology and particles related functions
impl Configuration {
    /// Check if the particles at indexes `i` and `j` are in the same molecule
    pub fn are_in_same_molecule(&self, i: usize, j: usize) -> bool {
        debug_assert_eq!(self.molecule_ids.len(), self.particles.len());
        self.molecule_ids[i] == self.molecule_ids[j]
    }

    /// Get an iterator over the molecules in the configuration.
    pub fn molecules(&self) -> MoleculeIter<'_> {
        let ptr = self.particles.as_ptr();
        let end = unsafe {
            ptr.add(self.particles.len())
        };
        MoleculeIter {
            bondings: self.bondings.iter(),
            ptr: ptr,
            end: end,
            _marker: PhantomData,
        }
    }

    /// Get an iterator over the molecules in the configuration.
    pub fn molecules_mut(&mut self) -> MoleculeIterMut<'_> {
        let ptr = self.particles.as_mut_ptr();
        let end = unsafe {
            ptr.add(self.particles.len())
        };
        MoleculeIterMut {
            bondings: self.bondings.iter(),
            ptr: ptr,
            end: end,
            _marker: PhantomData,
        }
    }

    /// Get the molecule at index `id`
    pub fn molecule(&self, id: usize) -> MoleculeRef<'_> {
        let bonding = &self.bondings[id];
        let particles = self.particles.slice(bonding.indexes());
        MoleculeRef::new(bonding, particles)
    }

    /// Get the molecule at index `id`
    pub fn molecule_mut(&mut self, id: usize) -> MoleculeRefMut<'_> {
        let bonding = &mut self.bondings[id];
        let particles = self.particles.slice_mut(bonding.indexes());
        MoleculeRefMut::new(bonding, particles)
    }

    /// Get the index of the molecule containing the particle `i`
    pub fn molecule_id(&self, i: usize) -> usize {
        self.molecule_ids[i]
    }

    /// Get the length of the shortest bond path to go from the particle `i` to
    /// the particle `j`. If the particles are not in the same molecule, the
    /// length is -1. Else, this length is 0 if `i == j`, 1 if there is a bond
    /// between `i` and `j`, etc.
    pub fn bond_path(&self, i: usize, j: usize) -> BondPath {
        assert!(i < self.size() && j < self.size());
        if !(self.are_in_same_molecule(i, j)) {
            BondPath::None
        } else if i == j {
            BondPath::SameParticle
        } else {
            let connect = self.molecule(self.molecule_id(i)).bond_distances(i, j);
            if connect.contains(BondDistances::ONE) {
                BondPath::OneBond
            } else if connect.contains(BondDistances::TWO) {
                BondPath::TwoBonds
            } else if connect.contains(BondDistances::THREE) {
                BondPath::ThreeBonds
            } else if connect.contains(BondDistances::FAR) {
                BondPath::Far
            } else {
                unreachable!();
            }
        }
    }

    /// Remove the molecule at index `i`
    pub fn remove_molecule(&mut self, molid: usize) {
        let molecule = self.bondings.remove(molid);
        let first = molecule.start();
        let size = molecule.size();

        for _ in 0..size {
            let _ = self.particles.remove(first);
            let _ = self.molecule_ids.remove(first);
        }

        for molecule in self.bondings.iter_mut().skip(molid) {
            molecule.translate_by(-(size as isize));
        }

        for molid in self.molecule_ids.iter_mut().skip(first) {
            *molid -= 1;
        }
    }

    /// Add a bond between the particles at indexes `i` and `j`. The particles
    /// should have been added to the configuration before calling this.
    ///
    /// # Warning
    ///
    /// If the bond is between two particles which are not in the same molecule,
    /// the two molecules are merged together by moving particles in the
    /// particles list, and thus invalidate any previously stored index. In
    /// particular, any bond, angle, dihedral or molecule is invalidated.
    ///
    /// This function will return the list of atomic permutations that where
    /// applied in order to ensure that molecules are contiguous in memory.
    pub fn add_bond(&mut self, mut particle_i: usize, mut particle_j: usize) -> Vec<Permutation> {
        assert!(particle_i <= self.particles.len());
        assert!(particle_j <= self.particles.len());
        assert_ne!(particle_i, particle_j);
        trace!(
            "Adding bond {}-{} between molecules {} and {}",
            particle_i,
            particle_j,
            self.molecule_ids[particle_i],
            self.molecule_ids[particle_j]
        );

        // Getting copy of the molecules before the merge
        let molid_i = self.molecule_ids[particle_i];
        let molid_j = self.molecule_ids[particle_j];
        let new_molid = min(molid_i, molid_j);
        let old_molid = max(molid_i, molid_j);
        let new_mol = self.bondings[new_molid].clone();
        let old_mol = self.bondings[old_molid].clone();
        let already_in_same_molecule = self.are_in_same_molecule(particle_i, particle_j);

        // Effective merge
        let delta = self.merge_molecules(molid_i, molid_j);

        let mut permutations = Vec::new();
        // If new_mol.last() + 1 == old_mol.first(), no one moved. Else,
        // we generate the permutations
        if !already_in_same_molecule && new_mol.end() != old_mol.start() {
            let size = old_mol.size();
            let first = old_mol.start();
            let second = new_mol.end();
            // Add permutation for the molecule we just moved around
            for i in 0..size {
                permutations.push(Permutation::new(first + i, second + i));
            }

            // Add permutations for molecules that where shifted to make
            // space for the just moved molecule.
            for connectivity in &self.bondings[new_molid + 1..old_molid] {
                for i in connectivity.indexes() {
                    permutations.push(Permutation::new(i - size, i));
                }
            }
        }

        // One of the `particle_i` or `particle_j` index is no longer valid, as
        // one molecule has been displaced.
        if molid_i == new_molid {
            particle_j -= delta; // j moved
        } else {
            particle_i -= delta; // i moved
        };

        assert_eq!(self.molecule_ids[particle_i], self.molecule_ids[particle_j]);
        self.bondings[self.molecule_ids[particle_i]].add_bond(particle_i, particle_j);
        return permutations;
    }

    /// Add a molecule to the configuration, putting the new particles at the
    /// end of the particles list
    pub fn add_molecule(&mut self, mut molecule: Molecule) {
        for particle in molecule.particles() {
            assert_ne!(*particle.kind, ParticleKind::invalid());
            if *particle.mass < 0.0 || f64::is_nan(*particle.mass) {
                warn_once!(
                    "Adding a particle ({}) with an invalid mass: {}",
                    particle.name, particle.mass
                );
            }
        }

        let mut bonding = molecule.bonding;
        bonding.translate_by(self.particles.len() as isize);

        self.molecule_ids.append(&mut vec![self.bondings.len(); bonding.size()]);
        self.bondings.push(bonding);
        self.particles.append(&mut molecule.particles);
    }

    /// Get the number of particles in this configuration
    #[inline]
    pub fn size(&self) -> usize {
        self.particles.len()
    }

    /// Check if this configuration contains any particle
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.particles.is_empty()
    }

    /// Return the center-of-mass of the configuration
    pub fn center_of_mass(&self) -> Vector3D {
        let mut total_mass = 0.0;
        let mut com = Vector3D::zero();
        for i in 0..self.size() {
            total_mass += self.particles.mass[i];
            com += self.particles.mass[i] * self.particles.position[i];
        }
        com / total_mass
    }

    /// Get the list of particles in this configuration, as a `ParticleSlice`.
    pub fn particles(&self) -> ParticleSlice<'_> {
        self.particles.as_slice()
    }

    /// Get the list of particles in this configuration, as a mutable
    /// `ParticleSliceMut`.
    pub fn particles_mut(&mut self) -> ParticleSliceMut<'_> {
        self.particles.as_mut_slice()
    }

    /// Merge the molecules at indexes `first` and `second` into one
    /// molecule. The molecule are merged into the one with the lower molecule
    /// index.
    ///
    /// For example, if we have
    /// ```text
    ///  0    1    2    3   # Molecules indexes
    /// H-H  H-H  H-H  H-H
    /// 0 1  2 3  4 5  6 7  # Particles indexes
    /// ```
    /// and call `merge_molecules(0, 3)` the result will be
    /// ```text
    /// 0 1 6 7  2 3  4 5  # Old indexes
    /// H-H-H-H  H-H  H-H
    /// 0 1 2 3  4 5  6 7  # New indexes
    /// ```
    ///
    /// This functions return the change in index of the first particle of the
    /// moved molecule, i.e. in this example `4`.
    fn merge_molecules(&mut self, first: usize, second: usize) -> usize {
        // Do not try to merge a molecule with itself
        if first == second {
            return 0;
        }

        let (new_molid, old_molid) = if first < second {
            (first, second)
        } else {
            (second, first)
        };

        let mut new_mol = self.bondings[new_molid].clone();
        let old_mol = self.bondings[old_molid].clone();

        if new_mol.end() == old_mol.start() {
            // Just update the molecules ids
            for i in old_mol.indexes() {
                self.molecule_ids[i] = new_molid;
            }
        } else {
            // Move the particles close together
            let mut new_index = new_mol.end();
            for i in old_mol.indexes() {
                // Remove particles from the old position, and insert it to the
                // new one. The indexes are valid during the movement, because
                // we insert a new particle for each particle removed.
                let particle = self.particles.remove(i);
                self.particles.insert(new_index, particle);

                // Update molecule_ids
                let _ = self.molecule_ids.remove(i);
                self.molecule_ids.insert(new_index, new_molid);

                new_index += 1;
            }
        }

        let mut old_mol = self.bondings[old_molid].clone();
        let size = old_mol.size() as isize;

        // translate all indexes in the molecules between new_mol and old_mol
        for molecule in
            self.bondings.iter_mut().skip(new_molid + 1).take(old_molid - new_molid - 1)
        {
            molecule.translate_by(size);
        }

        // Update molid for all particles after the old molecule
        for molid in self.molecule_ids.iter_mut().skip(old_mol.end()) {
            *molid -= 1;
        }

        let delta = old_mol.start() - new_mol.end();
        old_mol.translate_by(-(delta as isize));

        new_mol.merge_with(old_mol);
        self.bondings[new_molid] = new_mol;
        let _ = self.bondings.remove(old_molid);

        debug_assert!(check_molid_sorted(&self.molecule_ids), "Unsorted molecule ids {:?}", self.molecule_ids);

        return delta;
    }
}

/// Check that `molecule_ids` is sorted and only contains successive values
fn check_molid_sorted(molecule_ids: &[usize]) -> bool {
    let mut previous = 0;
    for &i in molecule_ids {
        if i == previous || i == previous + 1 {
            previous = i;
        } else {
            return false;
        }
    }
    return true;
}

/// `UnitCell` related functions
impl Configuration {
    /// Get the distance between the particles at indexes `i` and `j`
    #[inline]
    pub fn distance(&self, i: usize, j: usize) -> f64 {
        self.cell.distance(&self.particles.position[i], &self.particles.position[j])
    }

    /// Get the vector between the nearest image of particle `j` with respect to
    /// particle `i`.
    pub fn nearest_image(&self, i: usize, j: usize) -> Vector3D {
        let mut res = self.particles.position[i] - self.particles.position[j];
        self.cell.vector_image(&mut res);
        return res;
    }

    /// Get the angle between the particles `i`, `j` and `k`
    pub fn angle(&self, i: usize, j: usize, k: usize) -> f64 {
        self.cell.angle(
            &self.particles.position[i],
            &self.particles.position[j],
            &self.particles.position[k]
        )
    }

    /// Get the angle and the derivatives of the angle between the particles
    /// `i`, `j` and `k`
    pub fn angle_and_derivatives(
        &self,
        i: usize,
        j: usize,
        k: usize,
    ) -> (f64, Vector3D, Vector3D, Vector3D) {
      self.cell.angle_and_derivatives(
          &self.particles.position[i],
          &self.particles.position[j],
          &self.particles.position[k]
      )
    }

    /// Get the dihedral angle between the particles `i`, `j`, `k` and `m`
    pub fn dihedral(&self, i: usize, j: usize, k: usize, m: usize) -> f64 {
        self.cell.dihedral(
            &self.particles.position[i],
            &self.particles.position[j],
            &self.particles.position[k],
            &self.particles.position[m]
        )
    }

    /// Get the dihedral angle and the derivatives of the dihedral angle
    /// between the particles `i`, `j`, `k` and `m`
    pub fn dihedral_and_derivatives(
        &self,
        i: usize,
        j: usize,
        k: usize,
        m: usize,
    ) -> (f64, Vector3D, Vector3D, Vector3D, Vector3D) {
        self.cell.dihedral_and_derivatives(
            &self.particles.position[i],
            &self.particles.position[j],
            &self.particles.position[k],
            &self.particles.position[m]
        )
    }
}

/// An iterator over all the molecules in a `Configuration`
pub struct MoleculeIter<'a> {
    bondings: ::std::slice::Iter<'a, Bonding>,
    ptr: ParticlePtr,
    end: ParticlePtr,
    _marker: PhantomData<ParticleSlice<'a>>
}

unsafe impl<'a> Send for MoleculeIter<'a> {}

impl<'a> Iterator for MoleculeIter<'a> {
    type Item = MoleculeRef<'a>;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        if let Some(bonding) = self.bondings.next() {
            debug_assert!(self.ptr.name < self.end.name);
            let len = bonding.size();
            let slice = unsafe {
                let slice = ParticleSlice::from_raw_parts(self.ptr, len);
                self.ptr = self.ptr.add(len);
                slice
            };
            debug_assert!(self.ptr.name <= self.end.name);

            Some(MoleculeRef::new(bonding, slice))
        } else {
            None
        }
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        self.bondings.size_hint()
    }

    #[inline]
    fn count(self) -> usize {
        self.bondings.count()
    }
}

impl<'a> DoubleEndedIterator for MoleculeIter<'a> {
    #[inline]
    fn next_back(&mut self) -> Option<Self::Item> {
        if let Some(bonding) = self.bondings.next_back() {
            debug_assert!(self.ptr.name < self.end.name);
            let len = bonding.size();
            let slice = unsafe {
                self.end = self.end.sub(len);
                ParticleSlice::from_raw_parts(self.end, len)
            };
            debug_assert!(self.ptr.name <= self.end.name);

            Some(MoleculeRef::new(bonding, slice))
        } else {
            None
        }
    }
}

/// A mutable iterator over all the molecules in a `Configuration`
pub struct MoleculeIterMut<'a> {
    bondings: ::std::slice::Iter<'a, Bonding>,
    ptr: ParticlePtrMut,
    end: ParticlePtrMut,
    _marker: PhantomData<ParticleSliceMut<'a>>
}

unsafe impl<'a> Send for MoleculeIterMut<'a> {}

impl<'a> Iterator for MoleculeIterMut<'a> {
    type Item = MoleculeRefMut<'a>;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        if let Some(bonding) = self.bondings.next() {
            debug_assert!(self.ptr.name < self.end.name);
            let len = bonding.size();
            let slice = unsafe {
                let slice = ParticleSliceMut::from_raw_parts_mut(self.ptr, len);
                self.ptr = self.ptr.add(len);
                slice
            };
            debug_assert!(self.ptr.name <= self.end.name);

            Some(MoleculeRefMut::new(bonding, slice))
        } else {
            None
        }
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        self.bondings.size_hint()
    }

    #[inline]
    fn count(self) -> usize {
        self.bondings.count()
    }
}

impl<'a> DoubleEndedIterator for MoleculeIterMut<'a> {
    #[inline]
    fn next_back(&mut self) -> Option<Self::Item> {
        if let Some(bonding) = self.bondings.next_back() {
            debug_assert!(self.ptr.name < self.end.name);
            let len = bonding.size();
            let slice = unsafe {
                self.end = self.end.sub(len);
                ParticleSliceMut::from_raw_parts_mut(self.end, len)
            };
            debug_assert!(self.ptr.name <= self.end.name);

            Some(MoleculeRefMut::new(bonding, slice))
        } else {
            None
        }
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::Vector3D;
    use crate::{Angle, Bond, Dihedral, Particle, Molecule};
    use crate::BondPath;

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
    fn molecules() {
        let mut configuration = Configuration::new();

        let mut molecule = Molecule::new(particle("H"));
        molecule.add_particle_bonded_to(0, particle("O"));
        configuration.add_molecule(molecule);

        let mut molecule = Molecule::new(particle("0"));
        molecule.add_particle_bonded_to(0, particle("H"));
        configuration.add_molecule(molecule);

        assert_eq!(configuration.molecules().count(), 2);

        assert!(configuration.molecule(0).bonds().contains(&Bond::new(0, 1)));
        assert!(configuration.molecule(1).bonds().contains(&Bond::new(2, 3)));

        assert!(configuration.add_bond(1, 2).is_empty());
        assert_eq!(configuration.molecules().count(), 1);

        let molecule = configuration.molecule(0).clone();
        assert!(molecule.angles().contains(&Angle::new(0, 1, 2)));
        assert!(molecule.angles().contains(&Angle::new(1, 2, 3)));
        assert!(molecule.dihedrals().contains(&Dihedral::new(0, 1, 2, 3)));

        let molid = configuration.molecule_id(1);
        configuration.remove_molecule(molid);
        assert_eq!(configuration.molecules().count(), 0);
        assert_eq!(configuration.size(), 0);
    }

    mod iterators {
        use super::super::Configuration;
        use super::particle;
        use crate::Molecule;

        #[test]
        fn count() {
            let mut configuration = Configuration::new();
            assert_eq!(configuration.molecules().count(), 0);
            assert_eq!(configuration.molecules_mut().count(), 0);

            configuration.add_molecule(Molecule::new(particle("Ar")));
            configuration.add_molecule(Molecule::new(particle("Ar")));
            configuration.add_molecule(Molecule::new(particle("He")));

            assert_eq!(configuration.molecules().count(), 3);
            assert_eq!(configuration.molecules_mut().count(), 3);
        }

        #[test]
        fn next() {
            let mut configuration = Configuration::new();

            {
                let mut iter = configuration.molecules();
                assert!(iter.next().is_none());
                assert_eq!(iter.ptr.name, iter.end.name);
            }
            {
                let mut iter = configuration.molecules_mut();
                assert!(iter.next().is_none());
                assert_eq!(iter.ptr.name, iter.end.name);
            }

            configuration.add_molecule(Molecule::new(particle("Ar")));

            let mut molecule = Molecule::new(particle("H"));
            molecule.add_particle_bonded_to(0, particle("O"));
            molecule.add_particle_bonded_to(1, particle("H"));
            configuration.add_molecule(molecule);

            for (i, molecule) in configuration.molecules().enumerate() {
                if i == 0 {
                    assert_eq!(molecule.size(), 1);
                    assert_eq!(molecule.particles().name[0], "Ar");
                } else if i == 1 {
                    assert_eq!(molecule.size(), 3);
                    assert_eq!(molecule.particles().name[0], "H");
                    assert_eq!(molecule.particles().name[1], "O");
                    assert_eq!(molecule.particles().name[2], "H");
                }
            }

            for mut molecule in configuration.molecules_mut() {
                for name in molecule.particles_mut().name {
                    if *name == "Ar" {
                        *name = String::from("He");
                    }
                }
            }

            assert_eq!(configuration.molecule(0).particles().name[0], "He");

            {
                let mut iter = configuration.molecules();
                assert_eq!(iter.next().unwrap().size(), 1);
                assert_eq!(iter.next().unwrap().size(), 3);
                assert!(iter.next().is_none());
            }

            {
                let mut iter = configuration.molecules_mut();
                assert_eq!(iter.next().unwrap().size(), 1);
                assert_eq!(iter.next().unwrap().size(), 3);
                assert!(iter.next().is_none());
            }
        }

        #[test]
        fn next_back() {
            let mut configuration = Configuration::new();

            {
                let mut iter = configuration.molecules();
                assert_eq!(iter.ptr.name, iter.end.name);
                assert!(iter.next_back().is_none());
            }
            {
                let mut iter = configuration.molecules_mut();
                assert_eq!(iter.ptr.name, iter.end.name);
                assert!(iter.next_back().is_none());
            }

            configuration.add_molecule(Molecule::new(particle("Ar")));

            let mut molecule = Molecule::new(particle("H"));
            molecule.add_particle_bonded_to(0, particle("O"));
            molecule.add_particle_bonded_to(1, particle("H"));
            configuration.add_molecule(molecule);

            {
                let mut iter = configuration.molecules();
                assert_eq!(iter.next_back().unwrap().size(), 3);
                assert_eq!(iter.next().unwrap().size(), 1);
                assert!(iter.next().is_none());
                assert!(iter.next_back().is_none());
            }

            {
                let mut iter = configuration.molecules_mut();
                assert_eq!(iter.next_back().unwrap().size(), 3);
                assert_eq!(iter.next().unwrap().size(), 1);
                assert!(iter.next().is_none());
                assert!(iter.next_back().is_none());
            }
        }
    }

    #[test]
    fn bond_path() {
        let mut configuration = Configuration::new();

        let mut pentane = Molecule::new(particle("CH3"));
        pentane.add_particle_bonded_to(0, particle("CH2"));
        pentane.add_particle_bonded_to(1, particle("CH2"));
        pentane.add_particle_bonded_to(2, particle("CH2"));
        pentane.add_particle_bonded_to(3, particle("CH3"));

        configuration.add_molecule(pentane);
        configuration.add_molecule(Molecule::new(particle("Zn")));

        assert_eq!(configuration.bond_path(0, 0), BondPath::SameParticle);
        assert_eq!(configuration.bond_path(0, 1), BondPath::OneBond);
        assert_eq!(configuration.bond_path(0, 2), BondPath::TwoBonds);
        assert_eq!(configuration.bond_path(0, 3), BondPath::ThreeBonds);
        assert_eq!(configuration.bond_path(0, 4), BondPath::Far);
        assert_eq!(configuration.bond_path(0, 5), BondPath::None);
    }

    #[test]
    fn add_bond_permutations() {
        let mut configuration = Configuration::new();
        configuration.add_molecule(Molecule::new(particle("C")));
        configuration.add_molecule(Molecule::new(particle("H")));
        configuration.add_molecule(Molecule::new(particle("H")));
        configuration.add_molecule(Molecule::new(particle("H")));

        configuration.add_molecule(Molecule::new(particle("C")));
        configuration.add_molecule(Molecule::new(particle("H")));
        configuration.add_molecule(Molecule::new(particle("H")));
        configuration.add_molecule(Molecule::new(particle("H")));

        assert_eq!(configuration.add_bond(0, 3), vec![Permutation::new(3, 1), Permutation::new(1, 2), Permutation::new(2, 3)]);
        assert_eq!(configuration.add_bond(0, 3), vec![Permutation::new(3, 2), Permutation::new(2, 3)]);
        assert_eq!(configuration.add_bond(0, 3), vec![]);

        assert_eq!(configuration.add_bond(4, 5), vec![]);
        assert_eq!(configuration.add_bond(4, 7), vec![Permutation::new(7, 6), Permutation::new(6, 7)]);
        assert_eq!(configuration.add_bond(4, 7), vec![]);

        // This is a regression test for issue #76
        let mut configuration = Configuration::new();
        configuration.add_molecule(Molecule::new(particle("H")));
        configuration.add_molecule(Molecule::new(particle("H")));
        configuration.add_molecule(Molecule::new(particle("O")));

        assert_eq!(configuration.add_bond(0, 2), vec![Permutation::new(2, 1), Permutation::new(1, 2)]);
        assert_eq!(configuration.add_bond(2, 1), vec![]);
        assert_eq!(configuration.molecules().count(), 1);
    }

    #[test]
    fn particles() {
        let mut configuration = Configuration::new();
        configuration.add_molecule(Molecule::new(particle("O")));
        configuration.add_molecule(Molecule::new(particle("H")));
        configuration.add_molecule(Molecule::new(particle("H")));

        assert_eq!(configuration.size(), 3);
        assert_eq!(configuration.particles().name[0], "O");
        assert_eq!(configuration.particles().name[1], "H");
        assert_eq!(configuration.particles().name[2], "H");
    }

    #[test]
    fn distances() {
        let mut configuration = Configuration::new();
        configuration.cell = UnitCell::cubic(5.0);
        configuration.add_molecule(Molecule::new(particle("O")));
        configuration.add_molecule(Molecule::new(particle("H")));

        configuration.particles_mut().position[0] = Vector3D::new(9.0, 0.0, 0.0);
        configuration.particles_mut().position[1] = Vector3D::zero();
        assert_eq!(configuration.distance(0, 1), 1.0);

        configuration.cell = UnitCell::infinite();
        assert_eq!(configuration.distance(0, 1), 9.0);
    }

    #[test]
    fn hash() {
        let mut configuration = Configuration::new();
        // One helium
        configuration.add_molecule(Molecule::new(particle("He")));

        // Two water molecules
        let mut water = Molecule::new(particle("H"));
        water.add_particle_bonded_to(0, particle("O"));
        water.add_particle_bonded_to(1, particle("H"));
        configuration.add_molecule(water.clone());
        configuration.add_molecule(water);

        // Another helium
        configuration.add_molecule(Molecule::new(particle("He")));

        // A water molecules, with different atoms order
        let mut water = Molecule::new(particle("O"));
        water.add_particle_bonded_to(0, particle("H"));
        water.add_particle_bonded_to(0, particle("H"));
        configuration.add_molecule(water);

        assert_eq!(configuration.molecules().count(), 5);
        // The helium particles
        assert_eq!(
            configuration.molecule(0).hash(),
            configuration.molecule(3).hash()
        );

        // The water molecules
        assert_eq!(
            configuration.molecule(1).hash(),
            configuration.molecule(2).hash()
        );
        assert_ne!(
            configuration.molecule(1).hash(),
            configuration.molecule(4).hash()
        );
    }
}
