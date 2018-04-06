// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 Lumol's contributors — BSD license

//! The Configuration type definition

use std::cmp::{max, min};
use std::i8;

use types::{Vector3D, Zero};

use sys::{BondDistance, Molecule, ParticleKind, UnitCell};
use sys::{Particle, ParticleSlice, ParticleSliceMut, ParticleVec};
use sys::molecule_type;

/// Particles permutations:. Indexes are given in the `(old, new)` form.
pub type Permutations = Vec<(usize, usize)>;

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
    /// Molecules in the system
    molecules: Vec<Molecule>,
    /// Molecules indexes for all the particles
    molids: Vec<usize>,
}

impl Configuration {
    /// Create a new empty `Configuration`
    pub fn new() -> Configuration {
        Configuration {
            particles: ParticleVec::new(),
            molecules: Vec::new(),
            molids: Vec::new(),
            cell: UnitCell::new(),
        }
    }
}

/// Topology and particles related functions
impl Configuration {
    /// Get the type of the molecule at index `molid`. This type is a hash of
    /// the atoms names, and the set of bonds in the molecule. This means that
    /// two molecules will have the same type if and only if they contains the
    /// same atoms and the same bonds, **in the same order**.
    pub fn molecule_type(&self, molid: usize) -> u64 {
        let molecule = self.molecule(molid);
        molecule_type(molecule, self.particles.slice(molecule.into_iter()))
    }

    /// Get a list of molecules with `moltype` molecule type.
    pub fn molecules_with_moltype(&self, moltype: u64) -> Vec<usize> {
        let mut res = Vec::new();
        for i in 0..self.molecules().len() {
            if self.molecule_type(i) == moltype {
                res.push(i);
            }
        }
        return res;
    }

    /// Check if the particles at indexes `i` and `j` are in the same molecule
    #[inline]
    pub fn are_in_same_molecule(&self, i: usize, j: usize) -> bool {
        debug_assert_eq!(self.molids.len(), self.particles.len());
        self.molids[i] == self.molids[j]
    }

    /// Get the list of molecules in the configuration.
    #[inline]
    pub fn molecules(&self) -> &[Molecule] {
        &self.molecules
    }

    /// Get the molecule at index `id`
    #[inline]
    pub fn molecule(&self, id: usize) -> &Molecule {
        &self.molecules[id]
    }

    /// Get the index of the molecule containing the particle `i`
    #[inline]
    pub fn molid(&self, i: usize) -> usize {
        self.molids[i]
    }

    /// Get the length of the shortest bond path to go from the particle `i` to
    /// the particle `j`. If the particles are not in the same molecule, the
    /// length is -1. Else, this length is 0 if `i == j`, 1 if there is a bond
    /// between `i` and `j`, etc.
    pub fn bond_distance(&self, i: usize, j: usize) -> i8 {
        assert!(i < self.size() && j < self.size());
        if !(self.are_in_same_molecule(i, j)) {
            -1
        } else if i == j {
            0
        } else {
            let connect = self.molecule(self.molid(i)).bond_distance(i, j);
            if connect.contains(BondDistance::ONE) {
                1
            } else if connect.contains(BondDistance::TWO) {
                2
            } else if connect.contains(BondDistance::THREE) {
                3
            } else if connect.contains(BondDistance::FAR) {
                i8::MAX
            } else {
                unreachable!();
            }
        }
    }

    /// Remove the molecule at index `i`
    pub fn remove_molecule(&mut self, molid: usize) {
        let molecule = self.molecules.remove(molid);
        let first = molecule.start();
        let size = molecule.size();

        for _ in 0..size {
            let _ = self.particles.remove(first);
            let _ = self.molids.remove(first);
        }

        for molecule in self.molecules.iter_mut().skip(molid) {
            molecule.translate_by(-(size as isize));
        }

        for molid in self.molids.iter_mut().skip(first) {
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
    pub fn add_bond(&mut self, mut particle_i: usize, mut particle_j: usize) -> Permutations {
        assert!(particle_i <= self.particles.len());
        assert!(particle_j <= self.particles.len());
        assert_ne!(particle_i, particle_j);
        trace!(
            "Adding bond {}-{} between molecules {} and {}",
            particle_i,
            particle_j,
            self.molids[particle_i],
            self.molids[particle_j]
        );

        // Getting copy of the molecules before the merge
        let molid_i = self.molids[particle_i];
        let molid_j = self.molids[particle_j];
        let new_molid = min(molid_i, molid_j);
        let old_molid = max(molid_i, molid_j);
        let new_mol = self.molecules[new_molid].clone();
        let old_mol = self.molecules[old_molid].clone();
        let already_in_same_molecule = self.are_in_same_molecule(particle_i, particle_j);

        // Effective merge
        let delta = self.merge_molecules(molid_i, molid_j);

        let mut permutations = Permutations::new();
        // If new_mol.last() + 1 == old_mol.first(), no one moved. Else,
        // we generate the permutations
        if !already_in_same_molecule && new_mol.end() != old_mol.start() {
            let size = old_mol.size();
            let first = old_mol.start();
            let second = new_mol.end();
            // Add permutation for the molecule we just moved around
            for i in 0..size {
                permutations.push((first + i, second + i));
            }

            // Add permutations for molecules that where shifted to make
            // space for the just moved molecule.
            for molecule in &self.molecules[new_molid + 1..old_molid] {
                for i in molecule {
                    permutations.push((i - size, i));
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

        assert_eq!(self.molids[particle_i], self.molids[particle_j]);
        self.molecules[self.molids[particle_i]].add_bond(particle_i, particle_j);
        return permutations;
    }

    /// Removes particle at index `i` and any associated bonds, angle or dihedral
    pub fn remove_particle(&mut self, i: usize) {
        let id = self.molids[i];
        self.molecules[id].remove_particle(i);

        for molecule in self.molecules.iter_mut().skip(id + 1) {
            molecule.translate_by(-1);
        }

        let _ = self.particles.remove(i);
        let _ = self.molids.remove(i);
    }

    /// Insert a particle at the end of the internal list. The new particle
    /// must have a valid particle kind.
    pub fn add_particle(&mut self, particle: Particle) {
        assert_ne!(particle.kind, ParticleKind::invalid());
        self.particles.push(particle);
        self.molecules.push(Molecule::new(self.particles.len() - 1));
        self.molids.push(self.molecules.len() - 1);
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

    /// Return the center-of-mass of a molecule
    ///
    /// # Warning
    ///
    /// This function does not check for the particles' positions' nearest
    /// images. To use this function properly, make sure that all particles of
    /// the molecule are adjacent.
    pub fn molecule_com(&self, molid: usize) -> Vector3D {
        let mut total_mass = 0.0;
        let mut com = Vector3D::zero();
        for i in self.molecule(molid) {
            total_mass += self.particles.mass[i];
            com += self.particles.mass[i] * self.particles.position[i];
        }
        com / total_mass
    }

    /// Return the center-of-mass of the configuration
    ///
    /// # Warning
    ///
    /// This function does not check for the particles' positions' nearest
    /// images. To use this function properly, make sure that all particles of
    /// a molecule are adjacent.
    pub fn center_of_mass(&self) -> Vector3D {
        let mut total_mass = 0.0;
        let mut com = Vector3D::zero();
        for i in 0..self.size() {
            total_mass += self.particles.mass[i];
            com += self.particles.mass[i] * self.particles.position[i];
        }
        com / total_mass
    }

    /// Move all particles of a molecule such that the molecules center-of-mass
    /// position resides inside the simulation cell.
    ///
    /// # Note
    ///
    /// If the `CellShape` is `Infinite` there are no changes
    /// to the positions.
    pub fn wrap_molecule(&mut self, molid: usize) {
        let com = self.molecule_com(molid);
        let mut com_wrapped = com;
        self.cell.wrap_vector(&mut com_wrapped);
        let delta = com_wrapped - com;
        // iterate over all positions and move them accordingly
        for i in self.molecule(molid) {
            self.particles.position[i] += delta;
        }
    }

    /// Get the list of particles in this configuration, as a `ParticleSlice`.
    pub fn particles(&self) -> ParticleSlice {
        self.particles.as_slice()
    }

    /// Get the list of particles in this configuration, as a mutable
    /// `ParticleSliceMut`.
    pub fn particles_mut(&mut self) -> ParticleSliceMut {
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

        let mut new_mol = self.molecules[new_molid].clone();
        let old_mol = self.molecules[old_molid].clone();

        if new_mol.end() == old_mol.start() {
            // Just update the molecules ids
            for i in old_mol {
                self.molids[i] = new_molid;
            }
        } else {
            // Move the particles close together
            let mut new_index = new_mol.end();
            for i in old_mol {
                // Remove particles from the old position, and insert it to the
                // new one. The indexes are valid during the movement, because
                // we insert a new particle for each particle removed.
                let particle = self.particles.remove(i);
                self.particles.insert(new_index, particle);

                // Update molids
                let _ = self.molids.remove(i);
                self.molids.insert(new_index, new_molid);

                new_index += 1;
            }
        }

        let mut old_mol = self.molecules[old_molid].clone();
        let size = old_mol.size() as isize;

        // translate all indexes in the molecules between new_mol and old_mol
        for molecule in
            self.molecules.iter_mut().skip(new_molid + 1).take(old_molid - new_molid - 1)
        {
            molecule.translate_by(size);
        }

        // Update molid for all particles after the old molecule
        for molid in self.molids.iter_mut().skip(old_mol.end()) {
            *molid -= 1;
        }

        let delta = old_mol.start() - new_mol.end();
        old_mol.translate_by(-(delta as isize));

        new_mol.merge_with(old_mol);
        self.molecules[new_molid] = new_mol;
        let _ = self.molecules.remove(old_molid);

        debug_assert!(check_molid_sorted(&self.molids), "Unsorted molecule ids {:?}", self.molids);

        return delta;
    }
}

/// Check that molids is sorted and only contains successive values
fn check_molid_sorted(molids: &[usize]) -> bool {
    let mut previous = 0;
    for &i in molids {
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
        let a = self.particles.position[i];
        let b = self.particles.position[j];
        let c = self.particles.position[k];
        self.cell.angle(&a, &b, &c)
    }

    /// Get the angle and the derivatives of the angle between the particles
    /// `i`, `j` and `k`
    pub fn angle_and_derivatives(
        &self,
        i: usize,
        j: usize,
        k: usize,
    ) -> (f64, Vector3D, Vector3D, Vector3D) {
        let a = self.particles.position[i];
        let b = self.particles.position[j];
        let c = self.particles.position[k];
        self.cell.angle_and_derivatives(&a, &b, &c)
    }

    /// Get the dihedral angle between the particles `i`, `j`, `k` and `m`
    pub fn dihedral(&self, i: usize, j: usize, k: usize, m: usize) -> f64 {
        let a = self.particles.position[i];
        let b = self.particles.position[j];
        let c = self.particles.position[k];
        let d = self.particles.position[m];
        self.cell.dihedral(&a, &b, &c, &d)
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
        let a = self.particles.position[i];
        let b = self.particles.position[j];
        let c = self.particles.position[k];
        let d = self.particles.position[m];
        self.cell.dihedral_and_derivatives(&a, &b, &c, &d)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use sys::{Angle, Bond, Dihedral};
    use types::Vector3D;

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
        configuration.add_particle(particle("H"));
        configuration.add_particle(particle("O"));
        configuration.add_particle(particle("O"));
        configuration.add_particle(particle("H"));

        assert_eq!(configuration.add_bond(0, 1), vec![]);
        assert_eq!(configuration.add_bond(2, 3), vec![]);

        assert_eq!(configuration.molecules().len(), 2);

        let molecule = configuration.molecule(0).clone();
        assert!(molecule.bonds().contains(&Bond::new(0, 1)));
        let molecule = configuration.molecule(1).clone();
        assert!(molecule.bonds().contains(&Bond::new(2, 3)));

        assert_eq!(configuration.add_bond(1, 2), vec![]);
        assert_eq!(configuration.molecules().len(), 1);

        let molecule = configuration.molecule(0).clone();
        assert!(molecule.angles().contains(&Angle::new(0, 1, 2)));
        assert!(molecule.angles().contains(&Angle::new(1, 2, 3)));
        assert!(molecule.dihedrals().contains(&Dihedral::new(0, 1, 2, 3)));

        configuration.remove_particle(2);
        assert_eq!(configuration.molecules().len(), 1);

        let molid = configuration.molid(1);
        configuration.remove_molecule(molid);
        assert_eq!(configuration.molecules().len(), 0);
        assert_eq!(configuration.size(), 0);
    }

    #[test]
    fn add_bonds() {
        // This is a regression test for issue #76
        let mut configuration = Configuration::new();
        configuration.add_particle(particle("H"));
        configuration.add_particle(particle("H"));
        configuration.add_particle(particle("O"));

        assert_eq!(configuration.add_bond(0, 2), vec![(2, 1), (1, 2)]);
        assert_eq!(configuration.add_bond(2, 1), vec![]);
        assert_eq!(configuration.molecules().len(), 1);
    }

    #[test]
    fn bond_distance() {
        let mut configuration = Configuration::new();
        configuration.add_particle(particle("C"));
        configuration.add_particle(particle("C"));
        configuration.add_particle(particle("C"));
        configuration.add_particle(particle("C"));
        configuration.add_particle(particle("C"));
        configuration.add_particle(particle("Zn"));

        let _ = configuration.add_bond(0, 1);
        let _ = configuration.add_bond(1, 2);
        let _ = configuration.add_bond(2, 3);
        let _ = configuration.add_bond(3, 4);

        assert_eq!(configuration.bond_distance(0, 0), 0);
        assert_eq!(configuration.bond_distance(0, 1), 1);
        assert_eq!(configuration.bond_distance(0, 2), 2);
        assert_eq!(configuration.bond_distance(0, 3), 3);
        assert!(configuration.bond_distance(0, 4) > 3);
        assert_eq!(configuration.bond_distance(0, 5), -1);
    }

    #[test]
    fn molecules_permutations() {
        let mut configuration = Configuration::new();
        configuration.add_particle(particle("C"));
        configuration.add_particle(particle("H"));
        configuration.add_particle(particle("H"));
        configuration.add_particle(particle("H"));

        configuration.add_particle(particle("C"));
        configuration.add_particle(particle("H"));
        configuration.add_particle(particle("H"));
        configuration.add_particle(particle("H"));

        assert_eq!(configuration.add_bond(0, 3), vec![(3, 1), (1, 2), (2, 3)]);
        assert_eq!(configuration.add_bond(0, 3), vec![(3, 2), (2, 3)]);
        assert_eq!(configuration.add_bond(0, 3), vec![]);

        assert_eq!(configuration.add_bond(4, 5), vec![]);
        assert_eq!(configuration.add_bond(4, 7), vec![(7, 6), (6, 7)]);
        assert_eq!(configuration.add_bond(4, 7), vec![]);
    }

    #[test]
    fn particles() {
        let mut configuration = Configuration::new();
        configuration.add_particle(particle("O"));
        configuration.add_particle(particle("H"));
        configuration.add_particle(particle("H"));

        assert_eq!(configuration.size(), 3);
        assert_eq!(configuration.particles().name[0], "O");
        assert_eq!(configuration.particles().name[1], "H");
        assert_eq!(configuration.particles().name[2], "H");
    }

    #[test]
    fn distances() {
        let mut configuration = Configuration::new();
        configuration.cell = UnitCell::cubic(5.0);
        configuration.add_particle(particle("O"));
        configuration.add_particle(particle("H"));

        configuration.particles_mut().position[0] = Vector3D::new(9.0, 0.0, 0.0);
        configuration.particles_mut().position[1] = Vector3D::zero();
        assert_eq!(configuration.distance(0, 1), 1.0);

        configuration.cell = UnitCell::new();
        assert_eq!(configuration.distance(0, 1), 9.0);
    }

    #[test]
    fn center_of_mass() {
        let mut configuration = Configuration::new();
        configuration.cell = UnitCell::cubic(5.0);
        configuration.add_particle(particle("O"));
        configuration.add_particle(particle("O"));
        let _ = configuration.add_bond(0, 1);

        configuration.particles_mut().position[0] = Vector3D::new(1.0, 0.0, 0.0);
        configuration.particles_mut().position[1] = Vector3D::zero();
        assert_eq!(configuration.molecule_com(0), Vector3D::new(0.5, 0.0, 0.0));
        assert_eq!(configuration.center_of_mass(), Vector3D::new(0.5, 0.0, 0.0));
    }

    #[test]
    fn test_wrap_molecule() {
        let mut configuration = Configuration::new();
        configuration.cell = UnitCell::cubic(5.0);
        configuration.add_particle(particle("O"));
        configuration.add_particle(particle("O"));
        let _ = configuration.add_bond(0, 1);

        configuration.particles_mut().position[0] = Vector3D::new(-2.0, 0.0, 0.0);
        configuration.particles_mut().position[1] = Vector3D::zero();
        configuration.wrap_molecule(0);

        assert_eq!(configuration.particles().position[0], Vector3D::new(3.0, 0.0, 0.0));
        assert_eq!(configuration.particles().position[1], Vector3D::new(5.0, 0.0, 0.0));
        assert_eq!(configuration.molecule_com(0), Vector3D::new(4.0, 0.0, 0.0))
    }

    #[test]
    fn moltype() {
        let mut configuration = Configuration::new();
        // One helium
        configuration.add_particle(particle("He"));
        // Two water molecules
        configuration.add_particle(particle("H"));
        configuration.add_particle(particle("O"));
        configuration.add_particle(particle("H"));
        let _ = configuration.add_bond(1, 2);
        let _ = configuration.add_bond(2, 3);
        configuration.add_particle(particle("H"));
        configuration.add_particle(particle("O"));
        configuration.add_particle(particle("H"));
        let _ = configuration.add_bond(4, 5);
        let _ = configuration.add_bond(5, 6);
        // Another helium
        configuration.add_particle(particle("He"));
        // A water molecules, with different atoms order
        configuration.add_particle(particle("O"));
        configuration.add_particle(particle("H"));
        configuration.add_particle(particle("H"));
        let _ = configuration.add_bond(8, 9);
        let _ = configuration.add_bond(8, 10);

        assert_eq!(configuration.molecules().len(), 5);
        // The helium particles
        assert_eq!(configuration.molecule_type(0), configuration.molecule_type(3));

        // The water molecules
        assert!(configuration.molecule_type(1) == configuration.molecule_type(2));
        assert!(configuration.molecule_type(1) != configuration.molecule_type(4));
    }
}
