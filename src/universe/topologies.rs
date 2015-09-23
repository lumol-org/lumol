/*
 * Cymbalum, Molecular Simulation in Rust
 * Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
*/

//! Data about bonds and angles in the system.
use std::cmp::{min, max};
use std::collections::HashSet;
use std::collections::HashMap;

/// A `Bond` between the atoms at indexes `i` and `j`
///
/// This structure ensure unicity of a `Bond` representation by enforcing
/// `i < j`
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Bond {
    i: usize,
    j: usize,
}

impl Bond {
    /// Create a new Bond between the atoms at indexes `first` and `second`
    pub fn new(first: usize, second: usize) -> Bond {
        assert!(first != second);
        let i = min(first, second);
        let j = max(first, second);
        Bond{i: i, j: j}
    }

    /// Get the first particle in the bond
    #[inline] pub fn i(&self) -> usize {
        self.i
    }

    /// Get the second particle in the bond
    #[inline] pub fn j(&self) -> usize {
        self.j
    }
}

/// An `Angle` formed by the atoms at indexes `i`, `j` and `k`
///
/// This structure ensure unicity of the `Angle` representation by enforcing
/// `i < k`
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Angle {
    i: usize,
    j: usize,
    k: usize,
}

impl Angle {
    /// Create a new Angle between the atoms at indexes `first`, `second` and `third`
    pub fn new(first: usize, second: usize, third: usize) -> Angle {
        assert!(first != second);
        assert!(first != third);
        assert!(second != third);
        let i = min(first, third);
        let k = max(first, third);
        Angle{i: i, j: second, k: k}
    }

    /// Get the first particle in the angle
    #[inline] pub fn i(&self) -> usize {
        self.i
    }

    /// Get the second particle in the angle
    #[inline] pub fn j(&self) -> usize {
        self.j
    }

    /// Get the third particle in the angle
    #[inline] pub fn k(&self) -> usize {
        self.k
    }
}


/// A `Dihedral` angle formed by the atoms at indexes `i`, `j`, `k` and `m`
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Dihedral {
    i: usize,
    j: usize,
    k: usize,
    m: usize,
}

impl Dihedral {
    /// Create a new Dihedral between the atoms at indexes `first`, `second`, `third` and `fourth`
    pub fn new(first: usize, second: usize, third: usize, fourth: usize) -> Dihedral {
        assert!(first != second);
        assert!(second != third);
        assert!(third != fourth);
        let (i, j, k, m) = if max(first, second) < max(third, fourth) {
            (first, second, third, fourth)
        } else {
            (fourth, third, second, first)
        };
        Dihedral{i: i, j: j, k: k, m: m}
    }

    /// Get the first particle in the dihedral angle
    #[inline] pub fn i(&self) -> usize {
        self.i
    }

    /// Get the second particle in the dihedral angle
    #[inline] pub fn j(&self) -> usize {
        self.j
    }

    /// Get the third particle in the dihedral angle
    #[inline] pub fn k(&self) -> usize {
        self.k
    }

    /// Get the fourth particle in the dihedral angle
    #[inline] pub fn m(&self) -> usize {
        self.m
    }
}

/******************************************************************************/
/// The `Molecules` struct is a bidirectional HashMap associating particles indexes and
/// molecules indexes.
pub struct Molecules {
    /// Mapping particles indexes to molecules indexes
    mols: HashMap<usize, usize>,
    /// Mapping molecules indexes to list of particles indexes
    parts: HashMap<usize, Vec<usize>>,
    /// The last molecule we built.
    last_mol: usize,
}

impl Molecules {
    /// Create a new molecule set
    pub fn new() -> Molecules {
        Molecules{
            mols: HashMap::new(),
            parts: HashMap::new(),
            last_mol: 0
        }
    }

    /// Get the molecules as list of lists of particles  -> particles mapping
    #[inline] pub fn molecules<'a>(&'a self) -> Vec<&'a Vec<usize>> {
        self.parts.values().collect()
    }

    /// Get the molecules index of the particle `i`
    #[inline] pub fn molecule_id(&self, i:usize) -> usize {
        debug_assert!(self.mols.contains_key(&i));
        self.mols[&i]
    }

    /// Get the molecules as list of lists of particles  -> particles mapping
    #[inline] pub fn molecule_containing<'a>(&'a self, i:usize) -> &'a Vec<usize> {
        let mol = self.molecule_id(i);
        debug_assert!(self.mols.contains_key(&mol));
        &self.parts[&mol]
    }

    /// Merge the molecules containing the particles `i` and `j` in one molecule
    pub fn merge(&mut self, i: usize, j: usize) {
        debug_assert!(self.mols.contains_key(&i));
        debug_assert!(self.mols.contains_key(&j));
        let mol_i = self.mols[&i];
        let mol_j = self.mols[&j];
        let new_mol = min(mol_i, mol_j);
        let old_mol = max(mol_i, mol_j);
        for (_, mol) in self.mols.iter_mut() {
            if *mol == old_mol {
                *mol = new_mol;
            }
        }
        {
            debug_assert!(self.parts.contains_key(&old_mol));
            debug_assert!(self.parts.contains_key(&new_mol));
            let old_mol_parts = self.parts.get(&old_mol).unwrap().clone();
            let new_mol_parts = self.parts.get_mut(&new_mol).unwrap();
            for part in old_mol_parts.iter() {
                new_mol_parts.push(*part);
            }
        }
        self.parts.remove(&old_mol);
    }

    /// Add a new particle in the internal list, as a molecule with only one atom.
    pub fn add_particle(&mut self, i: usize) {
        self.mols.insert(i, self.last_mol);
        self.parts.insert(self.last_mol, vec![i]);
        self.last_mol += 1;
    }

    /// Remove a particle from the internal list.
    pub fn remove_particle(&mut self, i: usize) {
        self.mols.remove(&i);
        for (_, part_list) in self.parts.iter_mut() {
            match find(part_list, &i) {
                Some(index) => {part_list.swap_remove(index);}
                None => {}
            }
        }
    }
}

// TODO: unstable remove this functions and use Vec::position_elem
fn find(vec: &Vec<usize>, val: &usize) -> Option<usize> {
    for (i, v) in vec.iter().enumerate() {
        if v == val {
            return Some(i);
        }
    }
    return None;
}

/// The `Topology` of a system contains the list of molecules, bonds, angles and dihedral angles
/// in the system.
///
/// The only reliable source of information is the list of bonds. All the other data are cached
/// for efficiency, and rebuilt when needed.
pub struct Topology {
    /// Molecules in the system
    molecules: Molecules,
    /// All the bonds in the system
    bonds: HashSet<Bond>,
    /// All the angles in the system
    angles: HashSet<Angle>,
    /// All the dihedral angles in the system
    dihedrals: HashSet<Dihedral>,
}

impl Topology {
    /// Create a new empty Topology
    pub fn new() -> Topology {
        Topology{
            molecules: Molecules::new(),
            bonds: HashSet::new(),
            angles: HashSet::new(),
            dihedrals: HashSet::new(),
        }
    }

    /// Get the list of molecules in the system as list of lists of particles
    /// indexes.
    #[inline] pub fn molecules<'a>(&'a self) -> Vec<&'a Vec<usize>> {
        self.molecules.molecules()
    }

    /// Get the molecule containing the particle at index `i` as a list of
    /// particles indexes.
    #[inline] pub fn molecule_containing<'a>(&'a self, i:usize) -> &'a Vec<usize> {
        self.molecules.molecule_containing(i)
    }

    /// Check if the particles at indexes `i` and `j` are in the same molecule
    #[inline] pub fn are_in_same_molecule(&self, i: usize, j:usize) -> bool {
        self.molecules.molecule_id(i) == self.molecules.molecule_id(j)
    }

    /// Merge the molecules containing the atoms at indexes `i` and `j` in one molecule
    fn merge_molecules(&mut self, i: usize, j: usize) {
        self.molecules.merge(i, j);
    }

    /// Remove the molecule containing the particle at index `i`
    pub fn remove_molecule_containing(&mut self, i: usize) {
        let to_remove = self.molecules.molecule_containing(i).clone();
        for part in to_remove.iter() {
            self.remove_particle(*part);
        }
    }

    /// Add a bond between the particles at indexes `i` and `j`. This fails if
    /// the particles were not added before calling this function, preferably
    /// through the Universe::add_particle function.
    pub fn add_bond(&mut self, i: usize, j: usize) {
        if !self.are_in_same_molecule(i, j) {
            self.merge_molecules(i, j);
        };
        self.bonds.insert(Bond::new(i, j));
        self.rebuild();
    }

    /// Add a new particle by it's index in the current universe.
    #[inline] pub fn add_particle(&mut self, i: usize) {
        self.molecules.add_particle(i);
    }

    /// Removes particle at index `i` and any associated bonds, angle or dihedral
    pub fn remove_particle(&mut self, i: usize) {
        self.molecules.remove_particle(i);
        // Remove bonds containing the particle `i`
        let mut to_remove = Vec::new();
        for bond in self.bonds.iter() {
            if bond.i == i || bond.j == i {
                to_remove.push(bond.clone());
            }
        }
        for bond in to_remove.iter() {
            self.bonds.remove(bond);
        }
        self.rebuild();
    }

    /// Cleanup cached data
    fn cleanup(&mut self) {
        self.angles.clear();
        self.dihedrals.clear();
    }

    /// Rebuild the full list of angles and dihedral angles from the list of bonds
    fn rebuild(&mut self) {
        self.cleanup();
        for bond1 in self.bonds.iter() {
            // Find angles
            for bond2 in self.bonds.iter() {
                if bond1 == bond2 {
                    continue;
                }

                let angle = if bond1.i == bond2.j {
                    Angle::new(bond2.i, bond2.j, bond1.j)
                } else if bond1.j == bond2.i {
                    Angle::new(bond1.i, bond1.j, bond2.j)
                } else if bond1.j == bond2.j {
                    Angle::new(bond1.i, bond1.j, bond2.i)
                } else if bond1.i == bond2.i {
                    Angle::new(bond1.j, bond1.i, bond2.j)
                } else {
                    // We will not find any dihedral angle from these bonds
                    continue;
                };
                self.angles.insert(angle);

                // Find dihedral angles
                for bond3 in self.bonds.iter() {
                    if bond2 == bond3 {
                        continue;
                    }

                    let dihedral = if angle.k == bond3.i && angle.j != bond3.j {
                        Dihedral::new(angle.i, angle.j, angle.k, bond3.j)
                    } else if angle.k == bond3.i && angle.j != bond3.j {
                        Dihedral::new(angle.i, angle.j, angle.k, bond3.j)
                    } else if angle.i == bond3.j && angle.j != bond3.i {
                        Dihedral::new(bond3.i, angle.i, angle.j, angle.k)
                    } else if angle.i == bond3.i && angle.j != bond3.j {
                        Dihedral::new(bond3.j, angle.i, angle.j, angle.k)
                    } else if angle.k == bond3.i || angle.k == bond3.j {
                        // TODO: this is an improper dihedral
                        continue;
                    } else {
                        continue;
                    };
                    self.dihedrals.insert(dihedral);
                }
            }
        }
    }

    /// Get the internal list of bonds
    #[inline] pub fn bonds<'a>(&'a self) -> &'a HashSet<Bond> {
        &self.bonds
    }

    /// Get the internal list of angles
    #[inline] pub fn angles<'a>(&'a self) -> &'a HashSet<Angle> {
        &self.angles
    }

    /// Get the internal list of dihedrals
    #[inline] pub fn dihedrals<'a>(&'a self) -> &'a HashSet<Dihedral> {
        &self.dihedrals
    }
}


#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn unicity() {
        let bond = Bond::new(2, 1);
        assert_eq!(bond.i, 1);
        assert_eq!(bond.j, 2);

        let angle = Angle::new(8, 7, 6);
        assert_eq!(angle.i, 6);
        assert_eq!(angle.j, 7);
        assert_eq!(angle.k, 8);

        let dihedral = Dihedral::new(8, 7, 6, 0);
        assert_eq!(dihedral.i, 0);
        assert_eq!(dihedral.j, 6);
        assert_eq!(dihedral.k, 7);
        assert_eq!(dihedral.m, 8);

        let dihedral = Dihedral::new(0, 7, 6, 8);
        assert_eq!(dihedral.i, 0);
        assert_eq!(dihedral.j, 7);
        assert_eq!(dihedral.k, 6);
        assert_eq!(dihedral.m, 8);
    }

    fn find_in<T: PartialEq>(vec: &Vec<T>, val: &T) -> bool {
        for v in vec.iter() {
            if v == val {
                return true;
            }
        }
        return false;
    }

    #[test]
    fn molecules() {
        let mut mols = Molecules::new();
        mols.add_particle(0);
        mols.add_particle(1);
        mols.add_particle(2);
        mols.add_particle(3);
        mols.add_particle(4);

        mols.merge(0, 1);
        mols.merge(0, 2);
        mols.merge(3, 4);

        assert_eq!(mols.molecules().len(), 2);

        {
            let all_mols = &mut mols.molecules().clone();
            // Sorting because we have no order guarantee in the HashSet
            all_mols.sort_by(|a, b| a.len().cmp(&b.len()));
            assert_eq!(all_mols.to_vec(), vec![&vec![3, 4], &vec![0, 1, 2]]);

            assert_eq!(mols.molecule_id(4), 3);
            assert_eq!(mols.molecule_containing(4), &vec![3, 4]);
            assert_eq!(mols.molecule_id(2), 0);
            assert_eq!(mols.molecule_containing(2), &vec![0, 1, 2]);
        }

        mols.remove_particle(1);
        assert_eq!(mols.molecules().len(), 2);
        assert_eq!(mols.molecule_containing(4), &vec![3, 4]);
        assert_eq!(mols.molecule_containing(0), &vec![0, 2]);
    }

    #[test]
    fn topology() {
        // Create a 2D ethane by hand, like this one
        //       H    H               4    5
        //       |    |               |    |
        //   H - C -- C - H       3 - 0 -- 1 - 6
        //       |    |               |    |
        //       H    H               2    7

        let mut topology = Topology::new();
        for i in 0..8 {
            topology.add_particle(i);
        }
        topology.add_bond(0, 1);
        topology.add_bond(0, 2);
        topology.add_bond(0, 3);
        topology.add_bond(0, 4);
        topology.add_bond(1, 5);
        topology.add_bond(1, 6);
        topology.add_bond(1, 7);

        /**********************************************************************/
        let bonds = vec![Bond::new(0, 1), Bond::new(0, 2), Bond::new(0, 3),
                         Bond::new(0, 4), Bond::new(1, 5), Bond::new(1, 6),
                         Bond::new(1, 7)];

        assert_eq!(topology.bonds().len(), bonds.len());
        for bond in topology.bonds().iter() {
            assert!(find_in(&bonds, bond))
        }

        /**********************************************************************/
        let angles = vec![Angle::new(0, 1, 5), Angle::new(0, 1, 6), Angle::new(0, 1, 7),
                          Angle::new(1, 0, 2), Angle::new(1, 0, 3), Angle::new(1, 0, 4),
                          Angle::new(2, 0, 3), Angle::new(3, 0, 4), Angle::new(2, 0, 4),
                          Angle::new(5, 1, 6), Angle::new(6, 1, 7), Angle::new(5, 1, 7)];

        assert_eq!(topology.angles().len(), angles.len());
        for angle in topology.angles().iter() {
            assert!(find_in(&angles, angle))
        }

        /**********************************************************************/
        let dihedrals = vec![Dihedral::new(2, 0, 1, 5), Dihedral::new(2, 0, 1, 6),
                             Dihedral::new(2, 0, 1, 7), Dihedral::new(3, 0, 1, 5),
                             Dihedral::new(3, 0, 1, 6), Dihedral::new(3, 0, 1, 7),
                             Dihedral::new(4, 0, 1, 5), Dihedral::new(4, 0, 1, 6),
                             Dihedral::new(4, 0, 1, 7)];

        assert_eq!(topology.dihedrals().len(), dihedrals.len());
        for dihedral in topology.dihedrals().iter() {
            assert!(find_in(&dihedrals, dihedral))
        }

        {
            let molecules = topology.molecules();
            assert_eq!(molecules.len(), 1);
            let mol = topology.molecule_containing(0);
            assert_eq!(mol, molecules[0]);
        }

        assert!(topology.are_in_same_molecule(2, 5));
        assert!(topology.are_in_same_molecule(0, 7));
        assert!(topology.are_in_same_molecule(4, 2));

        topology.remove_particle(6);
        assert_eq!(topology.bonds().len(), 6);
        assert_eq!(topology.angles().len(), 9);
        assert_eq!(topology.dihedrals().len(), 6);

        topology.remove_molecule_containing(2);

        assert_eq!(topology.molecules()[0].len(), 0);
        assert_eq!(topology.bonds().len(), 0);
        assert_eq!(topology.angles().len(), 0);
        assert_eq!(topology.dihedrals().len(), 0);
    }
}
