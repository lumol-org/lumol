// Cymbalum, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

//! Data about bonds and angles in the system.
use std::cmp::{min, max};
use std::collections::HashSet;
use std::iter::IntoIterator;
use std::ops::Range;
use std::hash::{Hash, SipHasher, Hasher};

use types::Array2;
use system::Particle;

/// A `Bond` between the atoms at indexes `i` and `j`
///
/// This structure ensure unicity of a `Bond` representation by enforcing
/// `i < j`
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
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
    #[inline] pub fn i(&self) -> usize {self.i}

    /// Get the second particle in the bond
    #[inline] pub fn j(&self) -> usize {self.j}
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
    #[inline] pub fn i(&self) -> usize {self.i}

    /// Get the second particle in the angle
    #[inline] pub fn j(&self) -> usize {self.j}

    /// Get the third particle in the angle
    #[inline] pub fn k(&self) -> usize {self.k}
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
    #[inline] pub fn i(&self) -> usize {self.i}

    /// Get the second particle in the dihedral angle
    #[inline] pub fn j(&self) -> usize {self.j}

    /// Get the third particle in the dihedral angle
    #[inline] pub fn k(&self) -> usize {self.k}

    /// Get the fourth particle in the dihedral angle
    #[inline] pub fn m(&self) -> usize {self.m}
}

/******************************************************************************/
mod connect {
    #![allow(dead_code)]
    bitflags! {
        /// The `Connectivity` bitflag encode the topological distance between
        /// two particles in the molecule, i.e. the number of bonds between the
        /// particles.
        flags Connectivity: u16 {
            /// The particles are separated by one bond
            const CONNECT_12   = 0b0001,
            /// The particles are separated by two bonds
            const CONNECT_13   = 0b0010,
            /// The particles are separated by three bonds
            const CONNECT_14   = 0b0100,
            /// The particles are separated by more than three bonds
            const CONNECT_FAR  = 0b1000,
        }
    }

    impl Default for Connectivity {
        fn default() -> Connectivity {
            CONNECT_FAR
        }
    }
}

pub use self::connect::Connectivity;
pub use self::connect::{CONNECT_12, CONNECT_13, CONNECT_14, CONNECT_FAR};

/******************************************************************************/

#[derive(Debug, Clone, PartialEq)]
/// A molecule is the basic building block for a topology. It contains data
/// about the connectivity (bonds, angles, dihedrals) in the system.
pub struct Molecule {
    /// All the bonds in the molecule.
    bonds: HashSet<Bond>,
    /// All the angles in the molecule. This is rebuilt as needed from the bond
    /// list.
    angles: HashSet<Angle>,
    /// All the dihedral angles in the molecule. This is rebuilt as needed from
    /// the bond list.
    dihedrals: HashSet<Dihedral>,
    /// Matrix of connectivity in the molecule. The item at index `i, j` encode
    /// the connectivity between the particles `i + self.first` and
    /// `j + self.first`
    connections: Array2<Connectivity>,
    /// Index of the first atom in this molecule.
    first: usize,
    /// Index of the last (included) atom in this molecule.
    last: usize,
    /// Hashed value of the set of bonds in the system
    cached_hash: u64
}

impl Hash for Molecule {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.cached_hash.hash(state);
    }
}

impl Molecule {
    /// Create a new `Molecule` containing only the atom i
    pub fn new(i: usize) -> Molecule {
        Molecule {
            bonds: HashSet::new(),
            angles: HashSet::new(),
            dihedrals: HashSet::new(),
            connections: Array2::default((1, 1)),
            first: i,
            last: i,
            cached_hash: 0
        }
    }

    /// Get the number of atoms in the molecule
    pub fn size(&self) -> usize {
        self.last - self.first + 1
    }

    /// Get the first atom of this molecule
    pub fn first(&self) -> usize {
        self.first
    }

    /// Get the last atom of this molecule
    pub fn last(&self) -> usize {
        self.last
    }

    /// Cache the hash of the bonds
    fn rehash(&mut self) {
        let mut hasher = SipHasher::new();
        (self.last - self.first).hash(&mut hasher);

        let mut bonds = self.bonds.iter()
                              .map(|bond| Bond::new(bond.i() - self.first, bond.j() - self.first))
                              .collect::<Vec<_>>();
        bonds.sort();
        for bond in &bonds {
            bond.i().hash(&mut hasher);
            bond.j().hash(&mut hasher);
        }
        self.cached_hash = hasher.finish();
    }

    /// Cleanup cached data
    fn cleanup(&mut self) {
        self.angles.clear();
        self.dihedrals.clear();
    }

    /// Rebuild the full list of angles and dihedral angles from the list of bonds
    fn rebuild(&mut self) {
        self.cleanup();
        for bond1 in &self.bonds {
            // Find angles
            for bond2 in &self.bonds {
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
                for bond3 in &self.bonds {
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
        self.rebuild_connections();
        self.rehash();
    }

    /// Recompute the connectivity matrix from the bonds, angles and dihedrals
    /// in the system.
    fn rebuild_connections(&mut self) {
        let n = self.size();
        self.connections = Array2::default((n, n));

        // Getting needed variables for the `add_connect_term` closure
        let first = self.first;
        let connections = &mut self.connections;
        let mut add_connect_term = |i, j, term| {
            let old_connect = connections[(i - first, j - first)];
            connections[(i - first, j - first)] = old_connect | term;
        };

        for bond in &self.bonds {
            add_connect_term(bond.i, bond.j, CONNECT_12);
            add_connect_term(bond.j, bond.i, CONNECT_12);
        }

        for angle in &self.angles {
            add_connect_term(angle.i, angle.k, CONNECT_13);
            add_connect_term(angle.k, angle.i, CONNECT_13);
        }

        for dihedral in &self.dihedrals {
            add_connect_term(dihedral.i, dihedral.m, CONNECT_14);
            add_connect_term(dihedral.m, dihedral.i, CONNECT_14);
        }
    }

    /// Merge this molecule with `other`. The first particle in `other` should
    /// be the particle just after the last one in `self`.
    pub fn merge_with(&mut self, other: Molecule) {
        assert!(self.last + 1 == other.first);
        self.last = other.last;
        for bond in other.bonds() {
            self.bonds.insert(bond.clone());
        }

        // Get angles and dihedrals from the other molecule, there is no need to
        // rebuild these.
        for angle in other.angles() {
            self.angles.insert(angle.clone());
        }

        for dihedral in other.dihedrals() {
            self.dihedrals.insert(dihedral.clone());
        }

        self.rebuild_connections();
        self.rehash();
    }

    /// Translate all indexes in this molecule by `delta`.
    pub fn translate_by(&mut self, delta: isize) {
        if delta < 0 {
            // We should not create negative indexes
            assert!((delta.abs() as usize) < self.first);
        }

        // The wrapping_add are necessary here, and produce the right result,
        // thanks to integer overflow and the conversion below.
        let delta = delta as usize;
        self.first = self.first.wrapping_add(delta);
        self.last = self.last.wrapping_add(delta);

        let mut new_bonds = HashSet::new();
        for bond in &self.bonds {
            new_bonds.insert(Bond::new(
                bond.i().wrapping_add(delta),
                bond.j().wrapping_add(delta)
            ));
        }
        self.bonds = new_bonds;

        let mut new_angles = HashSet::new();
        for angle in &self.angles {
            new_angles.insert(Angle::new(
                angle.i().wrapping_add(delta),
                angle.j().wrapping_add(delta),
                angle.k().wrapping_add(delta)
            ));
        }
        self.angles = new_angles;

        let mut new_dihedrals = HashSet::new();
        for dihedral in &self.dihedrals {
            new_dihedrals.insert(Dihedral::new(
                dihedral.i().wrapping_add(delta),
                dihedral.j().wrapping_add(delta),
                dihedral.k().wrapping_add(delta),
                dihedral.m().wrapping_add(delta)
            ));
        }
        self.dihedrals = new_dihedrals;
    }

    /// Add a bond between the particles at indexes `i` and `j`. These particles
    /// are assumed to be in the molecule
    pub fn add_bond(&mut self, i: usize, j: usize)  {
        assert!(self.first <= i && i <= self.last);
        assert!(self.first <= j && j <= self.last);
        self.bonds.insert(Bond::new(i, j));
        self.rebuild();
    }

    /// Removes particle at index `i` and any associated bonds, angle or
    /// dihedral. This function also update the indexes for the
    /// bonds/angles/dihedral by remove 1 to all the values `> i`
    pub fn remove_particle(&mut self, i: usize) {
        assert!(self.first <= i && i <= self.last);
        // Remove bonds containing the particle `i`
        let mut new_bonds = HashSet::new();
        for bond in self.bonds() {
            if bond.i == i || bond.j == i {
                continue;
            }

            let mut new_bond = bond.clone();
            if new_bond.i > i {
                new_bond.i -= 1;
            }
            if new_bond.j > i {
                new_bond.j -= 1;
            }

            new_bonds.insert(new_bond);
        }

        self.bonds = new_bonds;
        self.last -= 1;
        self.rebuild();
    }

    /// Get the internal list of bonds
    #[inline] pub fn bonds(&self) -> &HashSet<Bond> {
        &self.bonds
    }

    /// Get the internal list of angles
    #[inline] pub fn angles(&self) -> &HashSet<Angle> {
        &self.angles
    }

    /// Get the internal list of dihedrals
    #[inline] pub fn dihedrals(&self) -> &HashSet<Dihedral> {
        &self.dihedrals
    }

    /// Get the connectivity betweent the particles `i` and `j`
    #[inline] pub fn connectivity(&self, i: usize, j: usize) -> Connectivity {
        assert!(self.first <= i && i <= self.last);
        assert!(self.first <= j && j <= self.last);
        return self.connections[(i - self.first, j - self.first)];
    }

    /// Get an iterator over the particles in the molecule
    #[inline] pub fn iter(&self) -> Range<usize> {
        self.into_iter()
    }
}

/// Get the molecule type of the given `molecule` containing the `particles`.
/// This type can be used to identify all the molecules containing the same
/// bonds and particles (see `System::moltype` for more information).
pub fn moltype(molecule: &Molecule, particles: &[Particle]) -> u64 {
    assert!(particles.len() == molecule.size());
    let mut hasher = SipHasher::new();
    molecule.cached_hash.hash(&mut hasher);
    for particle in particles {
        particle.name().hash(&mut hasher);
    }
    hasher.finish()
}

impl IntoIterator for Molecule {
    type Item = usize;
    type IntoIter = Range<usize>;

    fn into_iter(self) -> Range<usize> {
        Range{
            start: self.first,
            end: self.last + 1,
        }
    }
}

impl<'a> IntoIterator for &'a Molecule {
    type Item = usize;
    type IntoIter = Range<usize>;

    fn into_iter(self) -> Range<usize> {
        Range{
            start: self.first,
            end: self.last + 1,
        }
    }
}

/******************************************************************************/
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

    #[test]
    fn translate() {
        let mut molecule = Molecule::new(0);
        for i in 1..4 {
            molecule.merge_with(Molecule::new(i));
        }
        molecule.add_bond(0, 1);
        molecule.add_bond(1, 2);
        molecule.add_bond(2, 3);

        assert_eq!(molecule.first(), 0);
        assert_eq!(molecule.last(), 3);
        assert!(molecule.bonds().contains(&Bond::new(2, 3)));
        assert!(molecule.angles().contains(&Angle::new(1, 2, 3)));
        assert!(molecule.dihedrals().contains(&Dihedral::new(0, 1, 2, 3)));

        molecule.translate_by(5);
        assert_eq!(molecule.first(), 5);
        assert_eq!(molecule.last(), 8);
        assert!(molecule.bonds().contains(&Bond::new(7, 8)));
        assert!(molecule.angles().contains(&Angle::new(6, 7, 8)));
        assert!(molecule.dihedrals().contains(&Dihedral::new(5, 6, 7, 8)));

        molecule.translate_by(-3);
        assert_eq!(molecule.first(), 2);
        assert_eq!(molecule.last(), 5);
        assert!(molecule.bonds().contains(&Bond::new(4, 5)));
        assert!(molecule.angles().contains(&Angle::new(3, 4, 5)));
        assert!(molecule.dihedrals().contains(&Dihedral::new(2, 3, 4, 5)));
    }

    #[test]
    fn bonding() {
        // Create ethane like this
        //       H    H               4    5
        //       |    |               |    |
        //   H - C -- C - H       3 - 0 -- 1 - 6
        //       |    |               |    |
        //       H    H               2    7
        let mut molecule = Molecule::new(0);

        for i in 1..8 {
            molecule.merge_with(Molecule::new(i));
        }

        molecule.add_bond(0, 1);
        molecule.add_bond(0, 2);
        molecule.add_bond(0, 3);
        molecule.add_bond(0, 4);
        molecule.add_bond(1, 5);
        molecule.add_bond(1, 6);
        molecule.add_bond(1, 7);

        assert_eq!(molecule.first(), 0);
        assert_eq!(molecule.last(), 7);

        assert_eq!(molecule.size(), 8);

        /**********************************************************************/
        let bonds = vec![Bond::new(0, 1), Bond::new(0, 2), Bond::new(0, 3),
                         Bond::new(0, 4), Bond::new(1, 5), Bond::new(1, 6),
                         Bond::new(1, 7)];

        assert_eq!(molecule.bonds().len(), bonds.len());
        for bond in &bonds {
            assert!(molecule.bonds().contains(bond));
        }

        /**********************************************************************/
        let angles = vec![Angle::new(0, 1, 5), Angle::new(0, 1, 6), Angle::new(0, 1, 7),
                          Angle::new(1, 0, 2), Angle::new(1, 0, 3), Angle::new(1, 0, 4),
                          Angle::new(2, 0, 3), Angle::new(3, 0, 4), Angle::new(2, 0, 4),
                          Angle::new(5, 1, 6), Angle::new(6, 1, 7), Angle::new(5, 1, 7)];

        assert_eq!(molecule.angles().len(), angles.len());
        for angle in &angles {
            assert!(molecule.angles().contains(angle));
        }

        /**********************************************************************/
        let dihedrals = vec![Dihedral::new(2, 0, 1, 5), Dihedral::new(2, 0, 1, 6),
                             Dihedral::new(2, 0, 1, 7), Dihedral::new(3, 0, 1, 5),
                             Dihedral::new(3, 0, 1, 6), Dihedral::new(3, 0, 1, 7),
                             Dihedral::new(4, 0, 1, 5), Dihedral::new(4, 0, 1, 6),
                             Dihedral::new(4, 0, 1, 7)];

        assert_eq!(molecule.dihedrals().len(), dihedrals.len());
        for dihedral in &dihedrals {
            assert!(molecule.dihedrals().contains(dihedral));
        }

        /**********************************************************************/
        assert!(molecule.connectivity(0, 1).contains(CONNECT_12));
        assert!(molecule.connectivity(1, 0).contains(CONNECT_12));

        assert!(molecule.connectivity(0, 7).contains(CONNECT_13));
        assert!(molecule.connectivity(7, 0).contains(CONNECT_13));

        assert!(molecule.connectivity(3, 5).contains(CONNECT_14));
        assert!(molecule.connectivity(5, 3).contains(CONNECT_14));

        /**********************************************************************/

        molecule.remove_particle(6);
        assert_eq!(molecule.bonds().len(), 6);
        assert_eq!(molecule.angles().len(), 9);
        assert_eq!(molecule.dihedrals().len(), 6);
    }

    #[test]
    fn cyclic() {
        //   0 -- 1
        //   |    |
        //   3 -- 2
        let mut molecule = Molecule::new(0);

        for i in 1..4 {
            molecule.merge_with(Molecule::new(i));
        }
        molecule.add_bond(0, 1);
        molecule.add_bond(1, 2);
        molecule.add_bond(2, 3);
        molecule.add_bond(3, 0);

        assert!(molecule.connectivity(0, 3).contains(CONNECT_12));
        assert!(molecule.connectivity(0, 3).contains(CONNECT_14));

        assert!(molecule.angles.contains(&Angle::new(0, 3, 2)));
        assert!(molecule.angles.contains(&Angle::new(0, 1, 2)));
    }

    #[test]
    fn remove_particle() {
        let mut molecule = Molecule::new(0);
        for i in 1..4 {
            molecule.merge_with(Molecule::new(i));
        }

        assert_eq!(molecule.bonds().len(), 0);
        assert_eq!(molecule.size(), 4);

        molecule.add_bond(0, 1);
        molecule.add_bond(2, 3);
        assert_eq!(molecule.bonds().len(), 2);
        assert_eq!(molecule.size(), 4);

        molecule.remove_particle(1);
        assert_eq!(molecule.bonds().len(), 1);
        assert_eq!(molecule.size(), 3);

        molecule.merge_with(Molecule::new(3));
        assert_eq!(molecule.bonds().len(), 1);
        assert_eq!(molecule.size(), 4);
    }

    #[test]
    fn hash() {
        let mut molecule = Molecule::new(0);
        assert_eq!(molecule.cached_hash, 0);

        for i in 1..4 {
            molecule.merge_with(Molecule::new(i));
        }

        let hash = molecule.cached_hash;
        assert!(hash != 0);

        molecule.add_bond(0, 1);
        molecule.add_bond(2, 3);
        assert!(molecule.cached_hash != hash);
        let hash = molecule.cached_hash;

        molecule.remove_particle(1);
        assert!(molecule.cached_hash != hash);
        let hash = molecule.cached_hash;

        // Hash should be the same when translating the molecule
        molecule.translate_by(67);
        molecule.rehash();
        assert_eq!(molecule.cached_hash, hash);
    }
}
