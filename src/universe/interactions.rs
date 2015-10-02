/*
 * Cymbalum, Molecular Simulation in Rust
 * Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
*/
use std::collections::BTreeMap;
use std::cmp::max;

use ::potentials::{PairPotential, AnglePotential, DihedralPotential};

/// The Interaction type hold all data about the potentials in the system,
/// indexed by particle type.
pub struct Interactions {
    /// Pair potentials
    pairs: BTreeMap<(u16, u16), Vec<Box<PairPotential>>>,
    /// Bond potentials
    bonds: BTreeMap<(u16, u16), Vec<Box<PairPotential>>>,
    /// Angle potentials
    angles: BTreeMap<(u16, u16, u16), Vec<Box<AnglePotential>>>,
    /// Dihedral angles potentials
    dihedrals: BTreeMap<(u16, u16, u16, u16), Vec<Box<DihedralPotential>>>,
}

impl Interactions {
    pub fn new() -> Interactions {
        Interactions{
            pairs: BTreeMap::new(),
            bonds: BTreeMap::new(),
            angles: BTreeMap::new(),
            dihedrals: BTreeMap::new(),
        }
    }

    /// Add a pair interaction to the pair `(i, j)`
    pub fn add_pair<T>(&mut self, i: u16, j:u16, potential: T) where T: PairPotential + 'static {
        let (i, j) = sort_pair(i, j);
        if !self.pairs.contains_key(&(i, j)) {
            self.pairs.insert((i, j), Vec::new());
        }
        let pairs = self.pairs.get_mut(&(i, j)).unwrap();
        pairs.push(Box::new(potential));
    }

    /// Get all pair interactions corresponding to the pair `(i, j)`
    pub fn pairs<'a>(&'a self, i: u16, j:u16) -> Option<&'a Vec<Box<PairPotential>>> {
        let (i, j) = sort_pair(i, j);
        self.pairs.get(&(i, j))
    }

    /// Add a bonded interaction to the pair `(i, j)`
    pub fn add_bond<T>(&mut self, i: u16, j:u16, potential: T) where T: PairPotential + 'static {
        let (i, j) = sort_pair(i, j);
        if !self.bonds.contains_key(&(i, j)) {
            self.bonds.insert((i, j), Vec::new());
        }
        let bonds = self.bonds.get_mut(&(i, j)).unwrap();
        bonds.push(Box::new(potential));
    }

    /// Get all bonded interactions corresponding to the pair `(i, j)`
    pub fn bonds<'a>(&'a self, i: u16, j:u16) -> Option<&'a Vec<Box<PairPotential>>> {
        let (i, j) = sort_pair(i, j);
        self.bonds.get(&(i, j))
    }

    /// Add an angle interaction to the angle `(i, j, k)`
    pub fn add_angle<T>(&mut self, i: u16, j:u16, k:u16, potential: T) where T: AnglePotential + 'static {
        let (i, j, k) = sort_angle(i, j, k);
        if !self.angles.contains_key(&(i, j, k)) {
            self.angles.insert((i, j, k), Vec::new());
        }
        let angles = self.angles.get_mut(&(i, j, k)).unwrap();
        angles.push(Box::new(potential));
    }

    /// Get all angle interactions corresponding to the angle `(i, j, k)`
    pub fn angles<'a>(&'a self, i: u16, j:u16, k:u16) -> Option<&'a Vec<Box<AnglePotential>>> {
        let (i, j, k) = sort_angle(i, j, k);
        self.angles.get(&(i, j, k))
    }

    /// Add a dihedral interaction to the dihedral `(i, j, k, m)`
    pub fn add_dihedral<T>(&mut self, i: u16, j:u16, k:u16, m:u16, potential: T) where T: DihedralPotential + 'static {
        let (i, j, k, m) = sort_dihedral(i, j, k, m);
        if !self.dihedrals.contains_key(&(i, j, k, m)) {
            self.dihedrals.insert((i, j, k, m), Vec::new());
        }
        let dihedrals = self.dihedrals.get_mut(&(i, j, k, m)).unwrap();
        dihedrals.push(Box::new(potential));
    }

    /// Get all dihedral interactions corresponding to the dihedral `(i, j, k, m)`
    pub fn dihedrals<'a>(&'a self, i: u16, j:u16, k:u16, m:u16) -> Option<&'a Vec<Box<DihedralPotential>>> {
        let (i, j, k, m) = sort_dihedral(i, j, k, m);
        self.dihedrals.get(&(i, j, k, m))
    }
}

/// Sort pair indexes to get a cannonical representation
#[inline] fn sort_pair(i: u16, j:u16) -> (u16, u16) {
    if i < j {
        (i, j)
    } else {
        (j, i)
    }
}

/// Sort angle indexes to get a cannonical representation
#[inline] fn sort_angle(i: u16, j:u16, k:u16) -> (u16, u16, u16) {
    if i < k {
        (i, j, k)
    } else {
        (k, j, i)
    }
}

/// Sort dihedral indexes to get a cannonical representation
#[inline] fn sort_dihedral(i: u16, j:u16, k:u16, m:u16) -> (u16, u16, u16, u16) {
    if max(i, j) < max(k, m) {
        (i, j, k, m)
    } else {
        (m, k, j, i)
    }
}


#[cfg(test)]
mod test {
    use super::*;
    use ::potentials::Harmonic;

    #[test]
    fn pairs() {
        let mut interactions = Interactions::new();

        interactions.add_pair(0, 3, Harmonic{x0: 0.0, k: 0.0});
        assert_eq!(interactions.pairs(0, 3).unwrap().len(), 1);
        assert_eq!(interactions.pairs(3, 0).unwrap().len(), 1);

        interactions.add_pair(0, 0, Harmonic{x0: 0.0, k: 0.0});
        assert_eq!(interactions.pairs(0, 0).unwrap().len(), 1);
    }

    #[test]
    fn bonds() {
        let mut interactions = Interactions::new();

        interactions.add_bond(0, 3, Harmonic{x0: 0.0, k: 0.0});
        assert_eq!(interactions.bonds(0, 3).unwrap().len(), 1);
        assert_eq!(interactions.bonds(3, 0).unwrap().len(), 1);

        interactions.add_bond(0, 0, Harmonic{x0: 0.0, k: 0.0});
        assert_eq!(interactions.bonds(0, 0).unwrap().len(), 1);
    }

    #[test]
    fn angles() {
        let mut interactions = Interactions::new();

        interactions.add_angle(0, 3, 7, Harmonic{x0: 0.0, k: 0.0});
        assert_eq!(interactions.angles(0, 3, 7).unwrap().len(), 1);
        assert_eq!(interactions.angles(7, 3, 0).unwrap().len(), 1);

        interactions.add_angle(0, 0, 7, Harmonic{x0: 0.0, k: 0.0});
        assert_eq!(interactions.angles(0, 0, 7).unwrap().len(), 1);
        assert_eq!(interactions.angles(7, 0, 0).unwrap().len(), 1);

        interactions.add_angle(42, 42, 42, Harmonic{x0: 0.0, k: 0.0});
        assert_eq!(interactions.angles(42, 42, 42).unwrap().len(), 1);
    }

    #[test]
    fn dihedrals() {
        let mut interactions = Interactions::new();

        interactions.add_dihedral(0, 3, 7, 2, Harmonic{x0: 0.0, k: 0.0});
        assert_eq!(interactions.dihedrals(0, 3, 7, 2).unwrap().len(), 1);
        assert_eq!(interactions.dihedrals(2, 7, 3, 0).unwrap().len(), 1);

        interactions.add_dihedral(0, 0, 7, 2, Harmonic{x0: 0.0, k: 0.0});
        assert_eq!(interactions.dihedrals(0, 0, 7, 2).unwrap().len(), 1);
        assert_eq!(interactions.dihedrals(2, 7, 0, 0).unwrap().len(), 1);

        interactions.add_dihedral(0, 0, 9, 0, Harmonic{x0: 0.0, k: 0.0});
        assert_eq!(interactions.dihedrals(0, 0, 9, 0).unwrap().len(), 1);
        assert_eq!(interactions.dihedrals(0, 9, 0, 0).unwrap().len(), 1);

        interactions.add_dihedral(42, 42, 42, 42, Harmonic{x0: 0.0, k: 0.0});
        assert_eq!(interactions.dihedrals(42, 42, 42, 42).unwrap().len(), 1);
    }
}
