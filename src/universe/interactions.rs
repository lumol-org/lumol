/* Cymbalum, Molecular Simulation in Rust - Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
 */

//! The `Interactions` type for storing particles types to potentials
//! associations.

use std::collections::BTreeMap;
use std::cmp::max;

use potentials::{PairPotential, AnglePotential, DihedralPotential};
use potentials::{GlobalPotential, CoulombicPotential};
use potentials::PairRestriction;

/// Type associating a potential and a pair restriction
pub type PairInteraction = (Box<PairPotential>, PairRestriction);

/// The Interaction type hold all data about the potentials in the system,
/// indexed by particle type.
#[derive(Clone)]
pub struct Interactions {
    /// Pair potentials
    pairs: BTreeMap<(u16, u16), Vec<PairInteraction>>,
    /// Bond potentials
    bonds: BTreeMap<(u16, u16), Vec<Box<PairPotential>>>,
    /// Angle potentials
    angles: BTreeMap<(u16, u16, u16), Vec<Box<AnglePotential>>>,
    /// Dihedral angles potentials
    dihedrals: BTreeMap<(u16, u16, u16, u16), Vec<Box<DihedralPotential>>>,
    /// Coulombic potential solver
    coulomb: Option<Box<CoulombicPotential>>,
    /// Global potentials
    globals: Vec<Box<GlobalPotential>>,
}

impl Interactions {
    pub fn new() -> Interactions {
        Interactions{
            pairs: BTreeMap::new(),
            bonds: BTreeMap::new(),
            angles: BTreeMap::new(),
            dihedrals: BTreeMap::new(),
            coulomb: None,
            globals: Vec::new(),
        }
    }

    /// Add the `potential` pair interaction to the pair `(i, j)`
    pub fn add_pair(&mut self, i: u16, j:u16, potential: Box<PairPotential>) {
        self.add_pair_with_restriction(i, j, potential, PairRestriction::None);
    }

    /// Add the `potential` pair interaction to the pair `(i, j)`, with the
    /// restriction scheme `restrict`.
    pub fn add_pair_with_restriction(&mut self, i: u16, j:u16, potential: Box<PairPotential>, restrict: PairRestriction) {
        let (i, j) = sort_pair(i, j);
        let pairs = self.pairs.entry((i, j)).or_insert(Vec::new());
        pairs.push((potential, restrict));
    }

    /// Get all pair interactions corresponding to the pair `(i, j)`
    pub fn pairs(&self, i: u16, j:u16) -> Option<&Vec<PairInteraction>> {
        let (i, j) = sort_pair(i, j);
        self.pairs.get(&(i, j))
    }

    /// Add the `potential` bonded interaction to the pair `(i, j)`
    pub fn add_bond(&mut self, i: u16, j:u16, potential: Box<PairPotential>) {
        let (i, j) = sort_pair(i, j);
        let bonds = self.bonds.entry((i, j)).or_insert(Vec::new());
        bonds.push(potential);
    }

    /// Get all bonded interactions corresponding to the pair `(i, j)`
    pub fn bonds(&self, i: u16, j:u16) -> Option<&Vec<Box<PairPotential>>> {
        let (i, j) = sort_pair(i, j);
        self.bonds.get(&(i, j))
    }

    /// Add the `potential` angle interaction to the angle `(i, j, k)`
    pub fn add_angle(&mut self, i: u16, j:u16, k:u16, potential: Box<AnglePotential>) {
        let (i, j, k) = sort_angle(i, j, k);
        let angles = self.angles.entry((i, j, k)).or_insert(Vec::new());
        angles.push(potential);
    }

    /// Get all angle interactions corresponding to the angle `(i, j, k)`
    pub fn angles(&self, i: u16, j:u16, k:u16) -> Option<&Vec<Box<AnglePotential>>> {
        let (i, j, k) = sort_angle(i, j, k);
        self.angles.get(&(i, j, k))
    }

    /// Add the `potential` dihedral interaction to the dihedral `(i, j, k, m)`
    pub fn add_dihedral(&mut self, i: u16, j:u16, k:u16, m:u16, potential: Box<DihedralPotential>) {
        let (i, j, k, m) = sort_dihedral(i, j, k, m);
        let dihedrals = self.dihedrals.entry((i, j, k, m)).or_insert(Vec::new());
        dihedrals.push(potential);
    }

    /// Get all dihedral interactions corresponding to the dihedral `(i, j, k, m)`
    pub fn dihedrals(&self, i: u16, j:u16, k:u16, m:u16) -> Option<&Vec<Box<DihedralPotential>>> {
        let (i, j, k, m) = sort_dihedral(i, j, k, m);
        self.dihedrals.get(&(i, j, k, m))
    }

    /// Set the coulombic interaction for all pairs to `potential`
    pub fn set_coulomb(&mut self, potential: Box<CoulombicPotential>) {
        self.coulomb = Some(potential);
    }

    /// Get the coulombic interaction
    pub fn coulomb(&self) -> Option<&Box<CoulombicPotential>> {
        self.coulomb.as_ref()
    }

    /// Get the coulombic interaction as a mutable reference
    pub fn coulomb_mut(&mut self) -> Option<&mut Box<CoulombicPotential>> {
        self.coulomb.as_mut()
    }

    /// Add the `potential` global interaction
    pub fn add_global(&mut self, potential: Box<GlobalPotential>) {
        self.globals.push(potential);
    }

    /// Get all global interactions
    pub fn globals(&self) -> &[Box<GlobalPotential>] {
        &self.globals
    }

    /// Get all global interactions as mutable references
    pub fn globals_mut(&mut self) -> &mut[Box<GlobalPotential>] {
        &mut self.globals
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
    use potentials::{Harmonic, Wolf};

    #[test]
    fn pairs() {
        let mut interactions = Interactions::new();

        interactions.add_pair(0, 3, Box::new(Harmonic{x0: 0.0, k: 0.0}));
        assert_eq!(interactions.pairs(0, 3).unwrap().len(), 1);
        assert_eq!(interactions.pairs(3, 0).unwrap().len(), 1);

        interactions.add_pair(0, 0, Box::new(Harmonic{x0: 0.0, k: 0.0}));
        assert_eq!(interactions.pairs(0, 0).unwrap().len(), 1);
    }

    #[test]
    fn bonds() {
        let mut interactions = Interactions::new();

        interactions.add_bond(0, 3, Box::new(Harmonic{x0: 0.0, k: 0.0}));
        assert_eq!(interactions.bonds(0, 3).unwrap().len(), 1);
        assert_eq!(interactions.bonds(3, 0).unwrap().len(), 1);

        interactions.add_bond(0, 0, Box::new(Harmonic{x0: 0.0, k: 0.0}));
        assert_eq!(interactions.bonds(0, 0).unwrap().len(), 1);
    }

    #[test]
    fn angles() {
        let mut interactions = Interactions::new();

        interactions.add_angle(0, 3, 7, Box::new(Harmonic{x0: 0.0, k: 0.0}));
        assert_eq!(interactions.angles(0, 3, 7).unwrap().len(), 1);
        assert_eq!(interactions.angles(7, 3, 0).unwrap().len(), 1);

        interactions.add_angle(0, 0, 7, Box::new(Harmonic{x0: 0.0, k: 0.0}));
        assert_eq!(interactions.angles(0, 0, 7).unwrap().len(), 1);
        assert_eq!(interactions.angles(7, 0, 0).unwrap().len(), 1);

        interactions.add_angle(42, 42, 42, Box::new(Harmonic{x0: 0.0, k: 0.0}));
        assert_eq!(interactions.angles(42, 42, 42).unwrap().len(), 1);
    }

    #[test]
    fn dihedrals() {
        let mut interactions = Interactions::new();

        interactions.add_dihedral(0, 3, 7, 2, Box::new(Harmonic{x0: 0.0, k: 0.0}));
        assert_eq!(interactions.dihedrals(0, 3, 7, 2).unwrap().len(), 1);
        assert_eq!(interactions.dihedrals(2, 7, 3, 0).unwrap().len(), 1);

        interactions.add_dihedral(0, 0, 7, 2, Box::new(Harmonic{x0: 0.0, k: 0.0}));
        assert_eq!(interactions.dihedrals(0, 0, 7, 2).unwrap().len(), 1);
        assert_eq!(interactions.dihedrals(2, 7, 0, 0).unwrap().len(), 1);

        interactions.add_dihedral(0, 0, 9, 0, Box::new(Harmonic{x0: 0.0, k: 0.0}));
        assert_eq!(interactions.dihedrals(0, 0, 9, 0).unwrap().len(), 1);
        assert_eq!(interactions.dihedrals(0, 9, 0, 0).unwrap().len(), 1);

        interactions.add_dihedral(42, 42, 42, 42, Box::new(Harmonic{x0: 0.0, k: 0.0}));
        assert_eq!(interactions.dihedrals(42, 42, 42, 42).unwrap().len(), 1);
    }

    #[test]
    fn globals() {
        let mut interactions = Interactions::new();
        interactions.add_global(Box::new(Wolf::new(1.0)));
        assert_eq!(interactions.globals().len(), 1);
    }
}
