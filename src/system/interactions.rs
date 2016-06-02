// Cymbalum, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

//! The `Interactions` type for storing particles types to potentials
//! associations.

use std::collections::BTreeMap;
use std::cmp::max;
use std::cell::RefCell;

use potentials::{PairPotential, AnglePotential, DihedralPotential};
use potentials::{GlobalPotential, CoulombicPotential};
use potentials::PairRestriction;

use system::ParticleKind as Kind;

/// Type associating a potential and a pair restriction
pub type PairInteraction = (Box<PairPotential>, PairRestriction);

/// The Interaction type hold all data about the potentials in the system,
/// indexed by particle type.
#[derive(Clone)]
pub struct Interactions {
    /// Pair potentials
    pairs: BTreeMap<(Kind, Kind), Vec<PairInteraction>>,
    /// Bond potentials
    bonds: BTreeMap<(Kind, Kind), Vec<Box<PairPotential>>>,
    /// Angle potentials
    angles: BTreeMap<(Kind, Kind, Kind), Vec<Box<AnglePotential>>>,
    /// Dihedral angles potentials
    dihedrals: BTreeMap<(Kind, Kind, Kind, Kind), Vec<Box<DihedralPotential>>>,
    /// Coulombic potential solver
    coulomb: Option<RefCell<Box<CoulombicPotential>>>,
    /// Global potentials
    globals: Vec<RefCell<Box<GlobalPotential>>>,
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
    pub fn add_pair(&mut self, i: Kind, j:Kind, potential: Box<PairPotential>) {
        self.add_pair_with_restriction(i, j, potential, PairRestriction::None);
    }

    /// Add the `potential` pair interaction to the pair `(i, j)`, with the
    /// restriction scheme `restrict`.
    pub fn add_pair_with_restriction(&mut self, i: Kind, j:Kind, potential: Box<PairPotential>, restrict: PairRestriction) {
        let (i, j) = sort_pair(i, j);
        let pairs = self.pairs.entry((i, j)).or_insert(Vec::new());
        pairs.push((potential, restrict));
    }

    /// Get all pair interactions corresponding to the pair `(i, j)`
    pub fn pairs(&self, i: Kind, j:Kind) -> Option<&Vec<PairInteraction>> {
        let (i, j) = sort_pair(i, j);
        self.pairs.get(&(i, j))
    }

    /// Add the `potential` bonded interaction to the pair `(i, j)`
    pub fn add_bond(&mut self, i: Kind, j:Kind, potential: Box<PairPotential>) {
        let (i, j) = sort_pair(i, j);
        let bonds = self.bonds.entry((i, j)).or_insert(Vec::new());
        bonds.push(potential);
    }

    /// Get all bonded interactions corresponding to the pair `(i, j)`
    pub fn bonds(&self, i: Kind, j:Kind) -> Option<&Vec<Box<PairPotential>>> {
        let (i, j) = sort_pair(i, j);
        self.bonds.get(&(i, j))
    }

    /// Add the `potential` angle interaction to the angle `(i, j, k)`
    pub fn add_angle(&mut self, i: Kind, j:Kind, k:Kind, potential: Box<AnglePotential>) {
        let (i, j, k) = sort_angle(i, j, k);
        let angles = self.angles.entry((i, j, k)).or_insert(Vec::new());
        angles.push(potential);
    }

    /// Get all angle interactions corresponding to the angle `(i, j, k)`
    pub fn angles(&self, i: Kind, j:Kind, k:Kind) -> Option<&Vec<Box<AnglePotential>>> {
        let (i, j, k) = sort_angle(i, j, k);
        self.angles.get(&(i, j, k))
    }

    /// Add the `potential` dihedral interaction to the dihedral `(i, j, k, m)`
    pub fn add_dihedral(&mut self, i: Kind, j:Kind, k:Kind, m:Kind, potential: Box<DihedralPotential>) {
        let (i, j, k, m) = sort_dihedral(i, j, k, m);
        let dihedrals = self.dihedrals.entry((i, j, k, m)).or_insert(Vec::new());
        dihedrals.push(potential);
    }

    /// Get all dihedral interactions corresponding to the dihedral `(i, j, k, m)`
    pub fn dihedrals(&self, i: Kind, j:Kind, k:Kind, m:Kind) -> Option<&Vec<Box<DihedralPotential>>> {
        let (i, j, k, m) = sort_dihedral(i, j, k, m);
        self.dihedrals.get(&(i, j, k, m))
    }

    /// Set the coulombic interaction for all pairs to `potential`
    pub fn set_coulomb(&mut self, potential: Box<CoulombicPotential>) {
        self.coulomb = Some(RefCell::new(potential));
    }

    /// Get the coulombic interaction as a RefCell, because the
    /// `GlobalPotential` are allowed to mutate themself when computing energy.
    pub fn coulomb(&self) -> Option<&RefCell<Box<CoulombicPotential>>> {
        self.coulomb.as_ref()
    }

    /// Add the `potential` global interaction
    pub fn add_global(&mut self, potential: Box<GlobalPotential>) {
        self.globals.push(RefCell::new(potential));
    }

    /// Get all global interactions
    pub fn globals(&self) -> &[RefCell<Box<GlobalPotential>>] {
        &self.globals
    }
}

/// Sort pair indexes to get a cannonical representation
#[inline] fn sort_pair(i: Kind, j:Kind) -> (Kind, Kind) {
    if i < j {
        (i, j)
    } else {
        (j, i)
    }
}

/// Sort angle indexes to get a cannonical representation
#[inline] fn sort_angle(i: Kind, j:Kind, k:Kind) -> (Kind, Kind, Kind) {
    if i < k {
        (i, j, k)
    } else {
        (k, j, i)
    }
}

/// Sort dihedral indexes to get a cannonical representation
#[inline] fn sort_dihedral(i: Kind, j:Kind, k:Kind, m:Kind) -> (Kind, Kind, Kind, Kind) {
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
    use system::ParticleKind as Kind;

    #[test]
    fn pairs() {
        let mut interactions = Interactions::new();

        interactions.add_pair(Kind(0), Kind(3), Box::new(Harmonic{x0: 0.0, k: 0.0}));
        assert_eq!(interactions.pairs(Kind(0), Kind(3)).unwrap().len(), 1);
        assert_eq!(interactions.pairs(Kind(3), Kind(0)).unwrap().len(), 1);

        interactions.add_pair(Kind(0), Kind(0), Box::new(Harmonic{x0: 0.0, k: 0.0}));
        assert_eq!(interactions.pairs(Kind(0), Kind(0)).unwrap().len(), 1);
    }

    #[test]
    fn bonds() {
        let mut interactions = Interactions::new();

        interactions.add_bond(Kind(0), Kind(3), Box::new(Harmonic{x0: 0.0, k: 0.0}));
        assert_eq!(interactions.bonds(Kind(0), Kind(3)).unwrap().len(), 1);
        assert_eq!(interactions.bonds(Kind(3), Kind(0)).unwrap().len(), 1);

        interactions.add_bond(Kind(0), Kind(0), Box::new(Harmonic{x0: 0.0, k: 0.0}));
        assert_eq!(interactions.bonds(Kind(0), Kind(0)).unwrap().len(), 1);
    }

    #[test]
    fn angles() {
        let mut interactions = Interactions::new();

        interactions.add_angle(Kind(0), Kind(3), Kind(7), Box::new(Harmonic{x0: 0.0, k: 0.0}));
        assert_eq!(interactions.angles(Kind(0), Kind(3), Kind(7)).unwrap().len(), 1);
        assert_eq!(interactions.angles(Kind(7), Kind(3), Kind(0)).unwrap().len(), 1);

        interactions.add_angle(Kind(0), Kind(0), Kind(7), Box::new(Harmonic{x0: 0.0, k: 0.0}));
        assert_eq!(interactions.angles(Kind(0), Kind(0), Kind(7)).unwrap().len(), 1);
        assert_eq!(interactions.angles(Kind(7), Kind(0), Kind(0)).unwrap().len(), 1);

        interactions.add_angle(Kind(42), Kind(42), Kind(42), Box::new(Harmonic{x0: 0.0, k: 0.0}));
        assert_eq!(interactions.angles(Kind(42), Kind(42), Kind(42)).unwrap().len(), 1);
    }

    #[test]
    fn dihedrals() {
        let mut interactions = Interactions::new();

        interactions.add_dihedral(
            Kind(0), Kind(3), Kind(7), Kind(2),
            Box::new(Harmonic{x0: 0.0, k: 0.0})
        );
        assert_eq!(interactions.dihedrals(
            Kind(0), Kind(3), Kind(7), Kind(2)
        ).unwrap().len(), 1);
        assert_eq!(interactions.dihedrals(
            Kind(2), Kind(7), Kind(3), Kind(0)
        ).unwrap().len(), 1);

        interactions.add_dihedral(
            Kind(0), Kind(0), Kind(7), Kind(2),
            Box::new(Harmonic{x0: 0.0, k: 0.0})
        );
        assert_eq!(interactions.dihedrals(
            Kind(0), Kind(0), Kind(7), Kind(2)
        ).unwrap().len(), 1);
        assert_eq!(interactions.dihedrals(
            Kind(2), Kind(7), Kind(0), Kind(0)
        ).unwrap().len(), 1);

        interactions.add_dihedral(
            Kind(0), Kind(0), Kind(9), Kind(0),
            Box::new(Harmonic{x0: 0.0, k: 0.0}));
        assert_eq!(interactions.dihedrals(
            Kind(0), Kind(0), Kind(9), Kind(0)
        ).unwrap().len(), 1);
        assert_eq!(interactions.dihedrals(
            Kind(0), Kind(9), Kind(0), Kind(0)
        ).unwrap().len(), 1);

        interactions.add_dihedral(
            Kind(42), Kind(42), Kind(42), Kind(42),
            Box::new(Harmonic{x0: 0.0, k: 0.0})
        );
        assert_eq!(interactions.dihedrals(
            Kind(42), Kind(42), Kind(42), Kind(42)
        ).unwrap().len(), 1);
    }

    #[test]
    fn globals() {
        let mut interactions = Interactions::new();
        interactions.add_global(Box::new(Wolf::new(1.0)));
        assert_eq!(interactions.globals().len(), 1);
    }
}
