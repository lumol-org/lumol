// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

//! The `Interactions` type for storing particles types to potentials
//! associations.

use std::cmp::{max, min};
use std::collections::BTreeMap;
use std::f64;

use energy::{AnglePotential, BondPotential, DihedralPotential, PairInteraction};
use energy::{CoulombicPotential, GlobalPotential};
use sys::ParticleKind as Kind;

/// Normalize pair indexes to get a canonical representation
#[inline]
fn normalize_pair(i: Kind, j: Kind) -> (Kind, Kind) {
    if i < j {
        (i, j)
    } else {
        (j, i)
    }
}

/// Normalize angle indexes to get a canonical representation
#[inline]
fn normalize_angle(i: Kind, j: Kind, k: Kind) -> (Kind, Kind, Kind) {
    if i < k {
        (i, j, k)
    } else {
        (k, j, i)
    }
}

/// Normalize dihedral indexes to get a canonical representation
#[inline]
fn normalize_dihedral(i: Kind, j: Kind, k: Kind, m: Kind) -> (Kind, Kind, Kind, Kind) {
    let max_ij = max(i, j);
    let max_km = max(k, m);
    if max_ij == max_km {
        if min(i, j) < min(k, m) {
            (i, j, k, m)
        } else {
            (m, k, j, i)
        }
    } else if max_ij < max_km {
        (i, j, k, m)
    } else {
        (m, k, j, i)
    }
}

type PairKind = (Kind, Kind);
type BondKind = (Kind, Kind);
type AngleKind = (Kind, Kind, Kind);
type DihedralKind = (Kind, Kind, Kind, Kind);

/// The `Interaction` type hold all data about the potentials in the system.
///
/// Its main role is to store and provide access
#[derive(Clone)]
pub struct Interactions {
    /// Pair potentials
    pairs: BTreeMap<PairKind, Vec<PairInteraction>>,
    /// Bond potentials
    bonds: BTreeMap<BondKind, Vec<Box<BondPotential>>>,
    /// Angle potentials
    angles: BTreeMap<AngleKind, Vec<Box<AnglePotential>>>,
    /// Dihedral angles potentials
    dihedrals: BTreeMap<DihedralKind, Vec<Box<DihedralPotential>>>,
    /// Coulombic potential solver
    pub coulomb: Option<Box<CoulombicPotential>>,
    /// Global potentials
    pub globals: Vec<Box<GlobalPotential>>,
}

impl Interactions {
    /// Create a new empty `Interactions`
    pub fn new() -> Interactions {
        Interactions {
            pairs: BTreeMap::new(),
            bonds: BTreeMap::new(),
            angles: BTreeMap::new(),
            dihedrals: BTreeMap::new(),
            coulomb: None,
            globals: Vec::new(),
        }
    }

    /// Add the `potential` pair interaction for the pair `(i, j)`
    pub fn add_pair(&mut self, i: Kind, j: Kind, potential: PairInteraction) {
        let (i, j) = normalize_pair(i, j);
        let pairs = self.pairs.entry((i, j)).or_insert(Vec::new());
        pairs.push(potential);
    }

    /// Add the `potential` bonded interaction for the pair `(i, j)`
    pub fn add_bond(&mut self, i: Kind, j: Kind, potential: Box<BondPotential>) {
        let (i, j) = normalize_pair(i, j);
        let bonds = self.bonds.entry((i, j)).or_insert(Vec::new());
        bonds.push(potential);
    }

    /// Add the `potential` angle interaction for the angle `(i, j, k)`
    pub fn add_angle(&mut self, i: Kind, j: Kind, k: Kind, potential: Box<AnglePotential>) {
        let (i, j, k) = normalize_angle(i, j, k);
        let angles = self.angles.entry((i, j, k)).or_insert(Vec::new());
        angles.push(potential);
    }

    /// Add the `potential` dihedral interaction for the dihedral angle `(i, j,
    /// k, m)`
    pub fn add_dihedral(
        &mut self,
        i: Kind,
        j: Kind,
        k: Kind,
        m: Kind,
        potential: Box<DihedralPotential>,
    ) {
        let (i, j, k, m) = normalize_dihedral(i, j, k, m);
        let dihedrals = self.dihedrals.entry((i, j, k, m)).or_insert(Vec::new());
        dihedrals.push(potential);
    }
}

impl Interactions {
    /// Get all pair interactions corresponding to the pair `(i, j)`
    pub fn pairs(&self, i: Kind, j: Kind) -> &[PairInteraction] {
        let (i, j) = normalize_pair(i, j);
        self.pairs.get(&(i, j)).map_or(&[], |pairs| &**pairs)
    }

    /// Get all bonded interactions corresponding to the pair `(i, j)`
    pub fn bonds(&self, i: Kind, j: Kind) -> &[Box<BondPotential>] {
        let (i, j) = normalize_pair(i, j);
        self.bonds.get(&(i, j)).map_or(&[], |bonds| &**bonds)
    }

    /// Get all angle interactions corresponding to the angle `(i, j, k)`
    pub fn angles(&self, i: Kind, j: Kind, k: Kind) -> &[Box<AnglePotential>] {
        let (i, j, k) = normalize_angle(i, j, k);
        self.angles.get(&(i, j, k)).map_or(&[], |angles| &**angles)
    }

    /// Get all dihedral interactions corresponding to the dihedral `(i, j, k, m)`
    pub fn dihedrals(&self, i: Kind, j: Kind, k: Kind, m: Kind) -> &[Box<DihedralPotential>] {
        let (i, j, k, m) = normalize_dihedral(i, j, k, m);
        self.dihedrals.get(&(i, j, k, m)).map_or(&[], |dihedrals| &**dihedrals)
    }

    /// Get maximum cutoff from `coulomb`, `pairs` and `global` interactons.
    pub fn maximum_cutoff(&self) -> Option<f64> {
        // Coulomb potential, return cutoff
        let coulomb_cutoff = match self.coulomb {
            Some(ref coulomb) => {
                match coulomb.cutoff() {
                    Some(cutoff) => cutoff,
                    None => f64::NAN,
                }
            }
            None => f64::NAN,
        };

        // Go through global interactions, return maximum cutoff
        let global_cutoff = self.globals.iter()
                                .map(|i| i.cutoff())
                                .filter_map(|rc| rc)
                                .fold(f64::NAN, f64::max);

        let mut maximum_cutoff = f64::max(global_cutoff, coulomb_cutoff);

        // Pair interactions, return maximum cutoff
        let pairs_cutoff = self.pairs.values()
                               .flat_map(|i| i.iter().map(|pair| pair.cutoff()))
                               .fold(f64::NAN, f64::max);

        maximum_cutoff = f64::max(maximum_cutoff, pairs_cutoff);
        if maximum_cutoff.is_nan() {
            None
        } else {
            Some(maximum_cutoff)
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    use energy::{NullPotential, PairInteraction, Wolf};
    use sys::ParticleKind as Kind;

    #[test]
    fn normalizing_pairs() {
        assert_eq!(normalize_pair(Kind(0), Kind(1)), normalize_pair(Kind(1), Kind(0)));
        assert_eq!(normalize_pair(Kind(125), Kind(0)), normalize_pair(Kind(0), Kind(125)));

        assert_eq!(normalize_pair(Kind(0), Kind(1)), (Kind(0), Kind(1)));
        assert_eq!(normalize_pair(Kind(125), Kind(0)), (Kind(0), Kind(125)));
    }

    #[test]
    fn normalizing_angles() {
        // Checking that equivalent angles are normalized to the same tuple
        assert_eq!(
            normalize_angle(Kind(1), Kind(0), Kind(2)),
            normalize_angle(Kind(2), Kind(0), Kind(1))
        );
        assert_eq!(
            normalize_angle(Kind(15), Kind(4), Kind(8)),
            normalize_angle(Kind(8), Kind(4), Kind(15))
        );
        assert_eq!(
            normalize_angle(Kind(0), Kind(0), Kind(1)),
            normalize_angle(Kind(1), Kind(0), Kind(0))
        );
        assert_eq!(
            normalize_angle(Kind(10), Kind(10), Kind(1)),
            normalize_angle(Kind(1), Kind(10), Kind(10))
        );

        // Checking the actual normalized value
        assert_eq!(normalize_angle(Kind(1), Kind(0), Kind(2)), (Kind(1), Kind(0), Kind(2)));
        assert_eq!(normalize_angle(Kind(15), Kind(4), Kind(8)), (Kind(8), Kind(4), Kind(15)));
        assert_eq!(normalize_angle(Kind(0), Kind(0), Kind(1)), (Kind(0), Kind(0), Kind(1)));
        assert_eq!(normalize_angle(Kind(10), Kind(10), Kind(1)), (Kind(1), Kind(10), Kind(10)));
    }

    #[test]
    fn normalizing_dihedrals() {
        // Checking that equivalent dihedrals are normalized to the same tuple
        assert_eq!(
            normalize_dihedral(Kind(1), Kind(0), Kind(2), Kind(3)),
            normalize_dihedral(Kind(3), Kind(2), Kind(0), Kind(1))
        );
        assert_eq!(
            normalize_dihedral(Kind(10), Kind(8), Kind(2), Kind(3)),
            normalize_dihedral(Kind(3), Kind(2), Kind(8), Kind(10))
        );
        assert_eq!(
            normalize_dihedral(Kind(0), Kind(0), Kind(2), Kind(3)),
            normalize_dihedral(Kind(3), Kind(2), Kind(0), Kind(0))
        );
        assert_eq!(
            normalize_dihedral(Kind(1), Kind(5), Kind(0), Kind(0)),
            normalize_dihedral(Kind(0), Kind(0), Kind(5), Kind(1))
        );
        assert_eq!(
            normalize_dihedral(Kind(10), Kind(10), Kind(2), Kind(3)),
            normalize_dihedral(Kind(3), Kind(2), Kind(10), Kind(10))
        );
        assert_eq!(
            normalize_dihedral(Kind(1), Kind(5), Kind(10), Kind(10)),
            normalize_dihedral(Kind(10), Kind(10), Kind(5), Kind(1))
        );
        assert_eq!(
            normalize_dihedral(Kind(0), Kind(0), Kind(0), Kind(3)),
            normalize_dihedral(Kind(3), Kind(0), Kind(0), Kind(0))
        );
        assert_eq!(
            normalize_dihedral(Kind(1), Kind(5), Kind(5), Kind(5)),
            normalize_dihedral(Kind(5), Kind(5), Kind(5), Kind(1))
        );
        assert_eq!(
            normalize_dihedral(Kind(0), Kind(0), Kind(3), Kind(0)),
            normalize_dihedral(Kind(0), Kind(3), Kind(0), Kind(0))
        );
        assert_eq!(
            normalize_dihedral(Kind(5), Kind(1), Kind(5), Kind(5)),
            normalize_dihedral(Kind(5), Kind(5), Kind(1), Kind(5))
        );

        // Checking the actual normalized value
        assert_eq!(
            normalize_dihedral(Kind(1), Kind(0), Kind(2), Kind(3)),
            (Kind(1), Kind(0), Kind(2), Kind(3))
        );
        assert_eq!(
            normalize_dihedral(Kind(10), Kind(8), Kind(2), Kind(3)),
            (Kind(3), Kind(2), Kind(8), Kind(10))
        );
        assert_eq!(
            normalize_dihedral(Kind(0), Kind(0), Kind(2), Kind(3)),
            (Kind(0), Kind(0), Kind(2), Kind(3))
        );
        assert_eq!(
            normalize_dihedral(Kind(1), Kind(5), Kind(0), Kind(0)),
            (Kind(0), Kind(0), Kind(5), Kind(1))
        );
        assert_eq!(
            normalize_dihedral(Kind(10), Kind(10), Kind(2), Kind(3)),
            (Kind(3), Kind(2), Kind(10), Kind(10))
        );
        assert_eq!(
            normalize_dihedral(Kind(1), Kind(5), Kind(10), Kind(10)),
            (Kind(1), Kind(5), Kind(10), Kind(10))
        );
        assert_eq!(
            normalize_dihedral(Kind(0), Kind(0), Kind(0), Kind(3)),
            (Kind(0), Kind(0), Kind(0), Kind(3))
        );
        assert_eq!(
            normalize_dihedral(Kind(1), Kind(5), Kind(5), Kind(5)),
            (Kind(1), Kind(5), Kind(5), Kind(5))
        );
        assert_eq!(
            normalize_dihedral(Kind(0), Kind(0), Kind(3), Kind(0)),
            (Kind(0), Kind(0), Kind(3), Kind(0))
        );
        assert_eq!(
            normalize_dihedral(Kind(5), Kind(1), Kind(5), Kind(5)),
            (Kind(5), Kind(1), Kind(5), Kind(5))
        );
    }

    #[test]
    fn pairs() {
        let mut interactions = Interactions::new();
        let pair = PairInteraction::new(Box::new(NullPotential), 0.0);
        interactions.add_pair(Kind(0), Kind(1), pair.clone());
        assert_eq!(interactions.pairs(Kind(0), Kind(1)).len(), 1);
        assert_eq!(interactions.pairs(Kind(1), Kind(0)).len(), 1);

        assert_eq!(interactions.pairs(Kind(0), Kind(0)).len(), 0);
        interactions.add_pair(Kind(0), Kind(0), pair.clone());
        assert_eq!(interactions.pairs(Kind(0), Kind(0)).len(), 1);

        // 'out of bounds' kinds
        assert_eq!(interactions.pairs(Kind(55), Kind(55)).len(), 0);
    }

    #[test]
    fn bonds() {
        let mut interactions = Interactions::new();

        interactions.add_bond(Kind(0), Kind(1), Box::new(NullPotential));
        assert_eq!(interactions.bonds(Kind(0), Kind(1)).len(), 1);
        assert_eq!(interactions.bonds(Kind(1), Kind(0)).len(), 1);

        assert_eq!(interactions.bonds(Kind(0), Kind(0)).len(), 0);
        interactions.add_bond(Kind(0), Kind(0), Box::new(NullPotential));
        assert_eq!(interactions.bonds(Kind(0), Kind(0)).len(), 1);


        // 'out of bounds' kinds
        assert_eq!(interactions.bonds(Kind(55), Kind(55)).len(), 0);
    }

    #[test]
    fn angles() {
        let mut interactions = Interactions::new();

        interactions.add_angle(Kind(0), Kind(1), Kind(2), Box::new(NullPotential));
        assert_eq!(interactions.angles(Kind(0), Kind(1), Kind(2)).len(), 1);
        assert_eq!(interactions.angles(Kind(2), Kind(1), Kind(0)).len(), 1);

        assert_eq!(interactions.angles(Kind(2), Kind(2), Kind(1)).len(), 0);
        interactions.add_angle(Kind(2), Kind(2), Kind(1), Box::new(NullPotential));
        assert_eq!(interactions.angles(Kind(2), Kind(2), Kind(1)).len(), 1);
        assert_eq!(interactions.angles(Kind(1), Kind(2), Kind(2)).len(), 1);

        interactions.add_angle(Kind(3), Kind(3), Kind(3), Box::new(NullPotential));
        assert_eq!(interactions.angles(Kind(3), Kind(3), Kind(3)).len(), 1);

        // 'out of bounds' kinds
        assert_eq!(interactions.angles(Kind(55), Kind(55), Kind(55)).len(), 0);
    }

    #[test]
    fn dihedrals() {
        let mut interactions = Interactions::new();

        interactions.add_dihedral(Kind(0), Kind(1), Kind(2), Kind(3), Box::new(NullPotential));
        assert_eq!(interactions.dihedrals(Kind(0), Kind(1), Kind(2), Kind(3)).len(), 1);
        assert_eq!(interactions.dihedrals(Kind(3), Kind(2), Kind(1), Kind(0)).len(), 1);

        assert_eq!(interactions.dihedrals(Kind(0), Kind(0), Kind(1), Kind(2)).len(), 0);
        interactions.add_dihedral(Kind(0), Kind(0), Kind(1), Kind(2), Box::new(NullPotential));
        assert_eq!(interactions.dihedrals(Kind(0), Kind(0), Kind(1), Kind(2)).len(), 1);
        assert_eq!(interactions.dihedrals(Kind(2), Kind(1), Kind(0), Kind(0)).len(), 1);

        interactions.add_dihedral(Kind(0), Kind(0), Kind(1), Kind(0), Box::new(NullPotential));
        assert_eq!(interactions.dihedrals(Kind(0), Kind(0), Kind(1), Kind(0)).len(), 1);
        assert_eq!(interactions.dihedrals(Kind(0), Kind(1), Kind(0), Kind(0)).len(), 1);

        interactions.add_dihedral(Kind(4), Kind(4), Kind(4), Kind(4), Box::new(NullPotential));
        assert_eq!(interactions.dihedrals(Kind(4), Kind(4), Kind(4), Kind(4)).len(), 1);

        // 'out of bounds' kinds
        assert_eq!(interactions.dihedrals(Kind(55), Kind(55), Kind(55), Kind(55)).len(), 0);
    }

    #[test]
    fn test_maximum_cutoff() {
        let mut interactions = Interactions::new();
        // No cutoff
        assert_eq!(interactions.maximum_cutoff(), None);

        // Only pairs
        let pair = PairInteraction::new(Box::new(NullPotential), 10.0);
        interactions.add_pair(Kind(0), Kind(1), pair.clone());
        assert_eq!(interactions.maximum_cutoff(), Some(10.0));

        // Only globals
        let mut interactions = Interactions::new();
        interactions.globals.push(Box::new(Wolf::new(1.0)));
        assert_eq!(interactions.maximum_cutoff(), Some(1.0));

        interactions.globals.push(Box::new(Wolf::new(5.0)));
        assert_eq!(interactions.maximum_cutoff(), Some(5.0));

        // Only coulomb
        let mut interactions = Interactions::new();
        interactions.coulomb = Some(Box::new(Wolf::new(1.0)));
        assert_eq!(interactions.maximum_cutoff(), Some(1.0));

        // All potentials
        interactions.add_pair(Kind(0), Kind(1), pair.clone());
        assert_eq!(interactions.maximum_cutoff(), Some(10.0));
        interactions.globals.push(Box::new(Wolf::new(15.0)));
        assert_eq!(interactions.maximum_cutoff(), Some(15.0));
        interactions.globals.push(Box::new(Wolf::new(1.0)));
        assert_eq!(interactions.maximum_cutoff(), Some(15.0));
    }
}
