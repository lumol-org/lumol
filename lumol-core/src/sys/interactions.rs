// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

//! The `Interactions` type for storing particles types to potentials
//! associations.

use std::cmp::{max, min};
use std::collections::BTreeMap;
use std::f64;

use log::warn;

use crate::{AnglePotential, BondPotential, DihedralPotential, PairInteraction};
use crate::{CoulombicPotential, GlobalPotential};
use crate::ParticleKind;

pub type PairKind = (ParticleKind, ParticleKind);
pub type BondKind = (ParticleKind, ParticleKind);
pub type AngleKind = (ParticleKind, ParticleKind, ParticleKind);
pub type DihedralKind = (ParticleKind, ParticleKind, ParticleKind, ParticleKind);

/// Normalize pair indexes to get a canonical representation
#[inline]
fn normalize_pair((i, j): PairKind) -> PairKind {
    if i < j {
        (i, j)
    } else {
        (j, i)
    }
}

/// Normalize angle indexes to get a canonical representation
#[inline]
fn normalize_angle((i, j, k): AngleKind) -> AngleKind {
    if i < k {
        (i, j, k)
    } else {
        (k, j, i)
    }
}

/// Normalize dihedral indexes to get a canonical representation
#[inline]
fn normalize_dihedral((i, j, k, m): DihedralKind) -> DihedralKind {
    match (max(i, j), max(k, m)) {
        (ij, km) if ij == km => {
            if min(i, j) < min(k, m) {
                (i, j, k, m)
            } else {
                (m, k, j, i)
            }
        },
        (ij, km) if ij < km => (i, j, k, m),
        (_, _) => (m, k, j, i),
    }
}

/// The `Interaction` type hold all data about the potentials in the system.
///
/// Its main role is to store and provide access
#[derive(Clone)]
pub struct Interactions {
    /// Coulombic potential solver
    pub coulomb: Option<Box<dyn CoulombicPotential>>,
    /// Global potentials
    pub globals: Vec<Box<dyn GlobalPotential>>,
    /// Pair potentials
    pairs: BTreeMap<PairKind, PairInteraction>,
    /// Bond potentials
    bonds: BTreeMap<BondKind, Box<dyn BondPotential>>,
    /// Angle potentials
    angles: BTreeMap<AngleKind, Box<dyn AnglePotential>>,
    /// Dihedral angles potentials
    dihedrals: BTreeMap<DihedralKind, Box<dyn DihedralPotential>>,
    /// Association particles names to particle kinds
    kinds: BTreeMap<String, ParticleKind>,
}

impl Interactions {
    /// Create a new empty `Interactions`
    pub fn new() -> Interactions {
        Interactions {
            coulomb: None,
            globals: Vec::new(),
            pairs: BTreeMap::new(),
            bonds: BTreeMap::new(),
            angles: BTreeMap::new(),
            dihedrals: BTreeMap::new(),
            kinds: BTreeMap::new(),
        }
    }

    /// Get the existing kind associated with `name` or create a new one
    pub(crate) fn get_kind(&mut self, name: &str) -> ParticleKind {
        if let Some(&kind) = self.kinds.get(name) {
            kind
        } else {
            let kind = ParticleKind(self.kinds.len() as u32);
            let _ = self.kinds.insert(String::from(name), kind);
            kind
        }
    }

    /// Set the pair interaction `potential` for atoms with types `i` and `j`
    pub fn set_pair(&mut self, (i, j): (&str, &str), potential: PairInteraction) {
        let kind = normalize_pair((self.get_kind(i), self.get_kind(j)));
        if self.pairs.insert(kind, potential).is_some() {
            warn!("replaced pair potential for ({}, {})", i, j);
        }
    }

    /// Set the bond interaction `potential` for atoms with types `i` and `j`
    pub fn set_bond(&mut self, (i, j): (&str, &str), potential: Box<dyn BondPotential>) {
        let kind = normalize_pair((self.get_kind(i), self.get_kind(j)));
        if self.bonds.insert(kind, potential).is_some() {
            warn!("replaced bond potential for ({}, {})", i, j);
        }
    }

    /// Set the angle interaction `potential` for atoms with types `i`, `j`, and `k`
    pub fn set_angle(&mut self, (i, j, k): (&str, &str, &str), potential: Box<dyn AnglePotential>) {
        let kind = normalize_angle((self.get_kind(i), self.get_kind(j), self.get_kind(k)));
        if self.angles.insert(kind, potential).is_some() {
            warn!("replaced angle potential for ({}, {}, {})", i, j, k);
        }
    }

    /// Set the dihedral angle interaction `potential` for atoms with types
    /// `i`, `j`, `k`, and `m`.
    pub fn set_dihedral(&mut self, (i, j, k, m): (&str, &str, &str, &str), potential: Box<dyn DihedralPotential>) {
        let kind = (self.get_kind(i), self.get_kind(j), self.get_kind(k), self.get_kind(m));
        let kind = normalize_dihedral(kind);
        if self.dihedrals.insert(kind, potential).is_some() {
            warn!("replaced dihedral angle potential for ({}, {}, {}, {})", i, j, k, m);
        }
    }
}


impl Interactions {
    /// Get the pair interactions corresponding to the `pair`, if any exists.
    pub fn pair(&self, pair: PairKind) -> Option<&PairInteraction> {
        let kind = normalize_pair(pair);
        self.pairs.get(&kind)
    }

    /// Get the bond interactions corresponding to the `bond`, if any exists.
    pub fn bond(&self, bond: BondKind) -> Option<&dyn BondPotential> {
        let kind = normalize_pair(bond);
        self.bonds.get(&kind).map(|potential| &**potential)
    }

    /// Get the angle interactions corresponding to the `angle`, if any exists.
    pub fn angle(&self, angle: AngleKind) -> Option<&dyn AnglePotential> {
        let kind = normalize_angle(angle);
        self.angles.get(&kind).map(|potential| &**potential)
    }

    /// Get the dihedral interactions corresponding to the `dihedral`, if any exists.
    pub fn dihedral(&self, dihedral: DihedralKind) -> Option<&dyn DihedralPotential> {
        let kind = normalize_dihedral(dihedral);
        self.dihedrals.get(&kind).map(|potential| &**potential)
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
                                .filter_map(|i| i.cutoff())
                                .fold(f64::NAN, f64::max);

        let mut maximum_cutoff = f64::max(global_cutoff, coulomb_cutoff);

        // Pair interactions, return maximum cutoff
        let pairs_cutoff = self.pairs.values()
                               .map(|pair| pair.cutoff())
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

    use crate::{NullPotential, PairInteraction, Wolf};
    use crate::ParticleKind as Kind;

    #[test]
    fn normalizing_pairs() {
        assert_eq!(normalize_pair((Kind(0), Kind(1))), normalize_pair((Kind(1), Kind(0))));
        assert_eq!(normalize_pair((Kind(125), Kind(0))), normalize_pair((Kind(0), Kind(125))));

        assert_eq!(normalize_pair((Kind(0), Kind(1))), (Kind(0), Kind(1)));
        assert_eq!(normalize_pair((Kind(125), Kind(0))), (Kind(0), Kind(125)));
    }

    #[test]
    fn normalizing_angles() {
        // Checking that equivalent angles are normalized to the same tuple
        assert_eq!(
            normalize_angle((Kind(1), Kind(0), Kind(2))),
            normalize_angle((Kind(2), Kind(0), Kind(1)))
        );
        assert_eq!(
            normalize_angle((Kind(15), Kind(4), Kind(8))),
            normalize_angle((Kind(8), Kind(4), Kind(15)))
        );
        assert_eq!(
            normalize_angle((Kind(0), Kind(0), Kind(1))),
            normalize_angle((Kind(1), Kind(0), Kind(0)))
        );
        assert_eq!(
            normalize_angle((Kind(10), Kind(10), Kind(1))),
            normalize_angle((Kind(1), Kind(10), Kind(10)))
        );

        // Checking the actual normalized value
        assert_eq!(normalize_angle((Kind(1), Kind(0), Kind(2))), (Kind(1), Kind(0), Kind(2)));
        assert_eq!(normalize_angle((Kind(15), Kind(4), Kind(8))), (Kind(8), Kind(4), Kind(15)));
        assert_eq!(normalize_angle((Kind(0), Kind(0), Kind(1))), (Kind(0), Kind(0), Kind(1)));
        assert_eq!(normalize_angle((Kind(10), Kind(10), Kind(1))), (Kind(1), Kind(10), Kind(10)));
    }

    #[test]
    fn normalizing_dihedrals() {
        // Checking that equivalent dihedrals are normalized to the same tuple
        assert_eq!(
            normalize_dihedral((Kind(1), Kind(0), Kind(2), Kind(3))),
            normalize_dihedral((Kind(3), Kind(2), Kind(0), Kind(1)))
        );
        assert_eq!(
            normalize_dihedral((Kind(10), Kind(8), Kind(2), Kind(3))),
            normalize_dihedral((Kind(3), Kind(2), Kind(8), Kind(10)))
        );
        assert_eq!(
            normalize_dihedral((Kind(0), Kind(0), Kind(2), Kind(3))),
            normalize_dihedral((Kind(3), Kind(2), Kind(0), Kind(0)))
        );
        assert_eq!(
            normalize_dihedral((Kind(1), Kind(5), Kind(0), Kind(0))),
            normalize_dihedral((Kind(0), Kind(0), Kind(5), Kind(1)))
        );
        assert_eq!(
            normalize_dihedral((Kind(10), Kind(10), Kind(2), Kind(3))),
            normalize_dihedral((Kind(3), Kind(2), Kind(10), Kind(10)))
        );
        assert_eq!(
            normalize_dihedral((Kind(1), Kind(5), Kind(10), Kind(10))),
            normalize_dihedral((Kind(10), Kind(10), Kind(5), Kind(1)))
        );
        assert_eq!(
            normalize_dihedral((Kind(0), Kind(0), Kind(0), Kind(3))),
            normalize_dihedral((Kind(3), Kind(0), Kind(0), Kind(0)))
        );
        assert_eq!(
            normalize_dihedral((Kind(1), Kind(5), Kind(5), Kind(5))),
            normalize_dihedral((Kind(5), Kind(5), Kind(5), Kind(1)))
        );
        assert_eq!(
            normalize_dihedral((Kind(0), Kind(0), Kind(3), Kind(0))),
            normalize_dihedral((Kind(0), Kind(3), Kind(0), Kind(0)))
        );
        assert_eq!(
            normalize_dihedral((Kind(5), Kind(1), Kind(5), Kind(5))),
            normalize_dihedral((Kind(5), Kind(5), Kind(1), Kind(5)))
        );

        // Checking the actual normalized value
        assert_eq!(
            normalize_dihedral((Kind(1), Kind(0), Kind(2), Kind(3))),
            (Kind(1), Kind(0), Kind(2), Kind(3))
        );
        assert_eq!(
            normalize_dihedral((Kind(10), Kind(8), Kind(2), Kind(3))),
            (Kind(3), Kind(2), Kind(8), Kind(10))
        );
        assert_eq!(
            normalize_dihedral((Kind(0), Kind(0), Kind(2), Kind(3))),
            (Kind(0), Kind(0), Kind(2), Kind(3))
        );
        assert_eq!(
            normalize_dihedral((Kind(1), Kind(5), Kind(0), Kind(0))),
            (Kind(0), Kind(0), Kind(5), Kind(1))
        );
        assert_eq!(
            normalize_dihedral((Kind(10), Kind(10), Kind(2), Kind(3))),
            (Kind(3), Kind(2), Kind(10), Kind(10))
        );
        assert_eq!(
            normalize_dihedral((Kind(1), Kind(5), Kind(10), Kind(10))),
            (Kind(1), Kind(5), Kind(10), Kind(10))
        );
        assert_eq!(
            normalize_dihedral((Kind(0), Kind(0), Kind(0), Kind(3))),
            (Kind(0), Kind(0), Kind(0), Kind(3))
        );
        assert_eq!(
            normalize_dihedral((Kind(1), Kind(5), Kind(5), Kind(5))),
            (Kind(1), Kind(5), Kind(5), Kind(5))
        );
        assert_eq!(
            normalize_dihedral((Kind(0), Kind(0), Kind(3), Kind(0))),
            (Kind(0), Kind(0), Kind(3), Kind(0))
        );
        assert_eq!(
            normalize_dihedral((Kind(5), Kind(1), Kind(5), Kind(5))),
            (Kind(5), Kind(1), Kind(5), Kind(5))
        );
    }

    #[test]
    fn pairs() {
        let mut interactions = Interactions::new();
        let pair = PairInteraction::new(Box::new(NullPotential), 0.0);
        interactions.set_pair(("A", "B"), pair.clone());
        assert!(interactions.pair((Kind(0), Kind(1))).is_some());
        assert!(interactions.pair((Kind(1), Kind(0))).is_some());

        assert!(interactions.pair((Kind(0), Kind(0))).is_none());
        interactions.set_pair(("A", "A"), pair);
        assert!(interactions.pair((Kind(0), Kind(0))).is_some());

        // 'out of bounds' kinds
        assert!(interactions.pair((Kind(55), Kind(55))).is_none());
    }

    #[test]
    fn bonds() {
        let mut interactions = Interactions::new();

        interactions.set_bond(("A", "B"), Box::new(NullPotential));
        assert!(interactions.bond((Kind(0), Kind(1))).is_some());
        assert!(interactions.bond((Kind(1), Kind(0))).is_some());

        assert!(interactions.bond((Kind(0), Kind(0))).is_none());
        interactions.set_bond(("A", "A"), Box::new(NullPotential));
        assert!(interactions.bond((Kind(0), Kind(0))).is_some());


        // 'out of bounds' kinds
        assert!(interactions.bond((Kind(55), Kind(55))).is_none());
    }

    #[test]
    fn angles() {
        let mut interactions = Interactions::new();

        interactions.set_angle(("A", "B", "C"), Box::new(NullPotential));
        assert!(interactions.angle((Kind(0), Kind(1), Kind(2))).is_some());
        assert!(interactions.angle((Kind(2), Kind(1), Kind(0))).is_some());

        assert!(interactions.angle((Kind(2), Kind(2), Kind(1))).is_none());
        interactions.set_angle(("C", "C", "B"), Box::new(NullPotential));
        assert!(interactions.angle((Kind(2), Kind(2), Kind(1))).is_some());
        assert!(interactions.angle((Kind(1), Kind(2), Kind(2))).is_some());

        assert!(interactions.angle((Kind(3), Kind(3), Kind(3))).is_none());
        interactions.set_angle(("D", "D", "D"), Box::new(NullPotential));
        assert!(interactions.angle((Kind(3), Kind(3), Kind(3))).is_some());

        // 'out of bounds' kinds
        assert!(interactions.angle((Kind(55), Kind(55), Kind(55))).is_none());
    }

    #[test]
    fn dihedrals() {
        let mut interactions = Interactions::new();

        interactions.set_dihedral(("A", "B", "C", "D"), Box::new(NullPotential));
        assert!(interactions.dihedral((Kind(0), Kind(1), Kind(2), Kind(3))).is_some());
        assert!(interactions.dihedral((Kind(3), Kind(2), Kind(1), Kind(0))).is_some());

        assert!(interactions.dihedral((Kind(0), Kind(0), Kind(1), Kind(2))).is_none());
        interactions.set_dihedral(("A", "A", "B", "C"), Box::new(NullPotential));
        assert!(interactions.dihedral((Kind(0), Kind(0), Kind(1), Kind(2))).is_some());
        assert!(interactions.dihedral((Kind(2), Kind(1), Kind(0), Kind(0))).is_some());

        interactions.set_dihedral(("A", "A", "B", "A"), Box::new(NullPotential));
        assert!(interactions.dihedral((Kind(0), Kind(0), Kind(1), Kind(0))).is_some());
        assert!(interactions.dihedral((Kind(0), Kind(1), Kind(0), Kind(0))).is_some());

        assert!(interactions.dihedral((Kind(4), Kind(4), Kind(4), Kind(4))).is_none());
        interactions.set_dihedral(("E", "E", "E", "E"), Box::new(NullPotential));
        assert!(interactions.dihedral((Kind(4), Kind(4), Kind(4), Kind(4))).is_some());

        // 'out of bounds' kinds
        assert!(interactions.dihedral((Kind(55), Kind(55), Kind(55), Kind(55))).is_none());
    }

    #[test]
    fn test_maximum_cutoff() {
        let mut interactions = Interactions::new();
        // No cutoff
        assert_eq!(interactions.maximum_cutoff(), None);

        // Only pairs
        let pair = PairInteraction::new(Box::new(NullPotential), 10.0);
        interactions.set_pair(("A", "B"), pair.clone());
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
        interactions.set_pair(("A", "B"), pair);
        assert_eq!(interactions.maximum_cutoff(), Some(10.0));
        interactions.globals.push(Box::new(Wolf::new(15.0)));
        assert_eq!(interactions.maximum_cutoff(), Some(15.0));
        interactions.globals.push(Box::new(Wolf::new(1.0)));
        assert_eq!(interactions.maximum_cutoff(), Some(15.0));
    }
}
