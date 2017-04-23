// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 Lumol's contributors — BSD license

//! The `Interactions` type for storing particles types to potentials
//! associations.

use std::f64;
use std::collections::BTreeMap;
use std::cmp::{min, max};

use energy::{PairInteraction, BondPotential, AnglePotential, DihedralPotential};
use energy::{GlobalPotential, CoulombicPotential};

use sys::ParticleKind as Kind;

/// Manage associations between particles names and kinds
#[derive(Clone, Debug)]
struct ParticleKinds {
    /// Associating particles names to kinds
    kinds: BTreeMap<String, Kind>,
    /// Associating particles kinds to names
    names: BTreeMap<Kind, String>,
}

impl ParticleKinds {
    /// Create an new and empty particle names <=> kind map
    fn new() -> ParticleKinds {
        ParticleKinds {
            kinds: BTreeMap::new(),
            names: BTreeMap::new(),
        }
    }

    /// Get or create the kind corresponding to the given `name`
    fn kind(&mut self, name: &str) -> Kind {
        debug_assert_eq!(self.names.len(), self.kinds.len());
        if self.kinds.get(name).is_none() {
            let kind = Kind(self.kinds.len() as u32);
            let _ = self.kinds.insert(String::from(name), kind);
            let _ = self.names.insert(kind, String::from(name));
        }
        *self.kinds.get(name).expect("Internal error: unknown particle kind")
    }

    /// Get the name corresponding to a given kind, or `None` if this kind is
    /// unknown
    fn name(&self, kind: Kind) -> Option<String> {
        debug_assert_eq!(self.names.len(), self.kinds.len());
        self.names.get(&kind).cloned()
    }

    /// Get a list of all the known particles kinds.
    fn all_kinds(&self) -> Vec<Kind> {
        self.kinds.values().cloned().collect()
    }
}

/// Sort pair indexes to get a canonical representation
#[inline] fn normalize_pair(i: Kind, j: Kind) -> (Kind, Kind) {
    if i < j { (i, j) } else { (j, i) }
}

/// Sort angle indexes to get a canonical representation
#[inline] fn normalize_angle(i: Kind, j: Kind, k: Kind) -> (Kind, Kind, Kind) {
    if i < k {
        (i, j, k)
    } else {
        (k, j, i)
    }
}

/// Sort dihedral indexes to get a canonical representation
#[inline] fn normalize_dihedral(i: Kind, j: Kind, k: Kind, m: Kind) -> (Kind, Kind, Kind, Kind) {
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

/// The Interaction type hold all data about the potentials in the system,
/// indexed by particle type.
#[derive(Clone)]
pub struct Interactions {
    pub pairs_bis: BTreeMap<Kind, BTreeMap<Kind, Vec<PairInteraction>>>,
    /// Pair potentials
    pairs: BTreeMap<PairKind, Vec<PairInteraction>>,
    /// Bond potentials
    bonds: BTreeMap<BondKind, Vec<Box<BondPotential>>>,
    /// Angle potentials
    angles: BTreeMap<AngleKind, Vec<Box<AnglePotential>>>,
    /// Dihedral angles potentials
    dihedrals: BTreeMap<DihedralKind, Vec<Box<DihedralPotential>>>,
    /// Coulombic potential solver
    coulomb: Option<Box<CoulombicPotential>>,
    /// Global potentials
    globals: Vec<Box<GlobalPotential>>,
    /// Particles names <=> kind associations
    kinds: ParticleKinds,
}

impl Default for Interactions {
    fn default() -> Interactions {
        Interactions::new()
    }
}

impl Interactions {
    /// Create a new empty `Interactions`
    pub fn new() -> Interactions {
        Interactions{
            pairs_bis: BTreeMap::new(),
            pairs: BTreeMap::new(),
            bonds: BTreeMap::new(),
            angles: BTreeMap::new(),
            dihedrals: BTreeMap::new(),
            coulomb: None,
            globals: Vec::new(),
            kinds: ParticleKinds::new(),
        }
    }

    /// Get the kind corresponding to the a given particle `name`.
    pub fn get_kind(&mut self, name: &str) -> Kind {
        self.kinds.kind(name)
    }

    /// Get a list of all the known particles kinds.
    pub fn all_kinds(&self) -> Vec<Kind> {
        self.kinds.all_kinds()
    }

    /// Add the `potential` pair interaction for the pair `(i, j)`
    pub fn add_pair(&mut self, i: &str, j: &str, potential: PairInteraction) {
        let (i, j) = (self.get_kind(i), self.get_kind(j));
        let (i, j) = normalize_pair(i, j);
        let pairs = self.pairs.entry((i, j)).or_insert(Vec::new());
        pairs.push(potential.clone());

        let pairs = self.pairs_bis.entry(i).or_insert(BTreeMap::new())
                        .entry(j).or_insert(Vec::new());
        pairs.push(potential);
    }

    /// Add the `potential` bonded interaction for the pair `(i, j)`
    pub fn add_bond(&mut self, i: &str, j: &str, potential: Box<BondPotential>) {
        let (i, j) = (self.get_kind(i), self.get_kind(j));
        let (i, j) = normalize_pair(i, j);
        let bonds = self.bonds.entry((i, j)).or_insert(Vec::new());
        bonds.push(potential);
    }

    /// Add the `potential` angle interaction for the angle `(i, j, k)`
    pub fn add_angle(&mut self, i: &str, j: &str, k: &str, potential: Box<AnglePotential>) {
        let (i, j, k) = (self.get_kind(i), self.get_kind(j), self.get_kind(k));
        let (i, j, k) = normalize_angle(i, j, k);
        let angles = self.angles.entry((i, j, k)).or_insert(Vec::new());
        angles.push(potential);
    }

    /// Add the `potential` dihedral interaction for the dihedral angle `(i, j,
    /// k, m)`
    pub fn add_dihedral(&mut self, i: &str, j: &str, k: &str, m: &str, potential: Box<DihedralPotential>) {
        let (i, j, k, m) = (self.get_kind(i), self.get_kind(j), self.get_kind(k), self.get_kind(m));
        let (i, j, k, m) = normalize_dihedral(i, j, k, m);
        let dihedrals = self.dihedrals.entry((i, j, k, m)).or_insert(Vec::new());
        dihedrals.push(potential);
    }

    /// Set the coulombic interaction for all pairs to `potential`
    pub fn set_coulomb(&mut self, potential: Box<CoulombicPotential>) {
        self.coulomb = Some(potential);
    }

    /// Add the `potential` global interaction
    pub fn add_global(&mut self, potential: Box<GlobalPotential>) {
        self.globals.push(potential);
    }
}

static NO_PAIR_INTERACTION: &'static [PairInteraction] = &[];
static NO_BOND_INTERACTION: &'static [Box<BondPotential>] = &[];
static NO_ANGLE_INTERACTION: &'static [Box<AnglePotential>] = &[];
static NO_DIHEDRAL_INTERACTION: &'static [Box<DihedralPotential>] = &[];

impl Interactions {
    /// Get all pair interactions corresponding to the pair `(i, j)`
    pub fn pairs(&self, i: Kind, j: Kind) -> &[PairInteraction] {
        let (i, j) = normalize_pair(i, j);
        if let Some(val) = self.pairs.get(&(i, j)) {
            val
        } else {
            let name_i = self.kinds.name(i).unwrap_or(format!("kind {}", i));
            let name_j = self.kinds.name(j).unwrap_or(format!("kind {}", j));
            warn_once!(
                "No potential defined for the pair ({}, {})", name_i, name_j
            );
            NO_PAIR_INTERACTION
        }
    }

    pub fn pair_interactions(&self, i: Kind) -> Option<&BTreeMap<Kind, Vec<PairInteraction>>> {
        self.pairs_bis.get(&i)
    }

    /// Get all pair interactions
    ///
    /// Using this function, you lose information
    /// about which interaction belongs to which pair of particles.
    pub fn all_pairs(&self) -> Vec<&PairInteraction> {
        self.pairs.values()        // get an iterator over map values
                  .flat_map(|i| i) // flatten nested vectors
                  .collect()       // collect into Vec<&PairInteraction>
    }

    /// Get all bonded interactions corresponding to the pair `(i, j)`
    pub fn bonds(&self, i: Kind, j: Kind) -> &[Box<BondPotential>] {
        let (i, j) = normalize_pair(i, j);
        if let Some(val) = self.bonds.get(&(i, j)) {
            val
        } else {
            let name_i = self.kinds.name(i).unwrap_or(format!("kind {}", i));
            let name_j = self.kinds.name(j).unwrap_or(format!("kind {}", j));
            warn_once!(
                "No potential defined for the bond ({}, {})", name_i, name_j
            );
            NO_BOND_INTERACTION
        }
    }

    /// Get all angle interactions corresponding to the angle `(i, j, k)`
    pub fn angles(&self, i: Kind, j:Kind, k:Kind) -> &[Box<AnglePotential>] {
        let (i, j, k) = normalize_angle(i, j, k);
        if let Some(val) = self.angles.get(&(i, j, k)) {
            val
        } else {
            let name_i = self.kinds.name(i).unwrap_or(format!("kind {}", i));
            let name_j = self.kinds.name(j).unwrap_or(format!("kind {}", j));
            let name_k = self.kinds.name(k).unwrap_or(format!("kind {}", k));
            warn_once!(
                "No potential defined for the angle ({}, {}, {})",
                name_i, name_j, name_k
            );
            NO_ANGLE_INTERACTION
        }
    }

    /// Get all dihedral interactions corresponding to the dihedral `(i, j, k, m)`
    pub fn dihedrals(&self, i: Kind, j: Kind, k: Kind, m: Kind) -> &[Box<DihedralPotential>] {
        let (i, j, k, m) = normalize_dihedral(i, j, k, m);
        if let Some(val) = self.dihedrals.get(&(i, j, k, m)) {
            val
        } else {
            let name_i = self.kinds.name(i).unwrap_or(format!("kind {}", i));
            let name_j = self.kinds.name(j).unwrap_or(format!("kind {}", j));
            let name_k = self.kinds.name(k).unwrap_or(format!("kind {}", k));
            let name_m = self.kinds.name(m).unwrap_or(format!("kind {}", m));
            warn_once!(
                "No potential defined for the dihedral ({}, {}, {}, {})",
                name_i, name_j, name_k, name_m
            );
            NO_DIHEDRAL_INTERACTION
        }
    }

    /// Get the coulombic interaction
    pub fn coulomb(&self) -> Option<&Box<CoulombicPotential>> {
        self.coulomb.as_ref()
    }

    /// Get all global interactions
    pub fn globals(&self) -> &[Box<GlobalPotential>] {
        &self.globals
    }

    /// Get maximum cutoff from `coulomb`, `pairs` and `global` interactons.
    pub fn maximum_cutoff(&self) -> Option<f64> {
        // Coulomb potential, return cutoff
        let coulomb_cutoff = match self.coulomb {
            Some(ref coulomb) => match coulomb.cutoff() {
                Some(cutoff) => cutoff,
                None => f64::NAN,
            },
            None => f64::NAN,
        };

        // Go through global interactions, return maximum cutoff
        let global_cutoff = self.globals()
            .iter()
            .map(|i| i.cutoff())
            .filter_map(|rc| rc)
            .fold(f64::NAN, f64::max);

        let mut maximum_cutoff = f64::max(global_cutoff, coulomb_cutoff);

        // Pair interactions, return maximum cutoff
        let pairs_cutoff = self.all_pairs()
            .iter()
            .map(|i| i.cutoff())
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

    use energy::{Harmonic, Wolf, Ewald, SharedEwald, LennardJones, PairInteraction};
    use sys::ParticleKind as Kind;

    #[test]
    fn kinds() {
        let mut kinds = super::ParticleKinds::new();

        assert!(kinds.name(Kind(0)).is_none());
        assert_eq!(kinds.kind("JK"), Kind(0));
        assert_eq!(kinds.kind("H"), Kind(1));

        assert_eq!(kinds.name(Kind(0)), Some(String::from("JK")));
        assert_eq!(kinds.name(Kind(1)), Some(String::from("H")));
    }

    #[test]
    fn sorting_pairs() {
        assert_eq!(normalize_pair(Kind(0), Kind(1)), normalize_pair(Kind(1), Kind(0)));
        assert_eq!(normalize_pair(Kind(125), Kind(0)), normalize_pair(Kind(0), Kind(125)));

        assert_eq!(normalize_pair(Kind(0), Kind(1)), (Kind(0), Kind(1)));
        assert_eq!(normalize_pair(Kind(125), Kind(0)), (Kind(0), Kind(125)));
    }

    #[test]
    fn sorting_angles() {
        // Checking that equivalent angles are normalized to the same tuple
        assert_eq!(normalize_angle(Kind(1), Kind(0), Kind(2)),
                   normalize_angle(Kind(2), Kind(0), Kind(1)));
        assert_eq!(normalize_angle(Kind(15), Kind(4), Kind(8)),
                   normalize_angle(Kind(8), Kind(4), Kind(15)));
        assert_eq!(normalize_angle(Kind(0), Kind(0), Kind(1)),
                   normalize_angle(Kind(1), Kind(0), Kind(0)));
        assert_eq!(normalize_angle(Kind(10), Kind(10), Kind(1)),
                   normalize_angle(Kind(1), Kind(10), Kind(10)));

        // Checking the actual normalized value
        assert_eq!(normalize_angle(Kind(1), Kind(0), Kind(2)),
                   (Kind(1), Kind(0), Kind(2)));
        assert_eq!(normalize_angle(Kind(15), Kind(4), Kind(8)),
                   (Kind(8), Kind(4), Kind(15)));
        assert_eq!(normalize_angle(Kind(0), Kind(0), Kind(1)),
                   (Kind(0), Kind(0), Kind(1)));
        assert_eq!(normalize_angle(Kind(10), Kind(10), Kind(1)),
                   (Kind(1), Kind(10), Kind(10)));
    }

    #[test]
    fn sorting_dihedrals() {
        // Checking that equivalent dihedrals are normalized to the same tuple
        assert_eq!(normalize_dihedral(Kind(1), Kind(0), Kind(2), Kind(3)),
                   normalize_dihedral(Kind(3), Kind(2), Kind(0), Kind(1)));
        assert_eq!(normalize_dihedral(Kind(10), Kind(8), Kind(2), Kind(3)),
                   normalize_dihedral(Kind(3), Kind(2), Kind(8), Kind(10)));
        assert_eq!(normalize_dihedral(Kind(0), Kind(0), Kind(2), Kind(3)),
                   normalize_dihedral(Kind(3), Kind(2), Kind(0), Kind(0)));
        assert_eq!(normalize_dihedral(Kind(1), Kind(5), Kind(0), Kind(0)),
                   normalize_dihedral(Kind(0), Kind(0), Kind(5), Kind(1)));
        assert_eq!(normalize_dihedral(Kind(10), Kind(10), Kind(2), Kind(3)),
                   normalize_dihedral(Kind(3), Kind(2), Kind(10), Kind(10)));
        assert_eq!(normalize_dihedral(Kind(1), Kind(5), Kind(10), Kind(10)),
                   normalize_dihedral(Kind(10), Kind(10), Kind(5), Kind(1)));
        assert_eq!(normalize_dihedral(Kind(0), Kind(0), Kind(0), Kind(3)),
                   normalize_dihedral(Kind(3), Kind(0), Kind(0), Kind(0)));
        assert_eq!(normalize_dihedral(Kind(1), Kind(5), Kind(5), Kind(5)),
                   normalize_dihedral(Kind(5), Kind(5), Kind(5), Kind(1)));
        assert_eq!(normalize_dihedral(Kind(0), Kind(0), Kind(3), Kind(0)),
                   normalize_dihedral(Kind(0), Kind(3), Kind(0), Kind(0)));
        assert_eq!(normalize_dihedral(Kind(5), Kind(1), Kind(5), Kind(5)),
                   normalize_dihedral(Kind(5), Kind(5), Kind(1), Kind(5)));

        // Checking the actual normalized value
        assert_eq!(normalize_dihedral(Kind(1), Kind(0), Kind(2), Kind(3)),
                   (Kind(1), Kind(0), Kind(2), Kind(3)));
        assert_eq!(normalize_dihedral(Kind(10), Kind(8), Kind(2), Kind(3)),
                   (Kind(3), Kind(2), Kind(8), Kind(10)));
        assert_eq!(normalize_dihedral(Kind(0), Kind(0), Kind(2), Kind(3)),
                   (Kind(0), Kind(0), Kind(2), Kind(3)));
        assert_eq!(normalize_dihedral(Kind(1), Kind(5), Kind(0), Kind(0)),
                   (Kind(0), Kind(0), Kind(5), Kind(1)));
        assert_eq!(normalize_dihedral(Kind(10), Kind(10), Kind(2), Kind(3)),
                   (Kind(3), Kind(2), Kind(10), Kind(10)));
        assert_eq!(normalize_dihedral(Kind(1), Kind(5), Kind(10), Kind(10)),
                   (Kind(1), Kind(5), Kind(10), Kind(10)));
        assert_eq!(normalize_dihedral(Kind(0), Kind(0), Kind(0), Kind(3)),
                   (Kind(0), Kind(0), Kind(0), Kind(3)));
        assert_eq!(normalize_dihedral(Kind(1), Kind(5), Kind(5), Kind(5)),
                   (Kind(1), Kind(5), Kind(5), Kind(5)));
        assert_eq!(normalize_dihedral(Kind(0), Kind(0), Kind(3), Kind(0)),
                   (Kind(0), Kind(0), Kind(3), Kind(0)));
        assert_eq!(normalize_dihedral(Kind(5), Kind(1), Kind(5), Kind(5)),
                   (Kind(5), Kind(1), Kind(5), Kind(5)));
    }

    #[test]
    fn pairs() {
        let mut interactions = Interactions::new();
        let pair = PairInteraction::new(Box::new(Harmonic{x0: 0.0, k: 0.0}), 0.0);
        interactions.add_pair("H", "O", pair.clone());
        assert_eq!(interactions.pairs(Kind(0), Kind(1)).len(), 1);
        assert_eq!(interactions.pairs(Kind(1), Kind(0)).len(), 1);

        assert_eq!(interactions.pairs(Kind(0), Kind(0)).len(), 0);
        interactions.add_pair("H", "H", pair.clone());
        assert_eq!(interactions.pairs(Kind(0), Kind(0)).len(), 1);
        assert_eq!(interactions.all_pairs().len(), 2);
    }

    #[test]
    fn bonds() {
        let mut interactions = Interactions::new();

        interactions.add_bond("H", "O", Box::new(Harmonic{x0: 0.0, k: 0.0}));
        assert_eq!(interactions.bonds(Kind(0), Kind(1)).len(), 1);
        assert_eq!(interactions.bonds(Kind(1), Kind(0)).len(), 1);

        assert_eq!(interactions.bonds(Kind(0), Kind(0)).len(), 0);
        interactions.add_bond("H", "H", Box::new(Harmonic{x0: 0.0, k: 0.0}));
        assert_eq!(interactions.bonds(Kind(0), Kind(0)).len(), 1);
    }

    #[test]
    fn angles() {
        let mut interactions = Interactions::new();

        interactions.add_angle("H", "O", "C", Box::new(Harmonic{x0: 0.0, k: 0.0}));
        assert_eq!(interactions.angles(Kind(0), Kind(1), Kind(2)).len(), 1);
        assert_eq!(interactions.angles(Kind(2), Kind(1), Kind(0)).len(), 1);

        interactions.add_angle("C", "C", "O", Box::new(Harmonic{x0: 0.0, k: 0.0}));
        assert_eq!(interactions.angles(Kind(2), Kind(2), Kind(1)).len(), 1);
        assert_eq!(interactions.angles(Kind(1), Kind(2), Kind(2)).len(), 1);

        interactions.add_angle("N", "N", "N", Box::new(Harmonic{x0: 0.0, k: 0.0}));
        assert_eq!(interactions.angles(Kind(3), Kind(3), Kind(3)).len(), 1);
    }

    #[test]
    fn dihedrals() {
        let mut interactions = Interactions::new();

        interactions.add_dihedral("C", "O", "H", "N", Box::new(Harmonic{x0: 0.0, k: 0.0}));
        assert_eq!(interactions.dihedrals(Kind(0), Kind(1), Kind(2), Kind(3)).len(), 1);
        assert_eq!(interactions.dihedrals(Kind(3), Kind(2), Kind(1), Kind(0)).len(), 1);

        interactions.add_dihedral("C", "C", "O", "H", Box::new(Harmonic{x0: 0.0, k: 0.0}));
        assert_eq!(interactions.dihedrals(Kind(0), Kind(0), Kind(1), Kind(2)).len(), 1);
        assert_eq!(interactions.dihedrals(Kind(2), Kind(1), Kind(0), Kind(0)).len(), 1);

        interactions.add_dihedral("C", "C", "O", "C", Box::new(Harmonic{x0: 0.0, k: 0.0}));
        assert_eq!(interactions.dihedrals(Kind(0), Kind(0), Kind(1), Kind(0)).len(), 1);
        assert_eq!(interactions.dihedrals(Kind(0), Kind(1), Kind(0), Kind(0)).len(), 1);

        interactions.add_dihedral("S", "S", "S", "S", Box::new(Harmonic{x0: 0.0, k: 0.0}));
        assert_eq!(interactions.dihedrals(Kind(4), Kind(4), Kind(4), Kind(4)).len(), 1);
    }

    #[test]
    fn globals() {
        let mut interactions = Interactions::new();
        interactions.add_global(Box::new(Wolf::new(1.0)));
        assert_eq!(interactions.globals().len(), 1);
    }

    #[test]
    fn test_maximum_cutoff() {
        let mut interactions = Interactions::new();
        // No cutoff
        assert_eq!(interactions.maximum_cutoff(), None);

        // Only pairs
        let pair = PairInteraction::new(Box::new(LennardJones{sigma: 3.405, epsilon: 1.0}), 10.0);
        interactions.add_pair("A", "B", pair.clone());
        assert_eq!(interactions.maximum_cutoff(), Some(10.0));

        // Only globals
        let mut interactions = Interactions::new();
        interactions.add_global(Box::new(Wolf::new(1.0)));
        assert_eq!(interactions.maximum_cutoff(), Some(1.0));

        interactions.add_global(Box::new(SharedEwald::new(Ewald::new(5.0, 5))));
        assert_eq!(interactions.maximum_cutoff(), Some(5.0));

        // Only coulomb
        let mut interactions = Interactions::new();
        interactions.set_coulomb(Box::new(Wolf::new(1.0)));
        assert_eq!(interactions.maximum_cutoff(), Some(1.0));

        // All potentials
        interactions.add_pair("A", "B", pair.clone());
        assert_eq!(interactions.maximum_cutoff(), Some(10.0));
        interactions.add_global(Box::new(SharedEwald::new(Ewald::new(15.0, 5))));
        assert_eq!(interactions.maximum_cutoff(), Some(15.0));
        interactions.add_global(Box::new(Wolf::new(1.0)));
        assert_eq!(interactions.maximum_cutoff(), Some(15.0));
    }
}
