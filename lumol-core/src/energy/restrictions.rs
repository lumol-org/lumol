// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

//! Encoding restrictions in interactions.

/// Possible restrictions on the pair interactions.
///
/// In some force fields, interactions are restricted to a subset of the
/// possible pairs. For example, we may want to have Lennard-Jones interactions
/// between atoms in different molecules only, or scale the coulombic
/// interactions if the particles are at 1-4 distance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum PairRestriction {
    /// No pair should be excluded.
    None,
    /// Only apply the interaction to intra-molecular pairs.
    IntraMolecular,
    /// Only apply the interaction to inter-molecular pairs.
    InterMolecular,
    /// Only apply the interaction to pairs which are not in 1-2 position
    /// (separated by one bond).
    Exclude12,
    /// Only apply the interaction to pairs which are not in 1-2 or 1-3 position
    /// (separated by one or two bonds).
    Exclude13,
    /// Only apply the interaction to pairs which are not in 1-2, 1-3 or 1-4
    /// position (separated by one, two or three bonds).
    Exclude14,
    /// Only apply the interaction to pairs which are not in 1-2 or 1-3
    /// position, and scale the interaction for pairs in 1-4 position (separated
    /// by three bonds).
    Scale14(f64),
}

/// Shortest bond path between two particles in a system
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum BondPath {
    /// No bond path exists, the particles are not in the same molecule
    None,
    /// The two particles are the same one
    SameParticle,
    /// The two particles are separated by one bond
    OneBond,
    /// The two particles are separated by two bonds
    TwoBonds,
    /// The two particles are separated by three bonds
    ThreeBonds,
    /// The two particles are in the same molecule and separated by more than three bonds
    Far,
}

/// Restriction information attached to a pair of `Particles` in a `System`.
#[derive(Clone, Copy, Debug)]
pub struct RestrictionInfo {
    /// Is this pair excluded?
    pub excluded: bool,
    /// Scaling factor for the potential. This value is contained between 0 and
    /// 1.
    pub scaling: f64,
}

impl PairRestriction {
    /// Get the restriction at the given [bond `path`][path].
    ///
    /// [path]: ../sys/struct.System.html#method.bond_path
    ///
    /// # Example
    ///
    /// ```
    /// # use lumol_core::energy::{PairRestriction, BondPath};
    /// let restriction = PairRestriction::None;
    /// assert_eq!(restriction.information(BondPath::ThreeBonds).excluded, false);
    /// assert_eq!(restriction.information(BondPath::TwoBonds).scaling, 1.0);
    ///
    /// let restriction = PairRestriction::Exclude13;
    /// assert_eq!(restriction.information(BondPath::TwoBonds).excluded, true);
    /// assert_eq!(restriction.information(BondPath::ThreeBonds).excluded, false);
    ///
    /// let restriction = PairRestriction::Scale14(0.5);
    /// assert_eq!(restriction.information(BondPath::TwoBonds).excluded, true);
    /// assert_eq!(restriction.information(BondPath::ThreeBonds).excluded, false);
    /// assert_eq!(restriction.information(BondPath::TwoBonds).scaling, 1.0);
    /// assert_eq!(restriction.information(BondPath::ThreeBonds).scaling, 0.5);
    /// ```
    pub fn information(&self, path: BondPath) -> RestrictionInfo {
        let are_in_same_molecule = path != BondPath::None;
        let excluded = match *self {
            PairRestriction::None => false,
            PairRestriction::InterMolecular => are_in_same_molecule,
            PairRestriction::IntraMolecular => !are_in_same_molecule,
            PairRestriction::Exclude12 => path == BondPath::OneBond,
            PairRestriction::Exclude13 | PairRestriction::Scale14(..) => {
                path == BondPath::OneBond || path == BondPath::TwoBonds
            }
            PairRestriction::Exclude14 => {
                path == BondPath::OneBond || path == BondPath::TwoBonds || path == BondPath::ThreeBonds
            },
        };

        let scaling = if let PairRestriction::Scale14(scaling) = *self {
            if path == BondPath::ThreeBonds {
                scaling
            } else {
                1.0
            }
        } else {
            1.0
        };

        RestrictionInfo {
            excluded: excluded,
            scaling: scaling,
        }
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Molecule, Particle, System};

    fn testing_system() -> System {
        // Creating 2 pentane molecule
        //    C   C        C   C
        //   / \ / \      / \ / \
        //  C   C   C    C   C   C
        let mut system = System::new();
        let mut pentane = Molecule::new(Particle::new("CH3"));
        pentane.add_particle_bonded_to(0, Particle::new("CH2"));
        pentane.add_particle_bonded_to(1, Particle::new("CH2"));
        pentane.add_particle_bonded_to(2, Particle::new("CH2"));
        pentane.add_particle_bonded_to(3, Particle::new("CH3"));

        system.add_molecule(pentane.clone());
        system.add_molecule(pentane);

        return system;
    }

    #[test]
    fn none() {
        let restriction = PairRestriction::None;
        let system = testing_system();
        for i in 0..10 {
            for j in 0..10 {
                let path = system.bond_path(i, j);
                let info = restriction.information(path);
                assert!(!info.excluded);
                assert_eq!(info.scaling, 1.0);
            }
        }
    }

    #[test]
    fn intra() {
        let restriction = PairRestriction::IntraMolecular;
        let system = testing_system();
        for i in 0..10 {
            for j in 0..10 {
                let path = system.bond_path(i, j);
                let info = restriction.information(path);
                assert_eq!(info.excluded, !system.are_in_same_molecule(i, j));
                assert_eq!(info.scaling, 1.0);
            }
        }
    }

    #[test]
    fn inter() {
        let restriction = PairRestriction::InterMolecular;
        let system = testing_system();
        for i in 0..10 {
            for j in 0..10 {
                let path = system.bond_path(i, j);
                let info = restriction.information(path);
                assert_eq!(info.excluded, system.are_in_same_molecule(i, j));
                assert_eq!(info.scaling, 1.0);
            }
        }
    }

    #[test]
    fn exclude_12() {
        let restriction = PairRestriction::Exclude12;
        let system = testing_system();
        for i in 0..10 {
            for j in 0..10 {
                let path = system.bond_path(i, j);
                assert_eq!(restriction.information(path).scaling, 1.0);
            }
        }

        // Bonds
        assert!(restriction.information(system.bond_path(0, 1)).excluded);
        assert!(restriction.information(system.bond_path(1, 2)).excluded);
        assert!(restriction.information(system.bond_path(7, 8)).excluded);

        // Not excluded
        assert!(!restriction.information(system.bond_path(0, 8)).excluded);
        assert!(!restriction.information(system.bond_path(0, 3)).excluded);
        assert!(!restriction.information(system.bond_path(8, 2)).excluded);
    }

    #[test]
    fn exclude_13() {
        let restriction = PairRestriction::Exclude13;
        let system = testing_system();
        for i in 0..10 {
            for j in 0..10 {
                let path = system.bond_path(i, j);
                assert_eq!(restriction.information(path).scaling, 1.0);
            }
        }

        // Bonds
        assert!(restriction.information(system.bond_path(0, 1)).excluded);
        assert!(restriction.information(system.bond_path(7, 6)).excluded);

        // Angles
        assert!(restriction.information(system.bond_path(0, 2)).excluded);
        assert!(restriction.information(system.bond_path(1, 3)).excluded);
        assert!(restriction.information(system.bond_path(7, 9)).excluded);

        // Not excluded
        assert!(!restriction.information(system.bond_path(4, 5)).excluded);
        assert!(!restriction.information(system.bond_path(0, 3)).excluded);
        assert!(!restriction.information(system.bond_path(8, 2)).excluded);
    }

    #[test]
    fn exclude_14() {
        let restriction = PairRestriction::Exclude14;
        let system = testing_system();
        for i in 0..10 {
            for j in 0..10 {
                let path = system.bond_path(i, j);
                assert_eq!(restriction.information(path).scaling, 1.0);
            }
        }

        // Bonds
        assert!(restriction.information(system.bond_path(0, 1)).excluded);
        assert!(restriction.information(system.bond_path(7, 6)).excluded);

        // Angles
        assert!(restriction.information(system.bond_path(0, 2)).excluded);
        assert!(restriction.information(system.bond_path(1, 3)).excluded);
        assert!(restriction.information(system.bond_path(7, 9)).excluded);

        // Dihedrals
        assert!(restriction.information(system.bond_path(0, 3)).excluded);
        assert!(restriction.information(system.bond_path(1, 4)).excluded);
        assert!(restriction.information(system.bond_path(6, 9)).excluded);

        // Not excluded
        assert!(!restriction.information(system.bond_path(4, 5)).excluded);
        assert!(!restriction.information(system.bond_path(0, 4)).excluded);
        assert!(!restriction.information(system.bond_path(8, 2)).excluded);
    }

    #[test]
    fn scale_14() {
        let restriction = PairRestriction::Scale14(0.8);
        let system = testing_system();
        for i in 0..10 {
            for j in 0..10 {
                let path = system.bond_path(i, j);
                if path == BondPath::ThreeBonds {
                    assert_eq!(restriction.information(path).scaling, 0.8);
                } else {
                    assert_eq!(restriction.information(path).scaling, 1.0);
                }
            }
        }

        // Bonds
        assert!(restriction.information(system.bond_path(0, 1)).excluded);
        assert!(restriction.information(system.bond_path(7, 6)).excluded);

        // Angles
        assert!(restriction.information(system.bond_path(0, 2)).excluded);
        assert!(restriction.information(system.bond_path(1, 3)).excluded);
        assert!(restriction.information(system.bond_path(7, 9)).excluded);

        // Dihedrals are not excluded, just scaled
        assert!(!restriction.information(system.bond_path(0, 3)).excluded);
        assert!(!restriction.information(system.bond_path(1, 4)).excluded);
        assert!(!restriction.information(system.bond_path(6, 9)).excluded);

        // Not excluded
        assert!(!restriction.information(system.bond_path(4, 5)).excluded);
        assert!(!restriction.information(system.bond_path(0, 4)).excluded);
        assert!(!restriction.information(system.bond_path(8, 2)).excluded);
    }
}
