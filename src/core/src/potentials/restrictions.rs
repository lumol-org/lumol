// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

//! Encoding restrictions in interactions.
//!
//! In some force fields, some interactions are restricted to a subset of the
//! possible pairs. For example, we may want to have Lennard-Jones interactions
//! between atoms in different molecules only, or rescale the coulombic
//! interactions if the particles are at 1-4 distance.
use system::System;

/// Possible restrictions on the pair interactions.
///
/// Some pairs can be excluded from the energy computation; this enum lists the
/// possible case of exclusion.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum PairRestriction {
    /// No pair should be excluded
    None,
    /// Only apply the interaction to intra-molecular pairs
    IntraMolecular,
    /// Only apply the interaction to inter-molecular pairs
    InterMolecular,
    /// Only apply the interaction to pairs which are not in 1-2 position
    Exclude12,
    /// Only apply the interaction to pairs which are not in 1-2 or 1-3 position
    Exclude13,
    /// Only apply the interaction to pairs which are not in 1-2, 1-3 or 1-4 position
    Exclude14,
    /// Only apply the interaction to pairs which are not in 1-2 or 1-3 position,
    /// and scale the interaction for pairs in 1-4 position
    Scale14 {
        /// The scaling factor for interactions in 1-4 position
        scaling: f64
    },
}

/// Restriction information attached to a pair of Particles in a System
#[derive(Clone, Copy, Debug)]
pub struct RestrictionInfo {
    /// Is this pair excluded
    pub excluded: bool,
    /// Scaling for this restrictions
    pub scaling: f64
}

impl PairRestriction {
    /// Get the restrictions information for the pair `(i, j)` in the `system`
    pub fn informations(&self, system: &System, i: usize, j: usize) -> RestrictionInfo {
        let distance = system.bond_distance(i, j);
        let are_in_same_molecule = distance >= 0;
        let excluded = match *self {
            PairRestriction::None => false,
            PairRestriction::InterMolecular => are_in_same_molecule,
            PairRestriction::IntraMolecular => !are_in_same_molecule,
            PairRestriction::Exclude12 => distance == 1,
            PairRestriction::Exclude13 | PairRestriction::Scale14{..} => {
                distance == 1 || distance == 2
            },
            PairRestriction::Exclude14 => {
                distance == 1 || distance == 2 || distance == 3
            }
        };

        let scaling = if let PairRestriction::Scale14{scaling} = *self {
            if distance == 3 {scaling} else {1.0}
        } else {
            1.0
        };

        RestrictionInfo {
            excluded: excluded,
            scaling: scaling
        }
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use system::*;

    fn testing_system() -> System {
        // Creating 2 pentane molecule
        //    C   C        C   C
        //   / \ / \      / \ / \
        //  C   C   C    C   C   C
        let mut system = System::new();
        for _ in 0..10 {
            system.add_particle(Particle::new("C"));
        }

        let _ = system.add_bond(0, 1);
        let _ = system.add_bond(1, 2);
        let _ = system.add_bond(2, 3);
        let _ = system.add_bond(3, 4);

        let _ = system.add_bond(5, 6);
        let _ = system.add_bond(6, 7);
        let _ = system.add_bond(7, 8);
        let _ = system.add_bond(8, 9);

        return system;
    }

    #[test]
    fn none() {
        let restrictions = PairRestriction::None;
        let system = testing_system();
        for i in 0..10 {
            for j in 0..10 {
                let info = restrictions.informations(&system, i, j);
                assert_eq!(info.excluded, false);
                assert_eq!(info.scaling, 1.0);
            }
        }
    }

    #[test]
    fn intra() {
        let restrictions = PairRestriction::IntraMolecular;
        let system = testing_system();
        for i in 0..10 {
            for j in 0..10 {
                let info = restrictions.informations(&system, i, j);
                assert_eq!(info.excluded, !system.are_in_same_molecule(i, j));
                assert_eq!(info.scaling, 1.0);
            }
        }
    }

    #[test]
    fn inter() {
        let restrictions = PairRestriction::InterMolecular;
        let system = testing_system();
        for i in 0..10 {
            for j in 0..10 {
                let info = restrictions.informations(&system, i, j);
                assert_eq!(info.excluded, system.are_in_same_molecule(i, j));
                assert_eq!(info.scaling, 1.0);
            }
        }
    }

    #[test]
    fn exclude_12() {
        let restrictions = PairRestriction::Exclude12;
        let system = testing_system();
        for i in 0..10 {
            for j in 0..10 {
                assert_eq!(restrictions.informations(&system, i, j).scaling, 1.0);
            }
        }

        // Bonds
        assert_eq!(restrictions.informations(&system, 0, 1).excluded, true);
        assert_eq!(restrictions.informations(&system, 1, 2).excluded, true);
        assert_eq!(restrictions.informations(&system, 7, 8).excluded, true);

        // Not excluded
        assert_eq!(restrictions.informations(&system, 0, 8).excluded, false);
        assert_eq!(restrictions.informations(&system, 0, 3).excluded, false);
        assert_eq!(restrictions.informations(&system, 8, 2).excluded, false);
    }

    #[test]
    fn exclude_13() {
        let restrictions = PairRestriction::Exclude13;
        let system = testing_system();
        for i in 0..10 {
            for j in 0..10 {
                assert_eq!(restrictions.informations(&system, i, j).scaling, 1.0);
            }
        }

        // Bonds
        assert_eq!(restrictions.informations(&system, 0, 1).excluded, true);
        assert_eq!(restrictions.informations(&system, 7, 6).excluded, true);

        // Angles
        assert_eq!(restrictions.informations(&system, 0, 2).excluded, true);
        assert_eq!(restrictions.informations(&system, 1, 3).excluded, true);
        assert_eq!(restrictions.informations(&system, 7, 9).excluded, true);

        // Not excluded
        assert_eq!(restrictions.informations(&system, 4, 5).excluded, false);
        assert_eq!(restrictions.informations(&system, 0, 3).excluded, false);
        assert_eq!(restrictions.informations(&system, 8, 2).excluded, false);
    }

    #[test]
    fn exclude_14() {
        let restrictions = PairRestriction::Exclude14;
        let system = testing_system();
        for i in 0..10 {
            for j in 0..10 {
                assert_eq!(restrictions.informations(&system, i, j).scaling, 1.0);
            }
        }

        // Bonds
        assert_eq!(restrictions.informations(&system, 0, 1).excluded, true);
        assert_eq!(restrictions.informations(&system, 7, 6).excluded, true);

        // Angles
        assert_eq!(restrictions.informations(&system, 0, 2).excluded, true);
        assert_eq!(restrictions.informations(&system, 1, 3).excluded, true);
        assert_eq!(restrictions.informations(&system, 7, 9).excluded, true);

        // Dihedrals
        assert_eq!(restrictions.informations(&system, 0, 3).excluded, true);
        assert_eq!(restrictions.informations(&system, 1, 4).excluded, true);
        assert_eq!(restrictions.informations(&system, 6, 9).excluded, true);

        // Not excluded
        assert_eq!(restrictions.informations(&system, 4, 5).excluded, false);
        assert_eq!(restrictions.informations(&system, 0, 4).excluded, false);
        assert_eq!(restrictions.informations(&system, 8, 2).excluded, false);
    }

    #[test]
    fn scale_14() {
        let restrictions = PairRestriction::Scale14{scaling: 0.8};
        let system = testing_system();
        for i in 0..10 {
            for j in 0..10 {
                if system.bond_distance(i, j) == 3 {
                    assert_eq!(restrictions.informations(&system, i, j).scaling, 0.8);
                } else {
                    assert_eq!(restrictions.informations(&system, i, j).scaling, 1.0);
                }
            }
        }

        // Bonds
        assert_eq!(restrictions.informations(&system, 0, 1).excluded, true);
        assert_eq!(restrictions.informations(&system, 7, 6).excluded, true);

        // Angles
        assert_eq!(restrictions.informations(&system, 0, 2).excluded, true);
        assert_eq!(restrictions.informations(&system, 1, 3).excluded, true);
        assert_eq!(restrictions.informations(&system, 7, 9).excluded, true);

        // Dihedrals are not excluded, just scaled
        assert_eq!(restrictions.informations(&system, 0, 3).excluded, false);
        assert_eq!(restrictions.informations(&system, 1, 4).excluded, false);
        assert_eq!(restrictions.informations(&system, 6, 9).excluded, false);

        // Not excluded
        assert_eq!(restrictions.informations(&system, 4, 5).excluded, false);
        assert_eq!(restrictions.informations(&system, 0, 4).excluded, false);
        assert_eq!(restrictions.informations(&system, 8, 2).excluded, false);
    }
}
