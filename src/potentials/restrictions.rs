/* Cymbalum, Molecular Simulation in Rust - Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
 */

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
#[derive(Clone, Debug)]
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

impl PairRestriction {
    /// Check wether this restriction exclude the pair `(i, j)` from the
    /// interactions.
    pub fn is_excluded_pair(&self, system: &System, i: usize, j: usize) -> bool {
        match *self {
            PairRestriction::None => false,
            PairRestriction::InterMolecular => {
                system.are_in_same_molecule(i, j)
            },
            PairRestriction::IntraMolecular => {
                !system.are_in_same_molecule(i, j)
            },
            PairRestriction::Exclude12 => {
                let distance = system.shortest_path(i, j);
                distance == 2
            },
            PairRestriction::Exclude13 => {
                let distance = system.shortest_path(i, j);
                distance == 2 || distance == 3
            },
            PairRestriction::Exclude14 | PairRestriction::Scale14{..} => {
                let distance = system.shortest_path(i, j);
                distance == 2 || distance == 3 || distance == 4
            }
        }
    }

    /// Get the scaling factor for the pair `(i, j)`
    pub fn scaling(&self, system: &System, i: usize, j: usize) -> f64 {
        match *self {
            PairRestriction::None           | PairRestriction::Exclude12      |
            PairRestriction::Exclude13      | PairRestriction::Exclude14      |
            PairRestriction::InterMolecular | PairRestriction::IntraMolecular => 1.0,
            PairRestriction::Scale14{scaling} => {
                if !system.are_in_same_molecule(i, j) {
                    1.0
                } else {
                    if system.shortest_path(i, j) == 4 {
                        scaling
                    } else {
                        1.0
                    }
                }
            }
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

        system.add_bond(0, 1);
        system.add_bond(1, 2);
        system.add_bond(2, 3);
        system.add_bond(3, 4);

        system.add_bond(5, 6);
        system.add_bond(6, 7);
        system.add_bond(7, 8);
        system.add_bond(8, 9);

        return system;
    }

    #[test]
    fn none() {
        let restrictions = PairRestriction::None;
        let system = testing_system();
        for i in 0..10 {
            for j in 0..10 {
                assert_eq!(restrictions.is_excluded_pair(&system, i, j), false);
                assert_eq!(restrictions.scaling(&system, i, j), 1.0);
            }
        }
    }

    #[test]
    fn intra() {
        let restrictions = PairRestriction::IntraMolecular;
        let system = testing_system();
        for i in 0..10 {
            for j in 0..10 {
                assert_eq!(restrictions.scaling(&system, i, j), 1.0);
                assert_eq!(
                    restrictions.is_excluded_pair(&system, i, j),
                    !system.are_in_same_molecule(i, j)
                );
            }
        }
    }

    #[test]
    fn inter() {
        let restrictions = PairRestriction::InterMolecular;
        let system = testing_system();
        for i in 0..10 {
            for j in 0..10 {
                assert_eq!(restrictions.scaling(&system, i, j), 1.0);
                assert_eq!(
                    restrictions.is_excluded_pair(&system, i, j),
                    system.are_in_same_molecule(i, j)
                );
            }
        }
    }

    #[test]
    fn exclude_12() {
        let restrictions = PairRestriction::Exclude12;
        let system = testing_system();
        for i in 0..10 {
            for j in 0..10 {
                assert_eq!(restrictions.scaling(&system, i, j), 1.0);
            }
        }

        // Bonds
        assert_eq!(restrictions.is_excluded_pair(&system, 0, 1), true);
        assert_eq!(restrictions.is_excluded_pair(&system, 1, 2), true);
        assert_eq!(restrictions.is_excluded_pair(&system, 7, 8), true);

        // Not excluded
        assert_eq!(restrictions.is_excluded_pair(&system, 0, 8), false);
        assert_eq!(restrictions.is_excluded_pair(&system, 0, 3), false);
        assert_eq!(restrictions.is_excluded_pair(&system, 8, 2), false);
    }

    #[test]
    fn exclude_13() {
        let restrictions = PairRestriction::Exclude13;
        let system = testing_system();
        for i in 0..10 {
            for j in 0..10 {
                assert_eq!(restrictions.scaling(&system, i, j), 1.0);
            }
        }

        // Bonds
        assert_eq!(restrictions.is_excluded_pair(&system, 0, 1), true);
        assert_eq!(restrictions.is_excluded_pair(&system, 7, 6), true);

        // Angles
        assert_eq!(restrictions.is_excluded_pair(&system, 0, 2), true);
        assert_eq!(restrictions.is_excluded_pair(&system, 1, 3), true);
        assert_eq!(restrictions.is_excluded_pair(&system, 7, 9), true);

        // Not excluded
        assert_eq!(restrictions.is_excluded_pair(&system, 4, 5), false);
        assert_eq!(restrictions.is_excluded_pair(&system, 0, 3), false);
        assert_eq!(restrictions.is_excluded_pair(&system, 8, 2), false);
    }

    #[test]
    fn exclude_14() {
        let restrictions = PairRestriction::Exclude14;
        let system = testing_system();
        for i in 0..10 {
            for j in 0..10 {
                assert_eq!(restrictions.scaling(&system, i, j), 1.0);
            }
        }

        // Bonds
        assert_eq!(restrictions.is_excluded_pair(&system, 0, 1), true);
        assert_eq!(restrictions.is_excluded_pair(&system, 7, 6), true);

        // Angles
        assert_eq!(restrictions.is_excluded_pair(&system, 0, 2), true);
        assert_eq!(restrictions.is_excluded_pair(&system, 1, 3), true);
        assert_eq!(restrictions.is_excluded_pair(&system, 7, 9), true);

        // Dihedrals
        assert_eq!(restrictions.is_excluded_pair(&system, 0, 3), true);
        assert_eq!(restrictions.is_excluded_pair(&system, 1, 4), true);
        assert_eq!(restrictions.is_excluded_pair(&system, 6, 9), true);

        // Not excluded
        assert_eq!(restrictions.is_excluded_pair(&system, 4, 5), false);
        assert_eq!(restrictions.is_excluded_pair(&system, 0, 4), false);
        assert_eq!(restrictions.is_excluded_pair(&system, 8, 2), false);
    }

    #[test]
    fn scale_14() {
        let restrictions = PairRestriction::Scale14{scaling: 0.8};
        let system = testing_system();
        for i in 0..10 {
            for j in 0..10 {
                if system.shortest_path(i, j) == 4 {
                    assert_eq!(restrictions.scaling(&system, i, j), 0.8);
                } else {
                    assert_eq!(restrictions.scaling(&system, i, j), 1.0);
                }
            }
        }

        // Bonds
        assert_eq!(restrictions.is_excluded_pair(&system, 0, 1), true);
        assert_eq!(restrictions.is_excluded_pair(&system, 7, 6), true);

        // Angles
        assert_eq!(restrictions.is_excluded_pair(&system, 0, 2), true);
        assert_eq!(restrictions.is_excluded_pair(&system, 1, 3), true);
        assert_eq!(restrictions.is_excluded_pair(&system, 7, 9), true);

        // Dihedrals
        assert_eq!(restrictions.is_excluded_pair(&system, 0, 3), true);
        assert_eq!(restrictions.is_excluded_pair(&system, 1, 4), true);
        assert_eq!(restrictions.is_excluded_pair(&system, 6, 9), true);

        // Not excluded
        assert_eq!(restrictions.is_excluded_pair(&system, 4, 5), false);
        assert_eq!(restrictions.is_excluded_pair(&system, 0, 4), false);
        assert_eq!(restrictions.is_excluded_pair(&system, 8, 2), false);
    }
}
