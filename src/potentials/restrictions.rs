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

use universe::Universe;

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
    pub fn is_excluded_pair(&self, universe: &Universe, i: usize, j: usize) -> bool {
        match *self {
            PairRestriction::None => false,
            PairRestriction::InterMolecular => {
                universe.are_in_same_molecule(i, j)
            },
            PairRestriction::IntraMolecular => {
                !universe.are_in_same_molecule(i, j)
            },
            PairRestriction::Exclude12 => {
                let distance = universe.shortest_path(i, j);
                distance == 2
            },
            PairRestriction::Exclude13 => {
                let distance = universe.shortest_path(i, j);
                distance == 2 || distance == 3
            },
            PairRestriction::Exclude14 | PairRestriction::Scale14{..} => {
                let distance = universe.shortest_path(i, j);
                distance == 2 || distance == 3 || distance == 4
            }
        }
    }

    /// Get the scaling factor for the pair `(i, j)`
    pub fn scaling(&self, universe: &Universe, i: usize, j: usize) -> f64 {
        match *self {
            PairRestriction::None           | PairRestriction::Exclude12      |
            PairRestriction::Exclude13      | PairRestriction::Exclude14      |
            PairRestriction::InterMolecular | PairRestriction::IntraMolecular => 1.0,
            PairRestriction::Scale14{scaling} => {
                if !universe.are_in_same_molecule(i, j) {
                    1.0
                } else {
                    if universe.shortest_path(i, j) == 4 {
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
    use universe::*;

    fn testing_universe() -> Universe {
        // Creating 2 pentane molecule
        //    C   C        C   C
        //   / \ / \      / \ / \
        //  C   C   C    C   C   C
        let mut universe = Universe::new();
        for _ in 0..10 {
            universe.add_particle(Particle::new("C"));
        }

        universe.add_bond(0, 1);
        universe.add_bond(1, 2);
        universe.add_bond(2, 3);
        universe.add_bond(3, 4);

        universe.add_bond(5, 6);
        universe.add_bond(6, 7);
        universe.add_bond(7, 8);
        universe.add_bond(8, 9);

        return universe;
    }

    #[test]
    fn none() {
        let restrictions = PairRestriction::None;
        let universe = testing_universe();
        for i in 0..10 {
            for j in 0..10 {
                assert_eq!(restrictions.is_excluded_pair(&universe, i, j), false);
                assert_eq!(restrictions.scaling(&universe, i, j), 1.0);
            }
        }
    }

    #[test]
    fn intra() {
        let restrictions = PairRestriction::IntraMolecular;
        let universe = testing_universe();
        for i in 0..10 {
            for j in 0..10 {
                assert_eq!(restrictions.scaling(&universe, i, j), 1.0);
                assert_eq!(
                    restrictions.is_excluded_pair(&universe, i, j),
                    !universe.are_in_same_molecule(i, j)
                );
            }
        }
    }

    #[test]
    fn inter() {
        let restrictions = PairRestriction::InterMolecular;
        let universe = testing_universe();
        for i in 0..10 {
            for j in 0..10 {
                assert_eq!(restrictions.scaling(&universe, i, j), 1.0);
                assert_eq!(
                    restrictions.is_excluded_pair(&universe, i, j),
                    universe.are_in_same_molecule(i, j)
                );
            }
        }
    }

    #[test]
    fn exclude_12() {
        let restrictions = PairRestriction::Exclude12;
        let universe = testing_universe();
        for i in 0..10 {
            for j in 0..10 {
                assert_eq!(restrictions.scaling(&universe, i, j), 1.0);
            }
        }

        // Bonds
        assert_eq!(restrictions.is_excluded_pair(&universe, 0, 1), true);
        assert_eq!(restrictions.is_excluded_pair(&universe, 1, 2), true);
        assert_eq!(restrictions.is_excluded_pair(&universe, 7, 8), true);

        // Not excluded
        assert_eq!(restrictions.is_excluded_pair(&universe, 0, 8), false);
        assert_eq!(restrictions.is_excluded_pair(&universe, 0, 3), false);
        assert_eq!(restrictions.is_excluded_pair(&universe, 8, 2), false);
    }

    #[test]
    fn exclude_13() {
        let restrictions = PairRestriction::Exclude13;
        let universe = testing_universe();
        for i in 0..10 {
            for j in 0..10 {
                assert_eq!(restrictions.scaling(&universe, i, j), 1.0);
            }
        }

        // Bonds
        assert_eq!(restrictions.is_excluded_pair(&universe, 0, 1), true);
        assert_eq!(restrictions.is_excluded_pair(&universe, 7, 6), true);

        // Angles
        assert_eq!(restrictions.is_excluded_pair(&universe, 0, 2), true);
        assert_eq!(restrictions.is_excluded_pair(&universe, 1, 3), true);
        assert_eq!(restrictions.is_excluded_pair(&universe, 7, 9), true);

        // Not excluded
        assert_eq!(restrictions.is_excluded_pair(&universe, 4, 5), false);
        assert_eq!(restrictions.is_excluded_pair(&universe, 0, 3), false);
        assert_eq!(restrictions.is_excluded_pair(&universe, 8, 2), false);
    }

    #[test]
    fn exclude_14() {
        let restrictions = PairRestriction::Exclude14;
        let universe = testing_universe();
        for i in 0..10 {
            for j in 0..10 {
                assert_eq!(restrictions.scaling(&universe, i, j), 1.0);
            }
        }

        // Bonds
        assert_eq!(restrictions.is_excluded_pair(&universe, 0, 1), true);
        assert_eq!(restrictions.is_excluded_pair(&universe, 7, 6), true);

        // Angles
        assert_eq!(restrictions.is_excluded_pair(&universe, 0, 2), true);
        assert_eq!(restrictions.is_excluded_pair(&universe, 1, 3), true);
        assert_eq!(restrictions.is_excluded_pair(&universe, 7, 9), true);

        // Dihedrals
        assert_eq!(restrictions.is_excluded_pair(&universe, 0, 3), true);
        assert_eq!(restrictions.is_excluded_pair(&universe, 1, 4), true);
        assert_eq!(restrictions.is_excluded_pair(&universe, 6, 9), true);

        // Not excluded
        assert_eq!(restrictions.is_excluded_pair(&universe, 4, 5), false);
        assert_eq!(restrictions.is_excluded_pair(&universe, 0, 4), false);
        assert_eq!(restrictions.is_excluded_pair(&universe, 8, 2), false);
    }

    #[test]
    fn scale_14() {
        let restrictions = PairRestriction::Scale14{scaling: 0.8};
        let universe = testing_universe();
        for i in 0..10 {
            for j in 0..10 {
                if universe.shortest_path(i, j) == 4 {
                    assert_eq!(restrictions.scaling(&universe, i, j), 0.8);
                } else {
                    assert_eq!(restrictions.scaling(&universe, i, j), 1.0);
                }
            }
        }

        // Bonds
        assert_eq!(restrictions.is_excluded_pair(&universe, 0, 1), true);
        assert_eq!(restrictions.is_excluded_pair(&universe, 7, 6), true);

        // Angles
        assert_eq!(restrictions.is_excluded_pair(&universe, 0, 2), true);
        assert_eq!(restrictions.is_excluded_pair(&universe, 1, 3), true);
        assert_eq!(restrictions.is_excluded_pair(&universe, 7, 9), true);

        // Dihedrals
        assert_eq!(restrictions.is_excluded_pair(&universe, 0, 3), true);
        assert_eq!(restrictions.is_excluded_pair(&universe, 1, 4), true);
        assert_eq!(restrictions.is_excluded_pair(&universe, 6, 9), true);

        // Not excluded
        assert_eq!(restrictions.is_excluded_pair(&universe, 4, 5), false);
        assert_eq!(restrictions.is_excluded_pair(&universe, 0, 4), false);
        assert_eq!(restrictions.is_excluded_pair(&universe, 8, 2), false);
    }
}
