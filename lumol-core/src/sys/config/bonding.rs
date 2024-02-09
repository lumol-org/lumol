// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

use std::collections::HashSet;
use std::hash::{Hash, Hasher};
use std::ops::Range;

use crate::{Angle, Bond, BondDistances, Dihedral};
use crate::Array2;


/// The basic building block for a topology. A `Bonding` contains data about
/// the connectivity (bonds, angles, dihedrals) between particles in a single
/// molecule.
#[derive(Debug, Clone)]
pub struct Bonding {
    /// All the bonds in the molecule.
    bonds: HashSet<Bond>,
    /// All the angles in the molecule. Rebuilt as needed from the bond list.
    angles: HashSet<Angle>,
    /// All the dihedral angles in the molecule. Rebuilt as needed from the
    /// bond list.
    dihedrals: HashSet<Dihedral>,
    /// Matrix of bond distances in the molecule. The item at index `i, j`
    /// encode the bond distance between the particles `i + self.first` and
    /// `j + self.first`
    distances: Array2<BondDistances>,
    /// Range of atomic indexes in this molecule.
    range: Range<usize>,
}

impl Bonding {
    /// Create a new `Bonding` containing only the atom i
    pub fn new(i: usize) -> Bonding {
        Bonding {
            bonds: HashSet::new(),
            angles: HashSet::new(),
            dihedrals: HashSet::new(),
            distances: Array2::default((1, 1)),
            range: i..i + 1,
        }
    }

    /// Get the number of atoms in the molecule
    pub fn size(&self) -> usize {
        self.range.len()
    }

    /// Get the first atom of this molecule
    pub fn start(&self) -> usize {
        self.range.start
    }

    /// Get the index of the first atom after this molecule
    pub fn end(&self) -> usize {
        self.range.end
    }

    /// Does this molecule contains the particle `i`
    pub fn contains(&self, i: usize) -> bool {
        self.range.start <= i && i < self.range.end
    }

    /// Hash the bonds in this molecule
    pub(crate) fn hash<H: Hasher + Sized>(&self, hasher: &mut H) {
        let mut bonds = self.bonds.iter()
            .map(|bond| Bond::new(bond.i() - self.start(), bond.j() - self.start()))
            .collect::<Vec<_>>();

        bonds.sort_unstable();
        for bond in &bonds {
            bond.i().hash(hasher);
            bond.j().hash(hasher);
        }
    }

    /// Rebuild the full list of angles and dihedral angles from the list of bonds
    fn rebuild(&mut self) {
        self.angles.clear();
        self.dihedrals.clear();
        for bond1 in &self.bonds {
            // Find angles
            for bond2 in &self.bonds {
                if bond1 == bond2 {
                    continue;
                }

                let angle = if bond1.i() == bond2.j() {
                    Angle::new(bond2.i(), bond2.j(), bond1.j())
                } else if bond1.j() == bond2.i() {
                    Angle::new(bond1.i(), bond1.j(), bond2.j())
                } else if bond1.j() == bond2.j() {
                    Angle::new(bond1.i(), bond1.j(), bond2.i())
                } else if bond1.i() == bond2.i() {
                    Angle::new(bond1.j(), bond1.i(), bond2.j())
                } else {
                    // We will not find any dihedral angle from these bonds
                    continue;
                };
                let _ = self.angles.insert(angle);

                // Find dihedral angles
                for bond3 in &self.bonds {
                    if bond2 == bond3 {
                        continue;
                    }

                    let dihedral = if angle.k() == bond3.i() && angle.j() != bond3.j() {
                        Dihedral::new(angle.i(), angle.j(), angle.k(), bond3.j())
                    } else if angle.k() == bond3.j() && angle.j() != bond3.i() {
                        Dihedral::new(angle.i(), angle.j(), angle.k(), bond3.i())
                    } else if angle.i() == bond3.j() && angle.j() != bond3.i() {
                        Dihedral::new(bond3.i(), angle.i(), angle.j(), angle.k())
                    } else if angle.i() == bond3.i() && angle.j() != bond3.j() {
                        Dihedral::new(bond3.j(), angle.i(), angle.j(), angle.k())
                    } else {
                        // (angle.k == bond3.i || angle.k == bond3.j) is an
                        // improper dihedral.
                        continue;
                    };
                    let _ = self.dihedrals.insert(dihedral);
                }
            }
        }
        self.rebuild_connections();
    }

    /// Recompute the connectivity matrix from the bonds, angles and dihedrals
    /// in the system.
    fn rebuild_connections(&mut self) {
        let n = self.size();
        self.distances = Array2::default((n, n));

        let first = self.start();
        let distances = &mut self.distances;
        let mut add_distance_term = |i, j, term| {
            let old_distance = distances[(i - first, j - first)];
            distances[(i - first, j - first)] = old_distance | term;
        };

        for bond in &self.bonds {
            add_distance_term(bond.i(), bond.j(), BondDistances::ONE);
            add_distance_term(bond.j(), bond.i(), BondDistances::ONE);
        }

        for angle in &self.angles {
            add_distance_term(angle.i(), angle.k(), BondDistances::TWO);
            add_distance_term(angle.k(), angle.i(), BondDistances::TWO);
        }

        for dihedral in &self.dihedrals {
            add_distance_term(dihedral.i(), dihedral.m(), BondDistances::THREE);
            add_distance_term(dihedral.m(), dihedral.i(), BondDistances::THREE);
        }
    }

    /// Merge this molecule with `other`. The first particle in `other` should
    /// be the particle just after the last one in `self`.
    pub fn merge_with(&mut self, other: Bonding) {
        assert_eq!(self.range.end, other.range.start);
        self.range.end = other.range.end;
        for bond in other.bonds {
            let _ = self.bonds.insert(bond);
        }

        // Get angles and dihedrals from the other molecule, there is no need to
        // rebuild these.
        for angle in other.angles {
            let _ = self.angles.insert(angle);
        }

        for dihedral in other.dihedrals {
            let _ = self.dihedrals.insert(dihedral);
        }

        self.rebuild_connections();
    }

    /// Translate all indexes in this molecule by `delta`.
    pub fn translate_by(&mut self, delta: isize) {
        if delta < 0 {
            // We should not create negative indexes
            assert!(delta.unsigned_abs() < self.start());
        }

        // The wrapping_add are necessary here, and produce the right result,
        // thanks to integer overflow and the conversion below.
        let delta = delta as usize;
        self.range.start = self.range.start.wrapping_add(delta);
        self.range.end = self.range.end.wrapping_add(delta);

        let mut new_bonds = HashSet::new();
        for bond in &self.bonds {
            let _ = new_bonds.insert(
                Bond::new(bond.i().wrapping_add(delta), bond.j().wrapping_add(delta)),
            );
        }
        self.bonds = new_bonds;

        let mut new_angles = HashSet::new();
        for angle in &self.angles {
            let _ = new_angles.insert(Angle::new(
                angle.i().wrapping_add(delta),
                angle.j().wrapping_add(delta),
                angle.k().wrapping_add(delta),
            ));
        }
        self.angles = new_angles;

        let mut new_dihedrals = HashSet::new();
        for dihedral in &self.dihedrals {
            let _ = new_dihedrals.insert(Dihedral::new(
                dihedral.i().wrapping_add(delta),
                dihedral.j().wrapping_add(delta),
                dihedral.k().wrapping_add(delta),
                dihedral.m().wrapping_add(delta),
            ));
        }
        self.dihedrals = new_dihedrals;
    }

    /// Add a bond between the particles at indexes `i` and `j`. These particles
    /// are assumed to be in the molecule
    pub fn add_bond(&mut self, i: usize, j: usize) {
        assert!(self.contains(i));
        assert!(self.contains(j));
        assert_ne!(i, j);
        let _ = self.bonds.insert(Bond::new(i, j));
        self.rebuild();
    }

    /// Removes particle at index `i` and any associated bonds, angle or
    /// dihedral. This function also update the indexes for the
    /// bonds/angles/dihedral by remove 1 to all the values `> i`
    pub fn remove_particle(&mut self, i: usize) {
        assert!(self.contains(i));
        // Remove bonds containing the particle `i`
        let mut new_bonds = HashSet::new();
        for bond in self.bonds() {
            if bond.i() == i || bond.j() == i {
                continue;
            }

            let (mut bond_i, mut bond_j) = (bond.i(), bond.j());
            if bond_i > i {
                bond_i -= 1;
            }
            if bond_j > i {
                bond_j -= 1;
            }

            let _ = new_bonds.insert(Bond::new(bond_i, bond_j));
        }

        self.bonds = new_bonds;
        self.range.end -= 1;
        self.rebuild();
    }

    /// Get the internal list of bonds
    pub fn bonds(&self) -> &HashSet<Bond> {
        &self.bonds
    }

    /// Get the internal list of angles
    pub fn angles(&self) -> &HashSet<Angle> {
        &self.angles
    }

    /// Get the internal list of dihedrals
    pub fn dihedrals(&self) -> &HashSet<Dihedral> {
        &self.dihedrals
    }

    /// Get the all the possible bond paths the particles `i` and `j` in this molecule
    pub fn bond_distances(&self, i: usize, j: usize) -> BondDistances {
        assert!(self.contains(i) && self.contains(j));
        return self.distances[(i - self.start(), j - self.start())];
    }

    /// Get the indexes of the particles in this molecule. All atoms in the
    /// returned range are inside this molecule.
    pub fn indexes(&self) -> Range<usize> {
        self.range.clone()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::{Angle, Bond, BondDistances, Dihedral};

    #[test]
    fn translate_by() {
        let mut bonding = Bonding::new(0);
        for i in 1..4 {
            bonding.merge_with(Bonding::new(i));
        }
        bonding.add_bond(0, 1);
        bonding.add_bond(1, 2);
        bonding.add_bond(2, 3);

        assert_eq!(bonding.start(), 0);
        assert_eq!(bonding.end(), 4);
        assert!(bonding.bonds().contains(&Bond::new(2, 3)));
        assert!(bonding.angles().contains(&Angle::new(1, 2, 3)));
        assert!(bonding.dihedrals().contains(&Dihedral::new(0, 1, 2, 3)));

        bonding.translate_by(5);
        assert_eq!(bonding.start(), 5);
        assert_eq!(bonding.end(), 9);
        assert!(bonding.bonds().contains(&Bond::new(7, 8)));
        assert!(bonding.angles().contains(&Angle::new(6, 7, 8)));
        assert!(bonding.dihedrals().contains(&Dihedral::new(5, 6, 7, 8)));

        bonding.translate_by(-3);
        assert_eq!(bonding.start(), 2);
        assert_eq!(bonding.end(), 6);
        assert!(bonding.bonds().contains(&Bond::new(4, 5)));
        assert!(bonding.angles().contains(&Angle::new(3, 4, 5)));
        assert!(bonding.dihedrals().contains(&Dihedral::new(2, 3, 4, 5)));
    }

    #[test]
    fn bonding() {
        // Create ethane like this
        //       H    H               4    5
        //       |    |               |    |
        //   H - C -- C - H       3 - 0 -- 1 - 6
        //       |    |               |    |
        //       H    H               2    7
        let mut bonding = Bonding::new(0);

        for i in 1..8 {
            bonding.merge_with(Bonding::new(i));
        }

        bonding.add_bond(0, 1);
        bonding.add_bond(0, 2);
        bonding.add_bond(0, 3);
        bonding.add_bond(0, 4);
        bonding.add_bond(1, 5);
        bonding.add_bond(1, 6);
        bonding.add_bond(1, 7);

        assert_eq!(bonding.start(), 0);
        assert_eq!(bonding.end(), 8);

        assert_eq!(bonding.size(), 8);

        let bonds = vec![
            Bond::new(0, 1),
            Bond::new(0, 2),
            Bond::new(0, 3),
            Bond::new(0, 4),
            Bond::new(1, 5),
            Bond::new(1, 6),
            Bond::new(1, 7),
        ];

        assert_eq!(bonding.bonds().len(), bonds.len());
        for bond in &bonds {
            assert!(bonding.bonds().contains(bond));
        }

        let angles = vec![
            Angle::new(0, 1, 5),
            Angle::new(0, 1, 6),
            Angle::new(0, 1, 7),
            Angle::new(1, 0, 2),
            Angle::new(1, 0, 3),
            Angle::new(1, 0, 4),
            Angle::new(2, 0, 3),
            Angle::new(3, 0, 4),
            Angle::new(2, 0, 4),
            Angle::new(5, 1, 6),
            Angle::new(6, 1, 7),
            Angle::new(5, 1, 7),
        ];

        assert_eq!(bonding.angles().len(), angles.len());
        for angle in &angles {
            assert!(bonding.angles().contains(angle));
        }

        let dihedrals = vec![
            Dihedral::new(2, 0, 1, 5),
            Dihedral::new(2, 0, 1, 6),
            Dihedral::new(2, 0, 1, 7),
            Dihedral::new(3, 0, 1, 5),
            Dihedral::new(3, 0, 1, 6),
            Dihedral::new(3, 0, 1, 7),
            Dihedral::new(4, 0, 1, 5),
            Dihedral::new(4, 0, 1, 6),
            Dihedral::new(4, 0, 1, 7),
        ];

        assert_eq!(bonding.dihedrals().len(), dihedrals.len());
        for dihedral in &dihedrals {
            assert!(bonding.dihedrals().contains(dihedral));
        }

        assert!(bonding.bond_distances(0, 1).contains(BondDistances::ONE));
        assert!(bonding.bond_distances(1, 0).contains(BondDistances::ONE));

        assert!(bonding.bond_distances(0, 7).contains(BondDistances::TWO));
        assert!(bonding.bond_distances(7, 0).contains(BondDistances::TWO));

        assert!(bonding.bond_distances(3, 5).contains(BondDistances::THREE));
        assert!(bonding.bond_distances(5, 3).contains(BondDistances::THREE));

        bonding.remove_particle(6);
        assert_eq!(bonding.bonds().len(), 6);
        assert_eq!(bonding.angles().len(), 9);
        assert_eq!(bonding.dihedrals().len(), 6);
    }

    #[test]
    fn cyclic() {
        //   0 -- 1
        //   |    |
        //   3 -- 2
        let mut bonding = Bonding::new(0);

        for i in 1..4 {
            bonding.merge_with(Bonding::new(i));
        }
        bonding.add_bond(0, 1);
        bonding.add_bond(1, 2);
        bonding.add_bond(2, 3);
        bonding.add_bond(3, 0);

        assert!(bonding.bond_distances(0, 3).contains(BondDistances::ONE));
        assert!(bonding.bond_distances(0, 3).contains(BondDistances::THREE));

        assert!(bonding.angles.contains(&Angle::new(0, 3, 2)));
        assert!(bonding.angles.contains(&Angle::new(0, 1, 2)));
    }

    #[test]
    fn remove_particle() {
        let mut bonding = Bonding::new(0);
        for i in 1..4 {
            bonding.merge_with(Bonding::new(i));
        }

        assert_eq!(bonding.bonds().len(), 0);
        assert_eq!(bonding.size(), 4);

        bonding.add_bond(0, 1);
        bonding.add_bond(2, 3);
        assert_eq!(bonding.bonds().len(), 2);
        assert_eq!(bonding.size(), 4);

        bonding.remove_particle(1);
        assert_eq!(bonding.bonds().len(), 1);
        assert_eq!(bonding.size(), 3);

        bonding.merge_with(Bonding::new(3));
        assert_eq!(bonding.bonds().len(), 1);
        assert_eq!(bonding.size(), 4);
    }
}
