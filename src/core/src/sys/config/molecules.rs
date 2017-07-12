// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

//! Data about molecules in the system.
use std::collections::HashSet;
use std::iter::IntoIterator;
use std::ops::Range;
use std::hash::{Hash, Hasher};
use std::collections::hash_map::DefaultHasher;

use types::Array2;
use sys::{Particle, Bond, Angle, Dihedral, BondDistance};

#[derive(Debug, Clone)]
/// A molecule is the basic building block for a topology. It contains data
/// about the connectivity (bonds, angles, dihedrals) in the system.
pub struct Molecule {
    /// All the bonds in the molecule.
    bonds: HashSet<Bond>,
    /// All the angles in the molecule. This is rebuilt as needed from the bond
    /// list.
    angles: HashSet<Angle>,
    /// All the dihedral angles in the molecule. This is rebuilt as needed from
    /// the bond list.
    dihedrals: HashSet<Dihedral>,
    /// Matrix of bond distances in the molecule. The item at index `i, j`
    /// encode the bond distance between the particles `i + self.first` and
    /// `j + self.first`
    distances: Array2<BondDistance>,
    /// Range of atomic indexes in this molecule.
    range: Range<usize>,
    /// Hashed value of the set of bonds in the system
    cached_hash: u64
}

impl Hash for Molecule {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.cached_hash.hash(state);
    }
}

impl Molecule {
    /// Create a new `Molecule` containing only the atom i
    pub fn new(i: usize) -> Molecule {
        Molecule {
            bonds: HashSet::new(),
            angles: HashSet::new(),
            dihedrals: HashSet::new(),
            distances: Array2::default((1, 1)),
            range: i..i+1,
            cached_hash: 0
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

    /// Cache the hash of the bonds
    fn rehash(&mut self) {
        let mut hasher = DefaultHasher::new();
        self.range.len().hash(&mut hasher);

        let mut bonds = self.bonds.iter()
                              .map(|bond| Bond::new(bond.i() - self.start(), bond.j() - self.start()))
                              .collect::<Vec<_>>();
        bonds.sort();
        for bond in &bonds {
            bond.i().hash(&mut hasher);
            bond.j().hash(&mut hasher);
        }
        self.cached_hash = hasher.finish();
    }

    /// Cleanup cached data
    fn cleanup(&mut self) {
        self.angles.clear();
        self.dihedrals.clear();
    }

    /// Rebuild the full list of angles and dihedral angles from the list of bonds
    fn rebuild(&mut self) {
        self.cleanup();
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
                self.angles.insert(angle);

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
                    self.dihedrals.insert(dihedral);
                }
            }
        }
        self.rebuild_connections();
        self.rehash();
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
            add_distance_term(bond.i(), bond.j(), BondDistance::ONE);
            add_distance_term(bond.j(), bond.i(), BondDistance::ONE);
        }

        for angle in &self.angles {
            add_distance_term(angle.i(), angle.k(), BondDistance::TWO);
            add_distance_term(angle.k(), angle.i(), BondDistance::TWO);
        }

        for dihedral in &self.dihedrals {
            add_distance_term(dihedral.i(), dihedral.m(), BondDistance::THREE);
            add_distance_term(dihedral.m(), dihedral.i(), BondDistance::THREE);
        }
    }

    /// Merge this molecule with `other`. The first particle in `other` should
    /// be the particle just after the last one in `self`.
    pub fn merge_with(&mut self, other: Molecule) {
        assert_eq!(self.range.end, other.range.start);
        self.range.end = other.range.end;
        for bond in other.bonds() {
            self.bonds.insert(*bond);
        }

        // Get angles and dihedrals from the other molecule, there is no need to
        // rebuild these.
        for angle in other.angles() {
            self.angles.insert(*angle);
        }

        for dihedral in other.dihedrals() {
            self.dihedrals.insert(*dihedral);
        }

        self.rebuild_connections();
        self.rehash();
    }

    /// Translate all indexes in this molecule by `delta`.
    pub fn translate_by(&mut self, delta: isize) {
        if delta < 0 {
            // We should not create negative indexes
            assert!((delta.abs() as usize) < self.start());
        }

        // The wrapping_add are necessary here, and produce the right result,
        // thanks to integer overflow and the conversion below.
        let delta = delta as usize;
        self.range.start = self.range.start.wrapping_add(delta);
        self.range.end = self.range.end.wrapping_add(delta);

        let mut new_bonds = HashSet::new();
        for bond in &self.bonds {
            new_bonds.insert(Bond::new(
                bond.i().wrapping_add(delta),
                bond.j().wrapping_add(delta)
            ));
        }
        self.bonds = new_bonds;

        let mut new_angles = HashSet::new();
        for angle in &self.angles {
            new_angles.insert(Angle::new(
                angle.i().wrapping_add(delta),
                angle.j().wrapping_add(delta),
                angle.k().wrapping_add(delta)
            ));
        }
        self.angles = new_angles;

        let mut new_dihedrals = HashSet::new();
        for dihedral in &self.dihedrals {
            new_dihedrals.insert(Dihedral::new(
                dihedral.i().wrapping_add(delta),
                dihedral.j().wrapping_add(delta),
                dihedral.k().wrapping_add(delta),
                dihedral.m().wrapping_add(delta)
            ));
        }
        self.dihedrals = new_dihedrals;
    }

    /// Add a bond between the particles at indexes `i` and `j`. These particles
    /// are assumed to be in the molecule
    pub fn add_bond(&mut self, i: usize, j: usize)  {
        assert!(self.contains(i));
        assert!(self.contains(j));
        assert_ne!(i, j);
        self.bonds.insert(Bond::new(i, j));
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

            new_bonds.insert(Bond::new(bond_i, bond_j));
        }

        self.bonds = new_bonds;
        self.range.end -= 1;
        self.rebuild();
    }

    /// Get the internal list of bonds
    #[inline] pub fn bonds(&self) -> &HashSet<Bond> {
        &self.bonds
    }

    /// Get the internal list of angles
    #[inline] pub fn angles(&self) -> &HashSet<Angle> {
        &self.angles
    }

    /// Get the internal list of dihedrals
    #[inline] pub fn dihedrals(&self) -> &HashSet<Dihedral> {
        &self.dihedrals
    }

    /// Get the bond distance between the particles `i` and `j`
    #[inline] pub fn bond_distance(&self, i: usize, j: usize) -> BondDistance {
        assert!(self.contains(i) && self.contains(j));
        return self.distances[(i - self.start(), j - self.start())];
    }

    /// Get an iterator over the particles in the molecule
    #[inline] pub fn iter(&self) -> Range<usize> {
        self.into_iter()
    }
}

/// Get the molecule type of the given `molecule` containing the `particles`.
/// This type can be used to identify all the molecules containing the same
/// bonds and particles (see `System::molecule_type` for more information).
pub fn molecule_type(molecule: &Molecule, particles: &[Particle]) -> u64 {
    assert_eq!(particles.len(), molecule.size());
    let mut hasher = DefaultHasher::new();
    molecule.cached_hash.hash(&mut hasher);
    for particle in particles {
        particle.name().hash(&mut hasher);
    }
    hasher.finish()
}

impl IntoIterator for Molecule {
    type Item = usize;
    type IntoIter = Range<usize>;

    fn into_iter(self) -> Range<usize> {
        self.range
    }
}

impl<'a> IntoIterator for &'a Molecule {
    type Item = usize;
    type IntoIter = Range<usize>;

    fn into_iter(self) -> Range<usize> {
        self.range.clone()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use sys::{Bond, Angle, Dihedral, BondDistance};

    #[test]
    fn translate() {
        let mut molecule = Molecule::new(0);
        for i in 1..4 {
            molecule.merge_with(Molecule::new(i));
        }
        molecule.add_bond(0, 1);
        molecule.add_bond(1, 2);
        molecule.add_bond(2, 3);

        assert_eq!(molecule.start(), 0);
        assert_eq!(molecule.end(), 4);
        assert!(molecule.bonds().contains(&Bond::new(2, 3)));
        assert!(molecule.angles().contains(&Angle::new(1, 2, 3)));
        assert!(molecule.dihedrals().contains(&Dihedral::new(0, 1, 2, 3)));

        molecule.translate_by(5);
        assert_eq!(molecule.start(), 5);
        assert_eq!(molecule.end(), 9);
        assert!(molecule.bonds().contains(&Bond::new(7, 8)));
        assert!(molecule.angles().contains(&Angle::new(6, 7, 8)));
        assert!(molecule.dihedrals().contains(&Dihedral::new(5, 6, 7, 8)));

        molecule.translate_by(-3);
        assert_eq!(molecule.start(), 2);
        assert_eq!(molecule.end(), 6);
        assert!(molecule.bonds().contains(&Bond::new(4, 5)));
        assert!(molecule.angles().contains(&Angle::new(3, 4, 5)));
        assert!(molecule.dihedrals().contains(&Dihedral::new(2, 3, 4, 5)));
    }

    #[test]
    fn bonding() {
        // Create ethane like this
        //       H    H               4    5
        //       |    |               |    |
        //   H - C -- C - H       3 - 0 -- 1 - 6
        //       |    |               |    |
        //       H    H               2    7
        let mut molecule = Molecule::new(0);

        for i in 1..8 {
            molecule.merge_with(Molecule::new(i));
        }

        molecule.add_bond(0, 1);
        molecule.add_bond(0, 2);
        molecule.add_bond(0, 3);
        molecule.add_bond(0, 4);
        molecule.add_bond(1, 5);
        molecule.add_bond(1, 6);
        molecule.add_bond(1, 7);

        assert_eq!(molecule.start(), 0);
        assert_eq!(molecule.end(), 8);

        assert_eq!(molecule.size(), 8);

        /**********************************************************************/
        let bonds = vec![Bond::new(0, 1), Bond::new(0, 2), Bond::new(0, 3),
                         Bond::new(0, 4), Bond::new(1, 5), Bond::new(1, 6),
                         Bond::new(1, 7)];

        assert_eq!(molecule.bonds().len(), bonds.len());
        for bond in &bonds {
            assert!(molecule.bonds().contains(bond));
        }

        /**********************************************************************/
        let angles = vec![Angle::new(0, 1, 5), Angle::new(0, 1, 6), Angle::new(0, 1, 7),
                          Angle::new(1, 0, 2), Angle::new(1, 0, 3), Angle::new(1, 0, 4),
                          Angle::new(2, 0, 3), Angle::new(3, 0, 4), Angle::new(2, 0, 4),
                          Angle::new(5, 1, 6), Angle::new(6, 1, 7), Angle::new(5, 1, 7)];

        assert_eq!(molecule.angles().len(), angles.len());
        for angle in &angles {
            assert!(molecule.angles().contains(angle));
        }

        /**********************************************************************/
        let dihedrals = vec![Dihedral::new(2, 0, 1, 5), Dihedral::new(2, 0, 1, 6),
                             Dihedral::new(2, 0, 1, 7), Dihedral::new(3, 0, 1, 5),
                             Dihedral::new(3, 0, 1, 6), Dihedral::new(3, 0, 1, 7),
                             Dihedral::new(4, 0, 1, 5), Dihedral::new(4, 0, 1, 6),
                             Dihedral::new(4, 0, 1, 7)];

        assert_eq!(molecule.dihedrals().len(), dihedrals.len());
        for dihedral in &dihedrals {
            assert!(molecule.dihedrals().contains(dihedral));
        }

        /**********************************************************************/
        assert!(molecule.bond_distance(0, 1).contains(BondDistance::ONE));
        assert!(molecule.bond_distance(1, 0).contains(BondDistance::ONE));

        assert!(molecule.bond_distance(0, 7).contains(BondDistance::TWO));
        assert!(molecule.bond_distance(7, 0).contains(BondDistance::TWO));

        assert!(molecule.bond_distance(3, 5).contains(BondDistance::THREE));
        assert!(molecule.bond_distance(5, 3).contains(BondDistance::THREE));

        /**********************************************************************/

        molecule.remove_particle(6);
        assert_eq!(molecule.bonds().len(), 6);
        assert_eq!(molecule.angles().len(), 9);
        assert_eq!(molecule.dihedrals().len(), 6);
    }

    #[test]
    fn cyclic() {
        //   0 -- 1
        //   |    |
        //   3 -- 2
        let mut molecule = Molecule::new(0);

        for i in 1..4 {
            molecule.merge_with(Molecule::new(i));
        }
        molecule.add_bond(0, 1);
        molecule.add_bond(1, 2);
        molecule.add_bond(2, 3);
        molecule.add_bond(3, 0);

        assert!(molecule.bond_distance(0, 3).contains(BondDistance::ONE));
        assert!(molecule.bond_distance(0, 3).contains(BondDistance::THREE));

        assert!(molecule.angles.contains(&Angle::new(0, 3, 2)));
        assert!(molecule.angles.contains(&Angle::new(0, 1, 2)));
    }

    #[test]
    fn remove_particle() {
        let mut molecule = Molecule::new(0);
        for i in 1..4 {
            molecule.merge_with(Molecule::new(i));
        }

        assert_eq!(molecule.bonds().len(), 0);
        assert_eq!(molecule.size(), 4);

        molecule.add_bond(0, 1);
        molecule.add_bond(2, 3);
        assert_eq!(molecule.bonds().len(), 2);
        assert_eq!(molecule.size(), 4);

        molecule.remove_particle(1);
        assert_eq!(molecule.bonds().len(), 1);
        assert_eq!(molecule.size(), 3);

        molecule.merge_with(Molecule::new(3));
        assert_eq!(molecule.bonds().len(), 1);
        assert_eq!(molecule.size(), 4);
    }

    #[test]
    fn hash() {
        let mut molecule = Molecule::new(0);
        assert_eq!(molecule.cached_hash, 0);

        for i in 1..4 {
            molecule.merge_with(Molecule::new(i));
        }

        let hash = molecule.cached_hash;
        assert!(hash != 0);

        molecule.add_bond(0, 1);
        molecule.add_bond(2, 3);
        assert!(molecule.cached_hash != hash);
        let hash = molecule.cached_hash;

        molecule.remove_particle(1);
        assert!(molecule.cached_hash != hash);
        let hash = molecule.cached_hash;

        // Hash should be the same when translating the molecule
        molecule.translate_by(67);
        molecule.rehash();
        assert_eq!(molecule.cached_hash, hash);
    }
}
