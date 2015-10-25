/* Cymbalum, Molecular Simulation in Rust - Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
 */

//! `Universe` type definition and implementation.

use std::collections::HashMap;
use std::ops::{Index, IndexMut};
use std::slice;
use std::cmp::{min, max};
use std::iter::IntoIterator;

extern crate chemfiles;
use self::chemfiles::{Trajectory, Frame};

use potentials::{PairPotential, AnglePotential, DihedralPotential, GlobalPotential};
use types::{Vector3D, Matrix3};

use super::Particle;
use super::Molecule;
use super::UnitCell;
use super::interactions::Interactions;
use super::chemfiles::frame_to_universe;

pub type Permutations = Vec<(usize, usize)>;

/// The Universe type hold all the data about a system. This data contains:
///
///   - an unit cell, containing the system;
///   - a list of particles in the system;
///   - a list of molecules in the system;
///   - a list of interactions, associating particles kinds and potentials
///   - a hash map associating particles names and particles kinds.
///
/// In this implementation, the particles contained in a molecule are guaranted
/// to be contiguous in memory. This allow for faster access when iterating over
/// molecules, and easier molecule removal in the universe.
pub struct Universe {
    /// Unit cell of the universe
    cell: UnitCell,
    /// List of particles in the system
    particles: Vec<Particle>,
    /// Molecules in the universe
    molecules: Vec<Molecule>,
    /// Particles kinds, associating particles names and indexes
    kinds: HashMap<String, u16>,
    /// Interactions is a hash map associating particles kinds and potentials
    interactions: Interactions,
    /// Current step of the simulation
    step: u64,
}

impl Universe {
    /// Create a new empty Universe
    pub fn new() -> Universe {
        Universe{
            particles: Vec::new(),
            molecules: Vec::new(),
            kinds: HashMap::new(),
            interactions: Interactions::new(),
            cell: UnitCell::new(),
            step: 0,
        }
    }

    /// Read a trajectory file and create an Universe from it. For a list of
    /// supported formats, please refer to
    /// [Chemharp](http://chemfiles.readthedocs.org/en/latest/formats.html)
    /// documentation.
    pub fn from_file(path: &str) -> Result<Universe, chemfiles::Error> {
        let mut trajectory = try!(Trajectory::open(path));
        let mut frame = try!(Frame::new(0));
        try!(trajectory.read(&mut frame));
        return frame_to_universe(frame);
    }

    /// Do the same work that the `from_file` function, and guess bonds in the
    /// universe based on the distances between the particles.
    pub fn from_file_auto_bonds(path: &str) -> Result<Universe, chemfiles::Error> {
        let mut trajectory = try!(Trajectory::open(path));
        let mut frame = try!(Frame::new(0));
        try!(trajectory.read(&mut frame));
        try!(frame.guess_topology(true));
        return frame_to_universe(frame);
    }

    /// Create an empty universe with a specific UnitCell
    pub fn from_cell(cell: UnitCell) -> Universe {
        let mut universe = Universe::new();
        universe.set_cell(cell);
        return universe;
    }

    /// Get the current step of the universe
    #[inline] pub fn step(&self) -> u64 {
        self.step
    }

    /// Set the current step of the universe to `step`
    #[inline] pub fn set_step(&mut self, step: u64) {
        self.step = step;
    }

    /// Reset the step of the universe to 0
    pub fn reset_step(&mut self) {
        self.step = 0;
    }

    /// Increment the universe step
    pub fn increment_step(&mut self) {
        self.step += 1;
    }
}

/// Topology and particles related functions
impl Universe {
    /// Get the index of the molecule containing the particle `i`
    fn molecule_id(&self, i: usize) -> usize {
        for (idx, mol) in self.molecules.iter().enumerate() {
            if mol.first() <= i && i <= mol.last() {
                return idx;
            }
        }
        error!("Could not find the molecule id for particle {}", i);
        unreachable!()
    }

    /// Get the molecule containing the particle `i`
    pub fn molecule_containing(&self, i:usize) -> &Molecule {
        let id = self.molecule_id(i);
        return &self.molecules[id];
    }

    /// Check if the particles at indexes `i` and `j` are in the same molecule
    #[inline] pub fn are_in_same_molecule(&self, i: usize, j:usize) -> bool {
        self.molecule_id(i) == self.molecule_id(j)
    }

    /// Get the list of molecules in the universe.
    #[inline] pub fn molecules(&self) -> &Vec<Molecule> {
        &self.molecules
    }

    /// Remove the molecule containing the particle at index `i`
    pub fn remove_molecule_containing(&mut self, i: usize) {
        let id = self.molecule_id(i);
        let molecule = self.molecules.remove(id);
        let first = molecule.first();
        let size = molecule.size();

        for _ in 0..size {
            self.particles.remove(first);
        }

        for molecule in self.molecules.iter_mut().skip(id) {
            molecule.translate_by(-(size as isize));
        }
    }

    /// Add a bond between the particles at indexes `i` and `j`. The particles
    /// should have been added to the system before calling this.
    ///
    /// # Warning
    ///
    /// If the bond is between two particles which are not in the same molecule,
    /// the two molecules are merged together by deplacing particles in the
    /// particles list, and thus invalidate any previously stored index. In
    /// particular, any bond, angle, dihedral or molecule is invalidated.
    ///
    /// If molecules where merged, this function will return `Some(perms)`,
    /// where `perms` contains a list of `(old_index, new_index)` associations.
    /// If no molecules where merged, this function will return `None`.
    pub fn add_bond(&mut self, i: usize, j: usize) -> Option<Permutations> {
        assert!(i <= self.particles.len());
        assert!(j <= self.particles.len());
        debug!("Adding bond between the particles {} and {}, in molecules {} and {}", i, j, self.molecule_id(i), self.molecule_id(j));

        let (i, j, perms) = if !self.are_in_same_molecule(i, j) {
            // Getting copy of the molecules before the merge
            let mol_i = self.molecule_id(i);
            let mol_j = self.molecule_id(j);
            let new_mol_idx = min(mol_i, mol_j);
            let old_mol_idx = max(mol_i, mol_j);
            let new_mol = self.molecules[new_mol_idx].clone();
            let old_mol = self.molecules[old_mol_idx].clone();

            let mut perms = Permutations::new();
            let new_mol_last = new_mol.last();
            let mut old_mol_idx = old_mol.first();

            // If new_mol_last + 1 == old_mol_idx, no one move
            if new_mol_last + 1 != old_mol_idx {
                let size = old_mol.size();
                // Add permutation for the molecule we just moved around
                for i in (new_mol_last + 1)..(new_mol_last + size + 1) {
                    perms.push((old_mol_idx, i));
                    old_mol_idx += 1;
                }

                // Add permutations for molecules that where shifted to make space
                // for the just moved molecule.
                for molecule in self.molecules.iter().skip(new_mol_idx + 2).take(old_mol_idx - new_mol_idx - 1) {
                    for i in molecule {
                        perms.push((i - size, i));
                    }
                }
            }

            // Effective merge
            let delta = self.merge_molecules_containing(i, j);

            // One of the `i` or `j` index is no longer valid, as one molecule
            // has been displaced.
            let (i, j) = if mol_i == new_mol_idx {
                (i, j - delta) // j moved
            } else {
                debug_assert!(mol_j == new_mol_idx);
                (i - delta, j) // i moved
            };
            (i, j, Some(perms))
        } else {
            (i, j, None)
        };

        let id = self.molecule_id(i);
        let molecule = &mut self.molecules[id];
        molecule.add_bond(i, j);

        return perms;
    }

    /// Removes particle at index `i` and any associated bonds, angle or dihedral
    pub fn remove_particle(&mut self, i: usize) {
        let id = self.molecule_id(i);
        self.molecules[id].remove_particle(i);

        for molecule in self.molecules.iter_mut().skip(id + 1) {
            molecule.translate_by(-1);
        }

        self.particles.remove(i);
    }


    /// Insert a particle at the end of the internal list
    pub fn add_particle(&mut self, p: Particle) {
        let mut part = p;
        if part.kind() == u16::max_value() {
            // If no value have been precised, set one from the internal list
            // of particles kinds.
            let kind = self.get_kind(part.name());
            part.set_kind(kind);
        }
        self.particles.push(part);
        self.molecules.push(Molecule::new(self.particles.len() - 1));
    }

    /// Get the number of particles in this universe
    #[inline] pub fn size(&self) -> usize {self.particles.len()}

    /// Get an iterator over the `Particle` in this universe
    #[inline] pub fn iter(&self) -> slice::Iter<Particle> {
        self.particles.iter()
    }

    /// Get a mutable iterator over the `Particle` in this universe
    #[inline] pub fn iter_mut(&mut self) -> slice::IterMut<Particle> {
        self.particles.iter_mut()
    }

    /// Get or create the usize kind index for the name `name` of a particle
    fn get_kind(&mut self, name: &str) -> u16 {
        if self.kinds.contains_key(name) {
            return self.kinds[name];
        } else {
            let index = self.kinds.len() as u16;
            self.kinds.insert(name.to_string(), index);
            return index;
        }
    }

    /// Merge the molecules containing the atoms at indexes `i` and `j` in one
    /// molecule. The molecule are merged into the one with the lower index in
    /// `molecules`.
    ///
    /// For example, if we have
    /// ```
    /// H-H  H-H  H-H  H-H
    /// 0 1  2 3  4 5  6 7
    /// ```
    /// and call `merge_molecules_containing(1, 6)`, the molecules 0 and 2 will
    /// be merged into the molecule 0, and the result will be
    /// ```
    /// 0 1 6 7  2 3  4 5  # Old indexes
    /// H-H-H-H  H-H  H-H
    /// 0 1 2 3  4 5  6 7  # New indexes
    /// ```
    ///
    /// This functions return the deplacement of the move molecule, i.e. in this
    /// example `4`.
    fn merge_molecules_containing(&mut self, i: usize, j: usize) -> usize {
        let mol_i = self.molecule_id(i);
        let mol_j = self.molecule_id(j);

        // Move the particles so that we still have molecules contiguous in
        // memory. The molecules are merged in the one with the smaller index.
        let new_mol_idx = min(mol_i, mol_j);
        let old_mol_idx = max(mol_i, mol_j);

        let mut new_mol = self.molecules[new_mol_idx].clone();
        let old_mol = self.molecules[old_mol_idx].clone();
        assert!(new_mol.last() < old_mol.first());
        debug!("Merging molecules \n{:#?}\n ---\n{:#?}", new_mol, old_mol);

        if new_mol.last() + 1 != old_mol.first() {
            let mut new_index = new_mol.last() + 1;
            for i in old_mol {
                // Remove particles from the old position, and insert it to the
                // new one. The indexes are valid during the movement, because
                // we insert a new particle for each particle removed.
                let particle = self.particles.remove(i);
                self.particles.insert(new_index, particle);
                new_index += 1;
            }
        }

        let mut old_mol = self.molecules[max(mol_i, mol_j)].clone();
        let size = old_mol.size() as isize;

        // translate all indexed in the molecules between new_mol and old_mol
        for molecule in self.molecules.iter_mut().skip(new_mol_idx + 1).take(old_mol_idx - new_mol_idx - 1) {
            molecule.translate_by(size);
        }

        let delta = old_mol.first() - new_mol.last() - 1;
        old_mol.translate_by(- (delta as isize));

        new_mol.merge_with(old_mol);
        self.molecules[new_mol_idx] = new_mol;
        self.molecules.remove(old_mol_idx);

        return delta;
    }
}

static NO_PAIR_INTERACTION: &'static [Box<PairPotential>] = &[];
static NO_ANGLE_INTERACTION: &'static [Box<AnglePotential>] = &[];
static NO_DIHEDRAL_INTERACTION: &'static [Box<DihedralPotential>] = &[];

/// Potentials related functions
impl Universe {
    /// Get the list of pair interaction between the particles at indexes `i`
    /// and `j`.
    pub fn pair_potentials(&self, i: usize, j: usize) -> &[Box<PairPotential>] {
        let ikind = self.particles[i].kind();
        let jkind = self.particles[j].kind();
        match self.interactions.pairs(ikind, jkind) {
            Some(val) => &val,
            None => {
                let i = self.particles[i].name();
                let j = self.particles[j].name();
                // TODO: add and use the warn_once! macro
                warn!("No potential defined for the pair ({}, {})", i, j);
                NO_PAIR_INTERACTION
            }
        }
    }

    /// Get the list of bonded interaction between the particles at indexes `i`
    /// and `j`.
    pub fn bond_potentials(&self, i: usize, j: usize) -> &[Box<PairPotential>] {
        let ikind = self.particles[i].kind();
        let jkind = self.particles[j].kind();
        match self.interactions.bonds(ikind, jkind) {
            Some(val) => &val,
            None => {
                let i = self.particles[i].name();
                let j = self.particles[j].name();
                // TODO: add and use the warn_once! macro
                warn!("No potential defined for the bond ({}, {})", i, j);
                NO_PAIR_INTERACTION
            }
        }
    }

    /// Get the list of angle interaction between the particles at indexes `i`,
    /// `j` and `k`.
    pub fn angle_potentials(&self, i: usize, j: usize, k: usize) -> &[Box<AnglePotential>] {
        let ikind = self.particles[i].kind();
        let jkind = self.particles[j].kind();
        let kkind = self.particles[k].kind();

        match self.interactions.angles(ikind, jkind, kkind) {
            Some(val) => &val,
            None => {
                let i = self.particles[i].name();
                let j = self.particles[j].name();
                let k = self.particles[k].name();
                // TODO: add and use the warn_once! macro
                warn!("No potential defined for the angle ({}, {}, {})", i, j, k);
                NO_ANGLE_INTERACTION
            }
        }
    }

    /// Get the list of dihedral angles interaction between the particles at
    /// indexes `i`, `j`, `k` and `m`.
    pub fn dihedral_potentials(&self, i: usize, j: usize, k: usize, m: usize) -> &[Box<DihedralPotential>] {
        let ikind = self.particles[i].kind();
        let jkind = self.particles[j].kind();
        let kkind = self.particles[k].kind();
        let mkind = self.particles[m].kind();

        match self.interactions.dihedrals(ikind, jkind, kkind, mkind) {
            Some(val) => &val,
            None => {
                let i = self.particles[i].name();
                let j = self.particles[j].name();
                let k = self.particles[k].name();
                let m = self.particles[m].name();
                // TODO: add and use the warn_once! macro
                warn!("No potential defined for the dihedral ({}, {}, {}, {})", i, j, k, m);
                NO_DIHEDRAL_INTERACTION
            }
        }
    }

    /// Get all the global interactions
    pub fn global_potentials(&self) -> &Vec<Box<GlobalPotential>> {
        self.interactions.globals()
    }

    /// Add a pair interaction between the particles with names `names`
    pub fn add_pair_interaction(&mut self, i: &str, j: &str, pot: Box<PairPotential>) {
        let ikind = self.get_kind(i);
        let jkind = self.get_kind(j);

        self.interactions.add_pair(ikind, jkind, pot);
    }

    /// Add a bonded interaction between the particles with names `names`
    pub fn add_bond_interaction(&mut self, i: &str, j: &str, pot: Box<PairPotential>) {
        let ikind = self.get_kind(i);
        let jkind = self.get_kind(j);

        self.interactions.add_bond(ikind, jkind, pot);
    }

    /// Add an angle interaction between the particles with names `names`
    pub fn add_angle_interaction(&mut self, i: &str, j: &str, k: &str, pot: Box<AnglePotential>) {
        let ikind = self.get_kind(i);
        let jkind = self.get_kind(j);
        let kkind = self.get_kind(k);

        self.interactions.add_angle(ikind, jkind, kkind, pot);
    }

    /// Add an angle interaction between the particles with names `names`
    pub fn add_dihedral_interaction(&mut self, i: &str, j: &str, k: &str, m: &str, pot: Box<DihedralPotential>) {
        let ikind = self.get_kind(i);
        let jkind = self.get_kind(j);
        let kkind = self.get_kind(k);
        let mkind = self.get_kind(m);

        self.interactions.add_dihedral(ikind, jkind, kkind, mkind, pot);
    }

    /// Add a global interaction to the universe
    pub fn add_global_interaction(&mut self, potential: Box<GlobalPotential>) {
        self.interactions.add_global(potential);
    }
}

impl<'a> IntoIterator for &'a Universe {
    type Item = &'a Particle;
    type IntoIter = slice::Iter<'a, Particle>;

    #[inline]
    fn into_iter(self) -> slice::Iter<'a, Particle> {
        self.iter()
    }
}

impl<'a> IntoIterator for &'a mut Universe {
    type Item = &'a mut Particle;
    type IntoIter = slice::IterMut<'a, Particle>;

    #[inline]
    fn into_iter(self) -> slice::IterMut<'a, Particle> {
        self.iter_mut()
    }
}

/// UnitCell related functions
impl Universe {
    /// Get a reference to  the universe unit cell
    #[inline] pub fn cell(&self) -> &UnitCell {&self.cell}
    /// Get a mutable reference to  the universe unit cell
    #[inline] pub fn cell_mut(&mut self) -> &mut UnitCell {&mut self.cell}
    /// Set the universe unit cell
    #[inline] pub fn set_cell(&mut self, cell: UnitCell) {self.cell = cell;}

    /// Get the distance between the particles at indexes `i` and `j`
    #[inline] pub fn distance(&self, i: usize, j:usize) -> f64 {
        self.cell.distance(self.particles[i].position(), self.particles[j].position())
    }

    /// Wrap the vector i->j in the cell.
    pub fn wrap_vector(&self, i: usize, j:usize) -> Vector3D {
        let mut res = *self.particles[i].position() - *self.particles[j].position();
        self.cell.wrap_vector(&mut res);
        return res;
    }

    /// Get the angle between the particles `i`, `j` and `k`
    pub fn angle(&self, i: usize, j: usize, k: usize) -> f64 {
        let a = self.particles[i].position();
        let b = self.particles[j].position();
        let c = self.particles[k].position();
        self.cell.angle(a, b, c)
    }

    /// Get the angle and the derivatives of the angle between the particles
    /// `i`, `j` and `k`
    pub fn angle_and_derivatives(&self, i: usize, j: usize, k: usize) -> (f64, Vector3D, Vector3D, Vector3D) {
        let a = self.particles[i].position();
        let b = self.particles[j].position();
        let c = self.particles[k].position();
        self.cell.angle_and_derivatives(a, b, c)
    }

    /// Get the dihedral angle between the particles `i`, `j`, `k` and `m`
    pub fn dihedral(&self, i: usize, j: usize, k: usize, m: usize) -> f64 {
        let a = self.particles[i].position();
        let b = self.particles[j].position();
        let c = self.particles[k].position();
        let d = self.particles[m].position();
        self.cell.dihedral(a, b, c, d)
    }

    /// Get the dihedral angle and the derivatives of the dihedral angle
    /// between the particles `i`, `j`, `k` and `m`
    pub fn dihedral_and_derivatives(&self, i: usize, j: usize, k: usize, m: usize) -> (f64, Vector3D, Vector3D, Vector3D, Vector3D) {
        let a = self.particles[i].position();
        let b = self.particles[j].position();
        let c = self.particles[k].position();
        let d = self.particles[m].position();
        self.cell.dihedral_and_derivatives(a, b, c, d)
    }
}

/******************************************************************************/

use simulation::Compute;
use simulation::{PotentialEnergy, KineticEnergy, TotalEnergy};
use simulation::Temperature;
use simulation::Volume;
use simulation::{Virial, Stress, Pressure};

/// Functions to get pysical properties of an universe.
impl Universe {
    // TODO: This implementation recompute the properties each time. These can
    // be cached somehow.

    /// Get the kinetic energy of the system.
    pub fn kinetic_energy(&self) -> f64 {KineticEnergy.compute(self)}
    /// Get the potential energy of the system.
    pub fn potential_energy(&self) -> f64 {PotentialEnergy.compute(self)}
    /// Get the total energy of the system.
    pub fn total_energy(&self) -> f64 {TotalEnergy.compute(self)}
    /// Get the temperature of the system.
    pub fn temperature(&self) -> f64 {Temperature.compute(self)}

    /// Get the volume of the system.
    pub fn volume(&self) -> f64 {Volume.compute(self)}

    /// Get the tensorial virial of the system.
    pub fn virial(&self) -> Matrix3 {Virial.compute(self)}
    /// Get the pressure of the system, from the virial equation
    pub fn pressure(&self) -> f64 {Pressure.compute(self)}
    /// Get the stress tensor of the system, from the virial equation
    pub fn stress(&self) -> Matrix3 {Stress.compute(self)}
}

/******************************************************************************/
impl Index<usize> for Universe {
    type Output = Particle;
    #[inline]
    fn index(&self, index: usize) -> &Particle {
        &self.particles[index]
    }
}

impl IndexMut<usize> for Universe {
    #[inline]
    fn index_mut(&mut self, index: usize) -> &mut Particle {
        &mut self.particles[index]
    }
}

#[cfg(test)]
mod tests {
    use universe::*;
    use types::*;
    use potentials::*;

    #[test]
    fn step() {
        let mut universe = Universe::new();

        assert_eq!(universe.step(), 0);

        universe.increment_step();
        universe.increment_step();
        universe.increment_step();

        assert_eq!(universe.step(), 3);
        universe.reset_step();
        assert_eq!(universe.step(), 0);
    }

    #[test]
    fn cell() {
        let mut universe = Universe::from_cell(UnitCell::cubic(67.0));
        {
            let cell = universe.cell();
            assert_eq!(cell.a(), 67.0);
        }

        universe.set_cell(UnitCell::cubic(10.0));
        let cell = universe.cell();
        assert_eq!(cell.a(), 10.0);
    }

    #[test]
    fn molecules() {
        let mut universe = Universe::new();
        universe.add_particle(Particle::new("H"));
        universe.add_particle(Particle::new("O"));
        universe.add_particle(Particle::new("O"));
        universe.add_particle(Particle::new("H"));

        assert_eq!(universe.add_bond(0, 1), Some(vec![]));
        assert_eq!(universe.add_bond(2, 3), Some(vec![]));

        assert_eq!(universe.molecules().len(), 2);

        let molecule = universe.molecules()[0].clone();
        assert!(molecule.bonds().contains(&Bond::new(0, 1)));
        let molecule = universe.molecules()[1].clone();
        assert!(molecule.bonds().contains(&Bond::new(2, 3)));

        assert_eq!(universe.add_bond(1, 2), Some(vec![]));
        assert_eq!(universe.molecules().len(), 1);

        let molecule = universe.molecules()[0].clone();
        assert!(molecule.angles().contains(&Angle::new(0, 1, 2)));
        assert!(molecule.angles().contains(&Angle::new(1, 2, 3)));
        assert!(molecule.dihedrals().contains(&Dihedral::new(0, 1, 2, 3)));

        universe.remove_particle(2);
        assert_eq!(universe.molecules().len(), 1);

        universe.remove_molecule_containing(1);
        assert_eq!(universe.molecules().len(), 0);
        assert_eq!(universe.size(), 0);
    }

    #[test]
    fn molecules_permutations() {
        let mut universe = Universe::new();
        universe.add_particle(Particle::new("O"));
        universe.add_particle(Particle::new("H"));
        universe.add_particle(Particle::new("H"));

        assert_eq!(universe.add_bond(0, 2), Some(vec![(2, 1), (1, 2)]));
        assert_eq!(universe.add_bond(0, 2), Some(vec![]));
    }

    #[test]
    fn particles() {
        let mut universe = Universe::new();
        universe.add_particle(Particle::new("O"));
        universe.add_particle(Particle::new("H"));
        universe.add_particle(Particle::new("H"));

        assert_eq!(universe.size(), 3);
        assert_eq!(universe[0].name(), "O");
        assert_eq!(universe[1].name(), "H");
        assert_eq!(universe[2].name(), "H");
    }

    #[test]
    fn distances() {
        let mut universe = Universe::from_cell(UnitCell::cubic(5.0));
        universe.add_particle(Particle::new("O"));
        universe.add_particle(Particle::new("H"));

        universe[0].set_position(Vector3D::new(9.0, 0.0, 0.0));
        universe[1].set_position(Vector3D::new(0.0, 0.0, 0.0));
        assert_eq!(universe.distance(0, 1), 1.0);

        universe.set_cell(UnitCell::new());
        assert_eq!(universe.distance(0, 1), 9.0);
    }

    #[test]
    fn interactions() {
        let mut universe = Universe::new();
        universe.add_particle(Particle::new("He"));

        universe.add_pair_interaction("He", "He", Box::new(LennardJones{sigma: 0.3, epsilon: 2.0}));
        universe.add_pair_interaction("He", "He", Box::new(Harmonic{k: 100.0, x0: 1.1}));
        assert_eq!(universe.pair_potentials(0, 0).len(), 2);

        universe.add_bond_interaction("He", "He", Box::new(Harmonic{k: 100.0, x0: 1.1}));
        assert_eq!(universe.bond_potentials(0, 0).len(), 1);

        universe.add_angle_interaction("He", "He", "He", Box::new(Harmonic{k: 100.0, x0: 1.1}));
        assert_eq!(universe.angle_potentials(0, 0, 0).len(), 1);

        universe.add_dihedral_interaction("He", "He", "He", "He", Box::new(CosineHarmonic::new(0.3, 2.0)));
        assert_eq!(universe.dihedral_potentials(0, 0, 0, 0).len(), 1);

        universe.add_global_interaction(Box::new(Wolf::new(1.0)));
        assert_eq!(universe.global_potentials().len(), 1);
    }

    #[test]
    fn missing_interaction() {
        let mut universe = Universe::new();
        universe.add_particle(Particle::new("He"));
        assert_eq!(universe.pair_potentials(0, 0).len(), 0);
        assert_eq!(universe.bond_potentials(0, 0).len(), 0);
        assert_eq!(universe.angle_potentials(0, 0, 0).len(), 0);
        assert_eq!(universe.dihedral_potentials(0, 0, 0, 0).len(), 0);
    }
}
