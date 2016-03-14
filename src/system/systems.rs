// Cymbalum, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

//! The base type for simulation data in `cymbalum` is the `System` type.
//!
//! An `System` consists of a list of `Particle`; a list of `Molecule`
//! specifying how the particles are bonded together; an unit cell for boundary
//! conditions; and the interactions between these particles.
use std::collections::HashMap;
use std::ops::{Index, IndexMut};
use std::slice;
use std::cmp::{min, max};
use std::iter::IntoIterator;
use std::u8;
use std::cell::RefCell;

use potentials::{PairPotential, AnglePotential, DihedralPotential};
use potentials::{CoulombicPotential, GlobalPotential};
use potentials::PairRestriction;
use types::{Vector3D, Matrix3};

use super::Particle;
use super::Molecule;
use super::{CONNECT_12, CONNECT_13, CONNECT_14, CONNECT_FAR};
use super::UnitCell;
use super::interactions::{Interactions, PairInteraction};
use super::EnergyEvaluator;
use super::molecules::moltype;

pub type Permutations = Vec<(usize, usize)>;

/// The System type hold all the data about a simulated system.
///
/// This data contains:
///
///   - an unit cell, containing the system;
///   - a list of particles in the system;
///   - a list of molecules in the system;
///   - a list of interactions, associating particles kinds and potentials
///   - a hash map associating particles names and particles kinds.
///
/// In the implementation, the particles contained in a molecule are guaranted
/// to be contiguous in memory. This allow for faster access when iterating over
/// molecules, and easier molecule removal in the system.
#[derive(Clone)]
pub struct System {
    /// Unit cell of the system
    cell: UnitCell,
    /// List of particles in the system
    particles: Vec<Particle>,
    /// Molecules in the system
    molecules: Vec<Molecule>,
    /// Particles kinds, associating particles names and indexes
    kinds: HashMap<String, u16>,
    /// Interactions is a hash map associating particles kinds and potentials
    interactions: Interactions,
    /// Current step of the simulation
    step: u64,
}

impl System {
    /// Create a new empty System
    pub fn new() -> System {
        System{
            particles: Vec::new(),
            molecules: Vec::new(),
            kinds: HashMap::new(),
            interactions: Interactions::new(),
            cell: UnitCell::new(),
            step: 0,
        }
    }

    /// Create an empty system with a specific UnitCell
    pub fn from_cell(cell: UnitCell) -> System {
        let mut system = System::new();
        system.set_cell(cell);
        return system;
    }

    /// Get the current step of the system
    #[inline] pub fn step(&self) -> u64 {
        self.step
    }

    /// Increment the system step
    pub fn increment_step(&mut self) {
        self.step += 1;
    }
}

/// Topology and particles related functions
impl System {
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

    /// Get the type of the molecule at index `molid`. This type is a hash of
    /// the atoms names, and the set of bonds in the molecule. This means that
    /// two molecules will have the same type if and only if they contains the
    /// same atoms and the same bonds, **in the same order**.
    pub fn moltype(&self, molid: usize) -> u64 {
        let molecule = self.molecule(molid);
        moltype(molecule, &self.particles[molecule.into_iter()])
    }

    /// Get a list of molecules with `moltype` molecule type.
    pub fn molecules_with_moltype(&self, moltype: u64) -> Vec<usize> {
        let mut res = Vec::new();
        for i in 0..self.molecules().len() {
            if self.moltype(i) == moltype {
                res.push(i);
            }
        }
        return res;
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

    /// Get the list of molecules in the system.
    #[inline] pub fn molecules(&self) -> &[Molecule] {
        &self.molecules
    }

    /// Get the molecule at index `id`
    #[inline] pub fn molecule(&self, id: usize) -> &Molecule {
        &self.molecules[id]
    }

    /// Get the length of the shortest bond path to go from the particle `i` to
    /// the particle `j`. This length is 0 if there is no path from `i` to `j`,
    /// 1 if `i == j`, 2 if there is a bond between `i` and `j`, etc.
    pub fn shortest_path(&self, i: usize, j: usize) -> u8 {
        if !(self.are_in_same_molecule(i, j)) {
            0
        } else if i == j {
            1
        } else {
            let connect = self.molecule_containing(i).connectivity(i, j);
            if connect.contains(CONNECT_12) {
                2
            } else if connect.contains(CONNECT_13) {
                3
            } else if connect.contains(CONNECT_14) {
                4
            } else if connect.contains(CONNECT_FAR) {
                u8::MAX
            } else {
                unreachable!();
            }
        }
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

        let (i, j, perms) = if self.are_in_same_molecule(i, j) {
            (i, j, None)
        } else { // Do all the hard work
            // Getting copy of the molecules before the merge
            let mol_i = self.molecule_id(i);
            let mol_j = self.molecule_id(j);
            let new_mol_idx = min(mol_i, mol_j);
            let old_mol_idx = max(mol_i, mol_j);
            let new_mol = self.molecules[new_mol_idx].clone();
            let old_mol = self.molecules[old_mol_idx].clone();

            let mut perms = Permutations::new();
            let new_mol_last = new_mol.last();
            let old_mol_first = old_mol.first();

            // If new_mol_last + 1 == old_mol_first, no one move
            if new_mol_last + 1 != old_mol_first {
                let size = old_mol.size();
                let mut i = old_mol.first();
                // Add permutation for the molecule we just moved around
                for j in (new_mol_last + 1)..(new_mol_last + size + 1) {
                    perms.push((i, j));
                    i += 1;
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
                (i - delta, j) // i moved
            };
            (i, j, Some(perms))
        };

        if i == j {
            error!("Can not add a bond between a particle and itself.");
        }

        let id = self.molecule_id(i);
        self.molecules[id].add_bond(i, j);
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
        if part.kind == u16::max_value() {
            // If no value have been precised, set one from the internal list
            // of particles kinds.
            part.kind = self.get_kind(part.name());
        }
        self.particles.push(part);
        self.molecules.push(Molecule::new(self.particles.len() - 1));
    }

    /// Get the number of particles in this system
    #[inline] pub fn size(&self) -> usize {self.particles.len()}

    /// Get an iterator over the `Particle` in this system
    #[inline] pub fn iter(&self) -> slice::Iter<Particle> {
        self.particles.iter()
    }

    /// Get a mutable iterator over the `Particle` in this system
    #[inline] pub fn iter_mut(&mut self) -> slice::IterMut<Particle> {
        self.particles.iter_mut()
    }

    /// Get or create the usize kind index for the name `name` of a particle
    fn get_kind(&mut self, name: &str) -> u16 {
        if self.kinds.contains_key(name) {
            return self.kinds[name];
        } else {
            let index = self.kinds.len() as u16;
            self.kinds.insert(String::from(name), index);
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
        debug!("Merging molecules: {} <-- {}", new_mol_idx, old_mol_idx);
        trace!("The molecules contains:\n{:#?}\n ---\n{:#?}", new_mol, old_mol);

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

static NO_PAIR_INTERACTION: &'static [PairInteraction] = &[];
static NO_BOND_INTERACTION: &'static [Box<PairPotential>] = &[];
static NO_ANGLE_INTERACTION: &'static [Box<AnglePotential>] = &[];
static NO_DIHEDRAL_INTERACTION: &'static [Box<DihedralPotential>] = &[];

/// Potentials related functions
impl System {
    /// Get an helper struct to evaluate the energy of this system.
    pub fn energy_evaluator(&self) -> EnergyEvaluator {
        EnergyEvaluator::new(self)
    }

    /// Get the list of pair potential acting between the particles at indexes
    /// `i` and `j`.
    pub fn pair_potentials(&self, i: usize, j: usize) -> &[PairInteraction] {
        let ikind = self.particles[i].kind;
        let jkind = self.particles[j].kind;
        if let Some(val) = self.interactions.pairs(ikind, jkind) {
            &val
        } else {
            let i = self.particles[i].name();
            let j = self.particles[j].name();
            // TODO: add and use the warn_once! macro
            warn!("No potential defined for the pair ({}, {})", i, j);
            NO_PAIR_INTERACTION
        }
    }

    /// Get the list of bonded potential acting between the particles at indexes
    /// `i` and `j`.
    pub fn bond_potentials(&self, i: usize, j: usize) -> &[Box<PairPotential>] {
        let ikind = self.particles[i].kind;
        let jkind = self.particles[j].kind;
        if let Some(val) = self.interactions.bonds(ikind, jkind) {
            &val
        } else {
            let i = self.particles[i].name();
            let j = self.particles[j].name();
            // TODO: add and use the warn_once! macro
            warn!("No potential defined for the bond ({}, {})", i, j);
            NO_BOND_INTERACTION
        }
    }

    /// Get the list of angle interaction acting between the particles at
    /// indexes `i`, `j` and `k`.
    pub fn angle_potentials(&self, i: usize, j: usize, k: usize) -> &[Box<AnglePotential>] {
        let ikind = self.particles[i].kind;
        let jkind = self.particles[j].kind;
        let kkind = self.particles[k].kind;

        if let Some(val) =  self.interactions.angles(ikind, jkind, kkind) {
            &val
        } else {
            let i = self.particles[i].name();
            let j = self.particles[j].name();
            let k = self.particles[k].name();
            // TODO: add and use the warn_once! macro
            warn!("No potential defined for the angle ({}, {}, {})", i, j, k);
            NO_ANGLE_INTERACTION
        }
    }

    /// Get the list of dihedral angles interaction acting between the particles
    /// at indexes `i`, `j`, `k` and `m`.
    pub fn dihedral_potentials(&self, i: usize, j: usize, k: usize, m: usize) -> &[Box<DihedralPotential>] {
        let ikind = self.particles[i].kind;
        let jkind = self.particles[j].kind;
        let kkind = self.particles[k].kind;
        let mkind = self.particles[m].kind;

        if let Some(val) = self.interactions.dihedrals(ikind, jkind, kkind, mkind) {
            &val
        } else {
            let i = self.particles[i].name();
            let j = self.particles[j].name();
            let k = self.particles[k].name();
            let m = self.particles[m].name();
            // TODO: add and use the warn_once! macro
            warn!("No potential defined for the dihedral ({}, {}, {}, {})", i, j, k, m);
            NO_DIHEDRAL_INTERACTION
        }
    }

    /// Get the current coulombic solver
    pub fn coulomb_potential(&self) -> Option<&RefCell<Box<CoulombicPotential>>> {
        self.interactions.coulomb()
    }

    /// Get all the global potentials
    pub fn global_potentials(&self) -> &[RefCell<Box<GlobalPotential>>] {
        self.interactions.globals()
    }

    /// Add the `potential` pair potential between the particles with names
    /// `i` and `j`.
    pub fn add_pair_interaction(&mut self, i: &str, j: &str, potential: Box<PairPotential>) {
        let ikind = self.get_kind(i);
        let jkind = self.get_kind(j);

        self.interactions.add_pair(ikind, jkind, potential);
    }

    /// Add the `potential` pair potential between the particles with names
    /// `i` and `j`, using the `restriction` restriction scheme.
    pub fn add_pair_interaction_with_restriction(&mut self, i: &str, j: &str, potential: Box<PairPotential>, restriction: PairRestriction) {
        let ikind = self.get_kind(i);
        let jkind = self.get_kind(j);

        self.interactions.add_pair_with_restriction(ikind, jkind, potential, restriction);
    }

    /// Add the `potential` bonded potentials between the particles with names
    /// `i` and `j`
    pub fn add_bond_interaction(&mut self, i: &str, j: &str, potential: Box<PairPotential>) {
        let ikind = self.get_kind(i);
        let jkind = self.get_kind(j);

        self.interactions.add_bond(ikind, jkind, potential);
    }

    /// Add the `potential` angle potential between the particles with names `i`,
    /// `j`, and `k`
    pub fn add_angle_interaction(&mut self, i: &str, j: &str, k: &str, potential: Box<AnglePotential>) {
        let ikind = self.get_kind(i);
        let jkind = self.get_kind(j);
        let kkind = self.get_kind(k);

        self.interactions.add_angle(ikind, jkind, kkind, potential);
    }

    /// Add the `potential` dihedral angle potential between the particles with
    /// names names `i`, `j`, `k`, and `m`
    pub fn add_dihedral_interaction(&mut self, i: &str, j: &str, k: &str, m: &str, potential: Box<DihedralPotential>) {
        let ikind = self.get_kind(i);
        let jkind = self.get_kind(j);
        let kkind = self.get_kind(k);
        let mkind = self.get_kind(m);

        self.interactions.add_dihedral(ikind, jkind, kkind, mkind, potential);
    }

    /// Set the coulombic potential to `potential`
    pub fn set_coulomb_interaction(&mut self, potential: Box<CoulombicPotential>) {
        self.interactions.set_coulomb(potential);
    }

    /// Add a global interaction to the system
    pub fn add_global_interaction(&mut self, potential: Box<GlobalPotential>) {
        self.interactions.add_global(potential);
    }
}

impl<'a> IntoIterator for &'a System {
    type Item = &'a Particle;
    type IntoIter = slice::Iter<'a, Particle>;

    #[inline]
    fn into_iter(self) -> slice::Iter<'a, Particle> {
        self.iter()
    }
}

impl<'a> IntoIterator for &'a mut System {
    type Item = &'a mut Particle;
    type IntoIter = slice::IterMut<'a, Particle>;

    #[inline]
    fn into_iter(self) -> slice::IterMut<'a, Particle> {
        self.iter_mut()
    }
}

/// UnitCell related functions
impl System {
    /// Get a reference to  the system unit cell
    #[inline] pub fn cell(&self) -> &UnitCell {&self.cell}
    /// Get a mutable reference to  the system unit cell
    #[inline] pub fn cell_mut(&mut self) -> &mut UnitCell {&mut self.cell}
    /// Set the system unit cell
    #[inline] pub fn set_cell(&mut self, cell: UnitCell) {self.cell = cell;}

    /// Get the distance between the particles at indexes `i` and `j`
    #[inline] pub fn distance(&self, i: usize, j:usize) -> f64 {
        self.cell.distance(&self.particles[i].position, &self.particles[j].position)
    }

    /// Get the vector ri - rj wrapped in the cell.
    pub fn wraped_vector(&self, i: usize, j:usize) -> Vector3D {
        let mut res = self.particles[i].position - self.particles[j].position;
        self.cell.wrap_vector(&mut res);
        return res;
    }

    /// Get the angle between the particles `i`, `j` and `k`
    pub fn angle(&self, i: usize, j: usize, k: usize) -> f64 {
        let a = self.particles[i].position;
        let b = self.particles[j].position;
        let c = self.particles[k].position;
        self.cell.angle(&a, &b, &c)
    }

    /// Get the angle and the derivatives of the angle between the particles
    /// `i`, `j` and `k`
    pub fn angle_and_derivatives(&self, i: usize, j: usize, k: usize) -> (f64, Vector3D, Vector3D, Vector3D) {
        let a = self.particles[i].position;
        let b = self.particles[j].position;
        let c = self.particles[k].position;
        self.cell.angle_and_derivatives(&a, &b, &c)
    }

    /// Get the dihedral angle between the particles `i`, `j`, `k` and `m`
    pub fn dihedral(&self, i: usize, j: usize, k: usize, m: usize) -> f64 {
        let a = self.particles[i].position;
        let b = self.particles[j].position;
        let c = self.particles[k].position;
        let d = self.particles[m].position;
        self.cell.dihedral(&a, &b, &c, &d)
    }

    /// Get the dihedral angle and the derivatives of the dihedral angle
    /// between the particles `i`, `j`, `k` and `m`
    pub fn dihedral_and_derivatives(&self, i: usize, j: usize, k: usize, m: usize) -> (f64, Vector3D, Vector3D, Vector3D, Vector3D) {
        let a = self.particles[i].position;
        let b = self.particles[j].position;
        let c = self.particles[k].position;
        let d = self.particles[m].position;
        self.cell.dihedral_and_derivatives(&a, &b, &c, &d)
    }
}

/******************************************************************************/
use system::Compute;
use system::{PotentialEnergy, KineticEnergy, TotalEnergy};
use system::Forces;
use system::Temperature;
use system::Volume;
use system::{Virial, Stress, Pressure};
use system::{StressAtTemperature, PressureAtTemperature};

/// Functions to get pysical properties of a system.
impl System {
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
    /// Get the pressure of the system from the virial equation, at the system
    /// instananeous temperature.
    pub fn pressure(&self) -> f64 {Pressure.compute(self)}
    /// Get the stress tensor of the system from the virial equation, at the
    /// system instananeous temperature.
    pub fn stress(&self) -> Matrix3 {Stress.compute(self)}

    /// Get the pressure of the system from the virial equation, at the given
    /// `temperature`.
    pub fn pressure_at_temperature(&self, temperature: f64) -> f64 {
        PressureAtTemperature{temperature: temperature}.compute(self)
    }
    /// Get the stress tensor of the system from the virial equation, at the
    /// given `temperature`.
    pub fn stress_at_temperature(&self, temperature: f64) -> Matrix3 {
        StressAtTemperature{temperature: temperature}.compute(self)
    }

    /// Get the forces acting on all the particles in the system
    pub fn forces(&self) -> Vec<Vector3D> {
        Forces.compute(self)
    }
}

/******************************************************************************/
impl Index<usize> for System {
    type Output = Particle;
    #[inline]
    fn index(&self, index: usize) -> &Particle {
        &self.particles[index]
    }
}

impl IndexMut<usize> for System {
    #[inline]
    fn index_mut(&mut self, index: usize) -> &mut Particle {
        &mut self.particles[index]
    }
}

#[cfg(test)]
mod tests {
    use system::*;
    use types::*;
    use potentials::*;

    #[test]
    fn step() {
        let mut system = System::new();

        assert_eq!(system.step(), 0);

        system.increment_step();
        system.increment_step();
        system.increment_step();

        assert_eq!(system.step(), 3);
    }

    #[test]
    fn cell() {
        let mut system = System::from_cell(UnitCell::cubic(67.0));
        {
            let cell = system.cell();
            assert_eq!(cell.a(), 67.0);
        }

        system.set_cell(UnitCell::cubic(10.0));
        let cell = system.cell();
        assert_eq!(cell.a(), 10.0);
    }

    #[test]
    fn molecules() {
        let mut system = System::new();
        system.add_particle(Particle::new("H"));
        system.add_particle(Particle::new("O"));
        system.add_particle(Particle::new("O"));
        system.add_particle(Particle::new("H"));

        assert_eq!(system.add_bond(0, 1), Some(vec![]));
        assert_eq!(system.add_bond(2, 3), Some(vec![]));

        assert_eq!(system.molecules().len(), 2);

        let molecule = system.molecule(0).clone();
        assert!(molecule.bonds().contains(&Bond::new(0, 1)));
        let molecule = system.molecule(1).clone();
        assert!(molecule.bonds().contains(&Bond::new(2, 3)));

        assert_eq!(system.add_bond(1, 2), Some(vec![]));
        assert_eq!(system.molecules().len(), 1);

        let molecule = system.molecule(0).clone();
        assert!(molecule.angles().contains(&Angle::new(0, 1, 2)));
        assert!(molecule.angles().contains(&Angle::new(1, 2, 3)));
        assert!(molecule.dihedrals().contains(&Dihedral::new(0, 1, 2, 3)));

        system.remove_particle(2);
        assert_eq!(system.molecules().len(), 1);

        system.remove_molecule_containing(1);
        assert_eq!(system.molecules().len(), 0);
        assert_eq!(system.size(), 0);
    }

    #[test]
    fn shortest_path() {
        let mut system = System::new();
        system.add_particle(Particle::new("C"));
        system.add_particle(Particle::new("C"));
        system.add_particle(Particle::new("C"));
        system.add_particle(Particle::new("C"));
        system.add_particle(Particle::new("C"));
        system.add_particle(Particle::new("Zn"));

        system.add_bond(0, 1);
        system.add_bond(1, 2);
        system.add_bond(2, 3);
        system.add_bond(3, 4);

        assert_eq!(system.shortest_path(0, 0), 1);
        assert_eq!(system.shortest_path(0, 1), 2);
        assert_eq!(system.shortest_path(0, 2), 3);
        assert_eq!(system.shortest_path(0, 3), 4);
        assert!(system.shortest_path(0, 4) > 4);
        assert_eq!(system.shortest_path(0, 5), 0);
    }

    #[test]
    fn molecules_permutations() {
        let mut system = System::new();
        system.add_particle(Particle::new("N"));
        system.add_particle(Particle::new("H"));
        system.add_particle(Particle::new("H"));
        system.add_particle(Particle::new("H"));
        system.add_particle(Particle::new("N"));
        system.add_particle(Particle::new("H"));
        system.add_particle(Particle::new("H"));
        system.add_particle(Particle::new("H"));

        assert_eq!(system.add_bond(0, 3), Some(vec![(3, 1), (1, 2), (2, 3)]));
        assert_eq!(system.add_bond(0, 3), Some(vec![(3, 2), (2, 3)]));
        assert_eq!(system.add_bond(0, 3), Some(vec![]));

        assert_eq!(system.add_bond(4, 5), Some(vec![]));
        assert_eq!(system.add_bond(4, 7), Some(vec![(7, 6), (6, 7)]));
        assert_eq!(system.add_bond(4, 7), Some(vec![]));
    }

    #[test]
    fn particles() {
        let mut system = System::new();
        system.add_particle(Particle::new("O"));
        system.add_particle(Particle::new("H"));
        system.add_particle(Particle::new("H"));

        assert_eq!(system.size(), 3);
        assert_eq!(system[0].name(), "O");
        assert_eq!(system[1].name(), "H");
        assert_eq!(system[2].name(), "H");
    }

    #[test]
    fn distances() {
        let mut system = System::from_cell(UnitCell::cubic(5.0));
        system.add_particle(Particle::new("O"));
        system.add_particle(Particle::new("H"));

        system[0].position = Vector3D::new(9.0, 0.0, 0.0);
        system[1].position = Vector3D::zero();
        assert_eq!(system.distance(0, 1), 1.0);

        system.set_cell(UnitCell::new());
        assert_eq!(system.distance(0, 1), 9.0);
    }

    #[test]
    fn interactions() {
        let mut system = System::new();
        system.add_particle(Particle::new("He"));

        system.add_pair_interaction("He", "He", Box::new(LennardJones{sigma: 0.3, epsilon: 2.0}));
        system.add_pair_interaction("He", "He", Box::new(Harmonic{k: 100.0, x0: 1.1}));
        assert_eq!(system.pair_potentials(0, 0).len(), 2);

        system.add_bond_interaction("He", "He", Box::new(Harmonic{k: 100.0, x0: 1.1}));
        assert_eq!(system.bond_potentials(0, 0).len(), 1);

        system.add_angle_interaction("He", "He", "He", Box::new(Harmonic{k: 100.0, x0: 1.1}));
        assert_eq!(system.angle_potentials(0, 0, 0).len(), 1);

        system.add_dihedral_interaction("He", "He", "He", "He", Box::new(CosineHarmonic::new(0.3, 2.0)));
        assert_eq!(system.dihedral_potentials(0, 0, 0, 0).len(), 1);

        assert!(system.coulomb_potential().is_none());
        system.set_coulomb_interaction(Box::new(Wolf::new(1.0)));
        assert!(system.coulomb_potential().is_some());

        system.add_global_interaction(Box::new(Wolf::new(1.0)));
        assert_eq!(system.global_potentials().len(), 1);
    }

    #[test]
    fn missing_interaction() {
        let mut system = System::new();
        system.add_particle(Particle::new("He"));
        assert_eq!(system.pair_potentials(0, 0).len(), 0);
        assert_eq!(system.bond_potentials(0, 0).len(), 0);
        assert_eq!(system.angle_potentials(0, 0, 0).len(), 0);
        assert_eq!(system.dihedral_potentials(0, 0, 0, 0).len(), 0);
    }

    #[test]
    fn moltype() {
        let mut system = System::new();
        // One helium
        system.add_particle(Particle::new("He"));
        // Two water molecules
        system.add_particle(Particle::new("H"));
        system.add_particle(Particle::new("O"));
        system.add_particle(Particle::new("H"));
        system.add_bond(1, 2);
        system.add_bond(2, 3);
        system.add_particle(Particle::new("H"));
        system.add_particle(Particle::new("O"));
        system.add_particle(Particle::new("H"));
        system.add_bond(4, 5);
        system.add_bond(5, 6);
        // Another helium
        system.add_particle(Particle::new("He"));
        // A water molecules, with different atoms order
        system.add_particle(Particle::new("O"));
        system.add_particle(Particle::new("H"));
        system.add_particle(Particle::new("H"));
        system.add_bond(8, 9);
        system.add_bond(8, 10);

        assert_eq!(system.molecules().len(), 5);
        // The heliums
        assert_eq!(system.moltype(0), system.moltype(3));

        // The waters
        assert!(system.moltype(1) == system.moltype(2));
        assert!(system.moltype(1) != system.moltype(4));
    }
}
