// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

//! The base type for simulation data in `lumol` is the `System` type.
//!
//! An `System` consists of a list of `Particle`; a list of `Molecule`
//! specifying how the particles are bonded together; an unit cell for boundary
//! conditions; and the interactions between these particles.
use std::ops::{Index, IndexMut, Range};
use std::slice;
use std::cmp::{min, max};
use std::iter::IntoIterator;
use std::i8;
use std::collections::BTreeMap;

use energy::PairInteraction;
use energy::{BondPotential, AnglePotential, DihedralPotential};
use types::{Vector3D, Matrix3, Zero};

use super::{Particle, ParticleKind};
use super::Molecule;
use super::{CONNECT_12, CONNECT_13, CONNECT_14, CONNECT_FAR};
use super::UnitCell;
use super::interactions::Interactions;
use super::EnergyEvaluator;
use super::molecules::molecule_type;

/// Particles permutations:. Indexes are given in the `(old, new)` form.
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
/// In the implementation, the particles contained in a molecule are guaranteed
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
    /// Molecules indexes for all the particles
    molids: Vec<usize>,
    /// Interactions manages the associations between particles and potentials
    interactions: Interactions,
    /// Current step of the simulation
    step: u64,
    /// Externally managed temperature for the system, if the propagation
    /// algorithm does not update the velocities.
    external_temperature: Option<f64>
}

impl Default for System {
    fn default() -> System {
        System::new()
    }
}

impl System {
    /// Create a new empty System
    pub fn new() -> System {
        System {
            particles: Vec::new(),
            molecules: Vec::new(),
            molids: Vec::new(),
            interactions: Interactions::new(),
            cell: UnitCell::new(),
            step: 0,
            external_temperature: None
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
    /// Get the type of the molecule at index `molid`. This type is a hash of
    /// the atoms names, and the set of bonds in the molecule. This means that
    /// two molecules will have the same type if and only if they contains the
    /// same atoms and the same bonds, **in the same order**.
    pub fn molecule_type(&self, molid: usize) -> u64 {
        let molecule = self.molecule(molid);
        molecule_type(molecule, &self.particles[molecule.into_iter()])
    }

    /// Get a list of molecules with `moltype` molecule type.
    pub fn molecules_with_moltype(&self, moltype: u64) -> Vec<usize> {
        let mut res = Vec::new();
        for i in 0..self.molecules().len() {
            if self.molecule_type(i) == moltype {
                res.push(i);
            }
        }
        return res;
    }

    /// Check if the particles at indexes `i` and `j` are in the same molecule
    #[inline] pub fn are_in_same_molecule(&self, i: usize, j:usize) -> bool {
        debug_assert!(self.molids.len() == self.particles.len());
        self.molids[i] == self.molids[j]
    }

    /// Get the list of molecules in the system.
    #[inline] pub fn molecules(&self) -> &[Molecule] {
        &self.molecules
    }

    /// Get the list of molecules in the system.
    #[inline] pub fn molecules_mut(&mut self) -> &[Molecule] {
        &self.molecules
    }

    /// Get the molecule at index `id`
    #[inline] pub fn molecule(&self, id: usize) -> &Molecule {
        &self.molecules[id]
    }

    /// Get the index of the molecule containing the particle `i`
    #[inline] pub fn molid(&self, i: usize) -> usize {
        self.molids[i]
    }

    /// Get the particles of the molecule
    #[inline] pub fn particles_of_molecule(&self, id: usize) -> &[Particle] {
        let range = self.molecules[id].first() .. self.molecules[id].size();
        &self.particles[range]
    }

        /// Get the particles of the molecule
    #[inline] pub fn particles_of_molecule2(&self, molecule: &Molecule) 
        -> &[Particle] {
        let range = molecule.first() .. molecule.size();
        &self.particles[range]
    }

    /// Get the length of the shortest bond path to go from the particle `i` to
    /// the particle `j`. If the particles are not in the same molecule, the
    /// length is -1. Else, this length is 0 if `i == j`, 1 if there is a bond
    /// between `i` and `j`, etc.
    pub fn bond_distance(&self, i: usize, j: usize) -> i8 {
        assert!(i < self.size() && j < self.size());
        if !(self.are_in_same_molecule(i, j)) {
            -1
        } else if i == j {
            0
        } else {
            let connect = self.molecule(self.molid(i)).connectivity(i, j);
            if connect.contains(CONNECT_12) {
                1
            } else if connect.contains(CONNECT_13) {
                2
            } else if connect.contains(CONNECT_14) {
                3
            } else if connect.contains(CONNECT_FAR) {
                i8::MAX
            } else {
                unreachable!();
            }
        }
    }

    /// Remove the molecule at index `i`
    pub fn remove_molecule(&mut self, molid: usize) {
        let molecule = self.molecules.remove(molid);
        let first = molecule.start();
        let size = molecule.size();

        for _ in 0..size {
            let _ = self.particles.remove(first);
            let _ = self.molids.remove(first);
        }

        for molecule in self.molecules.iter_mut().skip(molid) {
            molecule.translate_by(-(size as isize));
        }

        for molid in self.molids.iter_mut().skip(first) {
            *molid -= 1;
        }
    }

    /// Add a bond between the particles at indexes `i` and `j`. The particles
    /// should have been added to the system before calling this.
    ///
    /// # Warning
    ///
    /// If the bond is between two particles which are not in the same molecule,
    /// the two molecules are merged together by moving particles in the
    /// particles list, and thus invalidate any previously stored index. In
    /// particular, any bond, angle, dihedral or molecule is invalidated.
    ///
    /// This function will return the list of atomic permutations that where
    /// applied in order to ensure that molecules are contiguous in memory.
    pub fn add_bond(&mut self, mut particle_i: usize, mut particle_j: usize) -> Permutations {
        assert!(particle_i <= self.particles.len());
        assert!(particle_j <= self.particles.len());
        trace!(
            "Adding bond {} --- {}, in molecules {} and {}",
            particle_i, particle_j, self.molids[particle_i], self.molids[particle_j]
        );

        // Getting copy of the molecules before the merge
        let molid_i = self.molids[particle_i];
        let molid_j = self.molids[particle_j];
        let new_molid = min(molid_i, molid_j);
        let old_molid = max(molid_i, molid_j);
        let new_mol = self.molecules[new_molid].clone();
        let old_mol = self.molecules[old_molid].clone();

        // Effective merge
        let delta = self.merge_molecules(molid_i, molid_j);

        let mut permutations = Permutations::new();
        // If new_mol.last() + 1 == old_mol.first(), no one moved. Else,
        // we generate the permutations
        if !self.are_in_same_molecule(particle_i, particle_j) && new_mol.end() != old_mol.start() {
            let size = old_mol.size();
            let first = old_mol.start();
            let second = new_mol.end();
            // Add permutation for the molecule we just moved around
            for i in 0..size {
                permutations.push((first + i, second + i));
            }

            // Add permutations for molecules that where shifted to make
            // space for the just moved molecule.
            for molecule in &self.molecules[new_molid + 1 .. old_molid] {
                for i in molecule {
                    permutations.push((i - size, i));
                }
            }
        }

        // One of the `particle_i` or `particle_j` index is no longer valid, as one
        // molecule has been displaced.
        if molid_i == new_molid {
            particle_j -= delta; // j moved
        } else {
            particle_i -= delta; // i moved
        };

        assert!(self.molids[particle_i] == self.molids[particle_j]);
        self.molecules[self.molids[particle_i]].add_bond(particle_i, particle_j);
        return permutations;
    }

    /// Removes particle at index `i` and any associated bonds, angle or dihedral
    pub fn remove_particle(&mut self, i: usize) {
        let id = self.molids[i];
        self.molecules[id].remove_particle(i);

        for molecule in self.molecules.iter_mut().skip(id + 1) {
            molecule.translate_by(-1);
        }

        let _ = self.particles.remove(i);
        let _ = self.molids.remove(i);
    }

    /// Insert a particle at the end of the internal list
    pub fn add_particle(&mut self, p: Particle) {
        let mut part = p;
        if part.kind == ParticleKind::default() {
            // If no value have been precised, set one from the internal list
            // of particles kinds.
            part.kind = self.interactions.get_kind(part.name());
        }
        self.particles.push(part);
        self.molecules.push(Molecule::new(self.particles.len() - 1));
        self.molids.push(self.molecules.len() - 1);
    }

    /// Get the number of particles in this system
    #[inline] pub fn size(&self) -> usize {self.particles.len()}

    /// Return the center-of-mass of a molecule
    pub fn molecule_com(&self, molid: usize) -> Vector3D {
        // iterate over all particles of molecule(molid)
        let total_mass = self.molecule(molid)
            .iter()
            .fold(0.0, |total_mass, pi| total_mass + self[pi].mass);
        let com = self.molecule(molid)
            .iter()
            .fold(Vector3D::zero(), |com , pi| {
                com + self[pi].mass * self[pi].position
            });
        com / total_mass
    }

    /// Return the center-of-mass of the system
    pub fn center_of_mass(&self) -> Vector3D {
        // iterate over all particles in the system
        let total_mass = self
            .iter()
            .fold(0.0, |total_mass, particle| total_mass + particle.mass);
        let com: Vector3D = self
            .iter()
            .fold(Vector3D::zero(), |com, particle| {
                com + particle.position * particle.mass
            });
        com / total_mass
    }

    /// Get an iterator over the `Particle` in this system
    #[inline] pub fn iter(&self) -> slice::Iter<Particle> {
        self.particles.iter()
    }

    /// Get a mutable iterator over the `Particle` in this system
    #[inline] pub fn iter_mut(&mut self) -> slice::IterMut<Particle> {
        self.particles.iter_mut()
    }

    /// Merge the molecules at indexes `first` and `second` into one
    /// molecule. The molecule are merged into the one with the lower molecule
    /// index.
    ///
    /// For example, if we have
    /// ```
    ///  0    1    2    3   # Molecules indexes
    /// H-H  H-H  H-H  H-H
    /// 0 1  2 3  4 5  6 7  # Particles indexes
    /// ```
    /// and call `merge_molecules(0, 3)` the result will be
    /// ```
    /// 0 1 6 7  2 3  4 5  # Old indexes
    /// H-H-H-H  H-H  H-H
    /// 0 1 2 3  4 5  6 7  # New indexes
    /// ```
    ///
    /// This functions return the change in index of the first particle of the
    /// moved molecule, i.e. in this example `4`.
    #[allow(block_in_if_condition_stmt)]
    fn merge_molecules(&mut self, first: usize, second: usize) -> usize {
        let (new_molid, old_molid) = if first < second {
            (first, second)
        } else {
            (second, first)
        };

        let mut new_mol = self.molecules[new_molid].clone();
        let old_mol = self.molecules[old_molid].clone();

        if new_mol.end() == old_mol.start() {
            // Just update the molecules ids
            for i in old_mol {
                self.molids[i] = new_molid;
            }
        } else {
            // Move the particles close together
            let mut new_index = new_mol.end();
            for i in old_mol {
                // Remove particles from the old position, and insert it to the
                // new one. The indexes are valid during the movement, because
                // we insert a new particle for each particle removed.
                let particle = self.particles.remove(i);
                self.particles.insert(new_index, particle);

                // Update molids
                let _ = self.molids.remove(i);
                self.molids.insert(new_index, new_molid);

                new_index += 1;
            }
        }

        let mut old_mol = self.molecules[old_molid].clone();
        let size = old_mol.size() as isize;

        // translate all indexes in the molecules between new_mol and old_mol
        for molecule in self.molecules.iter_mut().skip(new_molid + 1).take(old_molid - new_molid - 1) {
            molecule.translate_by(size);
        }

        // Update molid for all particles after the old molecule
        for molid in self.molids.iter_mut().skip(old_mol.end()) {
            *molid -= 1;
        }

        let delta = old_mol.start() - new_mol.end();
        old_mol.translate_by(- (delta as isize));

        new_mol.merge_with(old_mol);
        self.molecules[new_molid] = new_mol;
        let _ = self.molecules.remove(old_molid);

        // Check that self.molids is sorted and only contains successive values
        debug_assert!(self.molids.iter().fold((true, 0), |(is_valid, previous), &i| {
            if i == previous || i == previous + 1 {
                (is_valid, i)
            } else {
                (false, i)
            }
        }).0, "Unsorted molecule ids {:?}", self.molids);

        return delta;
    }

    /// Get the number of particles of each kind in the system
    pub fn composition(&self) -> BTreeMap<ParticleKind, usize> {
        let mut map = BTreeMap::new();
        for particle in &self.particles {
            *map.entry(particle.kind).or_insert(0) += 1;
        }
        return map;
    }

    /// Get a list of all the particles kinds in the system.
    pub fn particle_kinds(&self) -> Vec<ParticleKind> {
        self.interactions.all_kinds()
    }
}

/// Potentials related functions
impl System {
    /// Get an helper struct to evaluate the energy of this system.
    pub fn energy_evaluator(&self) -> EnergyEvaluator {
        EnergyEvaluator::new(self)
    }

    /// Access the interactions for this system
    pub fn interactions(&self) -> &Interactions {
        &self.interactions
    }

    /// Access the interactions for this system in a mutable way
    pub fn interactions_mut(&mut self) -> &mut Interactions {
        &mut self.interactions
    }

    /// Get the list of pair potential acting between the particles at indexes
    /// `i` and `j`.
    pub fn pair_potentials(&self, i: usize, j: usize) -> &[PairInteraction] {
        let ikind = self.particles[i].kind;
        let jkind = self.particles[j].kind;
        self.interactions.pairs(ikind, jkind)
    }

    /// Get the list of bonded potential acting between the particles at indexes
    /// `i` and `j`.
    pub fn bond_potentials(&self, i: usize, j: usize) -> &[Box<BondPotential>] {
        let ikind = self.particles[i].kind;
        let jkind = self.particles[j].kind;
        self.interactions.bonds(ikind, jkind)
    }

    /// Get the list of angle interaction acting between the particles at
    /// indexes `i`, `j` and `k`.
    pub fn angle_potentials(&self, i: usize, j: usize, k: usize) -> &[Box<AnglePotential>] {
        let ikind = self.particles[i].kind;
        let jkind = self.particles[j].kind;
        let kkind = self.particles[k].kind;
        self.interactions.angles(ikind, jkind, kkind)
    }

    /// Get the list of dihedral angles interaction acting between the particles
    /// at indexes `i`, `j`, `k` and `m`.
    pub fn dihedral_potentials(&self, i: usize, j: usize, k: usize, m: usize) -> &[Box<DihedralPotential>] {
        let ikind = self.particles[i].kind;
        let jkind = self.particles[j].kind;
        let kkind = self.particles[k].kind;
        let mkind = self.particles[m].kind;
        self.interactions.dihedrals(ikind, jkind, kkind, mkind)
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

/// `UnitCell` related functions
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

    /// Get the vector between the nearest image of particle `j` with respect to
    /// particle `i`.
    pub fn nearest_image(&self, i: usize, j:usize) -> Vector3D {
        let mut res = self.particles[i].position - self.particles[j].position;
        self.cell.vector_image(&mut res);
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
use sys::compute::Compute;
use sys::compute::{PotentialEnergy, KineticEnergy, TotalEnergy};
use sys::compute::Forces;
use sys::compute::Temperature;
use sys::compute::Volume;
use sys::compute::{Virial, Stress, Pressure};
use sys::compute::{StressAtTemperature, PressureAtTemperature};

/// Functions to get physical properties of a system.
impl System {
    /// Get the kinetic energy of the system.
    pub fn kinetic_energy(&self) -> f64 {KineticEnergy.compute(self)}
    /// Get the potential energy of the system.
    pub fn potential_energy(&self) -> f64 {PotentialEnergy.compute(self)}
    /// Get the total energy of the system.
    pub fn total_energy(&self) -> f64 {TotalEnergy.compute(self)}

    /// Use an external temperature for all the system properties. Calling this
    /// with `Some(temperature)` will replace all the computation of the
    /// temperature from the velocities with the given values. Calling it with
    /// `None` will use the velocities.
    ///
    /// The default is to use the velocities unless this function is called.
    pub fn external_temperature(&mut self, temperature: Option<f64>) {
        if let Some(temperature) = temperature {
            assert!(temperature >= 0.0, "External temperature must be positive");
        }
        self.external_temperature = temperature;
    }

    /// Get the temperature of the system.
    pub fn temperature(&self) -> f64 {
        match self.external_temperature {
            Some(value) => value,
            None => Temperature.compute(self)
        }
    }

    /// Get the volume of the system.
    pub fn volume(&self) -> f64 {Volume.compute(self)}

    /// Get the virial of the system as a tensor
    pub fn virial(&self) -> Matrix3 {Virial.compute(self)}
    /// Get the pressure of the system from the virial equation, at the system
    /// instantaneous temperature.
    pub fn pressure(&self) -> f64 {
        match self.external_temperature {
            Some(temperature) => {
                PressureAtTemperature{temperature: temperature}.compute(self)
            }
            None => Pressure.compute(self)
        }
    }
    /// Get the stress tensor of the system from the virial equation.
    pub fn stress(&self) -> Matrix3 {
        match self.external_temperature {
            Some(temperature) => {
                StressAtTemperature{temperature: temperature}.compute(self)
            }
            None => Stress.compute(self)
        }

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

impl Index<Range<usize>> for System {
    type Output = [Particle];
    #[inline]
    fn index(&self, index: Range<usize>) -> &[Particle] {
        // Index::index(&self.particles, index)
        &self.particles[index]
    }
}

#[cfg(test)]
mod tests {
    use sys::*;
    use types::*;

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

        assert_eq!(system.add_bond(0, 1), vec![]);
        assert_eq!(system.add_bond(2, 3), vec![]);

        assert_eq!(system.molecules().len(), 2);

        let molecule = system.molecule(0).clone();
        assert!(molecule.bonds().contains(&Bond::new(0, 1)));
        let molecule = system.molecule(1).clone();
        assert!(molecule.bonds().contains(&Bond::new(2, 3)));

        assert_eq!(system.add_bond(1, 2), vec![]);
        assert_eq!(system.molecules().len(), 1);

        let molecule = system.molecule(0).clone();
        assert!(molecule.angles().contains(&Angle::new(0, 1, 2)));
        assert!(molecule.angles().contains(&Angle::new(1, 2, 3)));
        assert!(molecule.dihedrals().contains(&Dihedral::new(0, 1, 2, 3)));

        system.remove_particle(2);
        assert_eq!(system.molecules().len(), 1);

        let molid = system.molid(1);
        system.remove_molecule(molid);
        assert_eq!(system.molecules().len(), 0);
        assert_eq!(system.size(), 0);
    }

    #[test]
    fn add_bonds() {
        // This is a regression test for issue #76
        let mut system = System::new();
        system.add_particle(Particle::new("H"));
        system.add_particle(Particle::new("H"));
        system.add_particle(Particle::new("O"));

        assert_eq!(system.add_bond(0, 2), vec![(2, 1), (1, 2)]);
        assert_eq!(system.add_bond(2, 1), vec![]);
        assert_eq!(system.molecules().len(), 1);
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

        let _ = system.add_bond(0, 1);
        let _ = system.add_bond(1, 2);
        let _ = system.add_bond(2, 3);
        let _ = system.add_bond(3, 4);

        assert_eq!(system.bond_distance(0, 0), 0);
        assert_eq!(system.bond_distance(0, 1), 1);
        assert_eq!(system.bond_distance(0, 2), 2);
        assert_eq!(system.bond_distance(0, 3), 3);
        assert!(system.bond_distance(0, 4) > 3);
        assert_eq!(system.bond_distance(0, 5), -1);
    }

    #[test]
    fn molecules_permutations() {
        let mut system = System::new();
        system.add_particle(Particle::new("C"));
        system.add_particle(Particle::new("H"));
        system.add_particle(Particle::new("H"));
        system.add_particle(Particle::new("H"));

        system.add_particle(Particle::new("C"));
        system.add_particle(Particle::new("H"));
        system.add_particle(Particle::new("H"));
        system.add_particle(Particle::new("H"));

        assert_eq!(system.add_bond(0, 3), vec![(3, 1), (1, 2), (2, 3)]);
        assert_eq!(system.add_bond(0, 3), vec![(3, 2), (2, 3)]);
        assert_eq!(system.add_bond(0, 3), vec![]);

        assert_eq!(system.add_bond(4, 5), vec![]);
        assert_eq!(system.add_bond(4, 7), vec![(7, 6), (6, 7)]);
        assert_eq!(system.add_bond(4, 7), vec![]);
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
    fn particles_from_range() {
        let mut system = System::new();
        let particles = [Particle::new("O"),
                         Particle::new("H"),
                         Particle::new("C")];

        system.add_particle(particles[0].clone());
        system.add_particle(particles[1].clone());
        system.add_particle(particles[2].clone());
        {
            let particles_range = &system[0..2];
            assert_eq!(particles_range[0].name(), particles[0].name());
            assert_eq!(particles_range[1].name(), particles[1].name());    
        }
        {
            let particles_range = &system[1..3];
            assert_eq!(particles_range[0].name(), particles[1].name());
            assert_eq!(particles_range[1].name(), particles[2].name());    
        }
        {
            let particles_range = &system[1..3];
            assert_eq!(particles_range[0].name(), particles[1].name());
            assert_eq!(particles_range[1].name(), particles[2].name());
        }
    }

    #[test]
    fn particle_kinds() {
        let mut system = System::new();
        system.add_particle(Particle::new("H"));
        system.add_particle(Particle::new("O"));
        system.add_particle(Particle::new("O"));
        system.add_particle(Particle::new("H"));
        system.add_particle(Particle::new("C"));
        system.add_particle(Particle::new("U"));

        let kinds = system.particle_kinds();
        assert_eq!(kinds.len(), 4);

        assert!(kinds.contains(&ParticleKind(0)));
        assert!(kinds.contains(&ParticleKind(1)));
        assert!(kinds.contains(&ParticleKind(2)));
        assert!(kinds.contains(&ParticleKind(3)));
    }

    #[test]
    fn composition() {
        let mut system = System::new();
        system.add_particle(Particle::new("H"));
        system.add_particle(Particle::new("O"));
        system.add_particle(Particle::new("O"));
        system.add_particle(Particle::new("H"));
        system.add_particle(Particle::new("C"));
        system.add_particle(Particle::new("U"));
        system.add_particle(Particle::new("H"));

        let compo = system.composition();
        assert_eq!(compo.len(), 4);
        assert_eq!(compo[&ParticleKind(0)], 3);
        assert_eq!(compo[&ParticleKind(1)], 2);
        assert_eq!(compo[&ParticleKind(2)], 1);
        assert_eq!(compo[&ParticleKind(3)], 1);
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
    fn center_of_mass() {
        let mut system = System::from_cell(UnitCell::cubic(5.0));
        system.add_particle(Particle::new("O"));
        system.add_particle(Particle::new("O"));
        let _ = system.add_bond(0, 1);
        system[0].position = Vector3D::new(9.0, 0.0, 0.0);
        system[1].position = Vector3D::zero();
        assert_eq!(system.molecule_com(0), Vector3D::new(4.5, 0.0, 0.0));
        assert_eq!(system.center_of_mass(), Vector3D::new(4.5, 0.0, 0.0));
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
        let _ = system.add_bond(1, 2);
        let _ = system.add_bond(2, 3);
        system.add_particle(Particle::new("H"));
        system.add_particle(Particle::new("O"));
        system.add_particle(Particle::new("H"));
        let _ = system.add_bond(4, 5);
        let _ = system.add_bond(5, 6);
        // Another helium
        system.add_particle(Particle::new("He"));
        // A water molecules, with different atoms order
        system.add_particle(Particle::new("O"));
        system.add_particle(Particle::new("H"));
        system.add_particle(Particle::new("H"));
        let _ = system.add_bond(8, 9);
        let _ = system.add_bond(8, 10);

        assert_eq!(system.molecules().len(), 5);
        // The helium particles
        assert_eq!(system.molecule_type(0), system.molecule_type(3));

        // The water molecules
        assert!(system.molecule_type(1) == system.molecule_type(2));
        assert!(system.molecule_type(1) != system.molecule_type(4));
    }

    #[test]
    #[should_panic]
    fn negative_external_temperature() {
        let mut system = System::new();
        system.external_temperature(Some(-1.0));
    }
}
