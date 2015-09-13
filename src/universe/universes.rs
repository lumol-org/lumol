/*
 * Cymbalum, Molecular Simulation in Rust
 * Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
*/

//! `Universe` type definition and implementation.

use std::collections::HashMap;
use std::ops::{Index, IndexMut};
use std::slice;

use std::io::prelude::*;
use std::io;
use std::fs::File;
use std::num;

use ::potentials::{PairPotential, AnglePotential, DihedralPotential};
use ::types::{Vector3D, Matrix3};

use super::Particle;
use super::Topology;
use super::UnitCell;
use super::interactions::Interactions;

/// The Universe type hold all the data about a system. This data contains:
///
///   - an unit cell, containing the system;
///   - a list of particles in the system;
///   - a list of interactions, associating particles kinds and potentials
///   - a hash map associating particles names and particles kinds.
pub struct Universe {
    /// Unit cell of the universe
    cell: UnitCell,
    /// List of particles in the system
    particles: Vec<Particle>,
    /// Topology of the universe
    topology: Topology,
    /// Particles kinds, associating particles names and indexes
    kinds: HashMap<String, u16>,
    /// Interactions is a hash map associating particles kinds and potentials
    interactions: Interactions,
}

/// Possible error causes where reading XYZ file
#[derive(Debug)]
pub enum ReadError {
    IoError(io::Error),
    ParseIntError(num::ParseIntError),
    ParseFloatError(num::ParseFloatError),
    XYZFormatError{err: &'static str}, // Badly formated file
}

impl From<num::ParseIntError> for ReadError {
    fn from(err: num::ParseIntError) -> ReadError {
        ReadError::ParseIntError(err)
    }
}

impl From<num::ParseFloatError> for ReadError {
    fn from(err: num::ParseFloatError) -> ReadError {
        ReadError::ParseFloatError(err)
    }
}

impl From<io::Error> for ReadError {
    fn from(err: io::Error) -> ReadError {
        ReadError::IoError(err)
    }
}

impl Universe {
    /// Create a new empty Universe
    pub fn new() -> Universe {
        Universe{
            particles: Vec::new(),
            topology: Topology::new(),
            kinds: HashMap::new(),
            interactions: Interactions::new(),
            cell: UnitCell::new(),
        }
    }

    /// Read an XYZ file and create an Universe from it.
    pub fn from_file(path: &str) -> Result<Universe, ReadError> {
        // TODO: use chemharp for implementation
        let mut file = try!(File::open(path));
        let mut content = String::new();
        try!(file.read_to_string(&mut content));
        let lines : Vec<&str> = content.lines_any().collect();
        let natoms = try!(lines[0].parse::<usize>());
        if lines.len() != natoms + 2 {
            return Err(ReadError::XYZFormatError{err: "Bad atom number in file."});
        }

        let mut universe = Universe::new();
        for line in lines[2..].iter() {
            let splited: Vec<&str> = line.split_whitespace().collect();
            if splited.len() < 4 {
                return Err(ReadError::XYZFormatError{err: "Line is too short."});
            }
            let mut part = Particle::new(splited[0]);
            let x = try!(splited[1].parse::<f64>());
            let y = try!(splited[2].parse::<f64>());
            let z = try!(splited[3].parse::<f64>());
            part.set_position(Vector3D::new(x, y, z));
            universe.add_particle(part);
        }
        return Ok(universe);
    }

    /// Create an empty universe with a specific UnitCell
    pub fn from_cell(cell: UnitCell) -> Universe {
        let mut universe = Universe::new();
        universe.set_cell(cell);
        return universe;
    }

    /// Get a reference to  the universe unit cell
    #[inline] pub fn cell<'a>(&'a self) -> &'a UnitCell {&self.cell}
    /// Set the universe unit cell
    #[inline] pub fn set_cell(&mut self, cell: UnitCell) {self.cell = cell;}

    /// Get a reference to the universe topology
    #[inline] pub fn topology<'a>(&'a self) -> &'a Topology {&self.topology}
    /// Get a mutable reference to the universe topology
    #[inline] pub fn topology_mut<'a>(&'a mut self) -> &'a mut Topology {&mut self.topology}

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
        let index = self.size() - 1;
        self.topology.add_particle(index);
    }
    /// Get the number of particles in this universe
    #[inline] pub fn size(&self) -> usize {self.particles.len()}

    /// Get the list of pair interaction between the particles at indexes `i`
    /// and `j`.
    pub fn pair_potentials<'a>(&'a self, i: usize, j: usize) -> &'a Vec<Box<PairPotential>> {
        let ikind = self.particles[i].kind();
        let jkind = self.particles[j].kind();
        match self.interactions.pairs(ikind, jkind) {
            Some(val) => &val,
            None => {
                let i = self.particles[i].name();
                let j = self.particles[j].name();
                error!("No potential defined for the pair ({}, {})", i, j);
                panic!();
            }
        }
    }

    /// Get the list of bonded interaction between the particles at indexes `i`
    /// and `j`.
    pub fn bond_potentials<'a>(&'a self, i: usize, j: usize) -> &'a Vec<Box<PairPotential>> {
        let ikind = self.particles[i].kind();
        let jkind = self.particles[j].kind();
        match self.interactions.bonds(ikind, jkind) {
            Some(val) => &val,
            None => {
                let i = self.particles[i].name();
                let j = self.particles[j].name();
                error!("No potential defined for the bond ({}, {})", i, j);
                panic!();
            }
        }
    }

    /// Get the list of angle interaction between the particles at indexes `i`,
    /// `j` and `k`.
    pub fn angle_potentials<'a>(&'a self, i: usize, j: usize, k: usize) -> &'a Vec<Box<AnglePotential>> {
        let ikind = self.particles[i].kind();
        let jkind = self.particles[j].kind();
        let kkind = self.particles[k].kind();

        match self.interactions.angles(ikind, jkind, kkind) {
            Some(val) => &val,
            None => {
                let i = self.particles[i].name();
                let j = self.particles[j].name();
                let k = self.particles[k].name();
                error!("No potential defined for the angle ({}, {}, {})", i, j, k);
                panic!();
            }
        }
    }

    /// Get the list of dihedral angles interaction between the particles at
    /// indexes `i`, `j`, `k` and `m`.
    pub fn dihedral_potentials<'a>(&'a self, i: usize, j: usize, k: usize, m: usize) -> &'a Vec<Box<DihedralPotential>> {
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
                error!("No potential defined for the dihedral ({}, {}, {}, {})", i, j, k, m);
                panic!();
            }
        }
    }

    /// Add a pair interaction between the particles with names `names`
    pub fn add_pair_interaction<T>(&mut self, i: &str, j: &str, pot: T)
    where T: PairPotential + 'static {
        let ikind = self.get_kind(i);
        let jkind = self.get_kind(j);

        self.interactions.add_pair(ikind, jkind, pot);
    }

    /// Add a bonded interaction between the particles with names `names`
    pub fn add_bond_interaction<T>(&mut self, i: &str, j: &str, pot: T)
    where T: PairPotential + 'static {
        let ikind = self.get_kind(i);
        let jkind = self.get_kind(j);

        self.interactions.add_bond(ikind, jkind, pot);
    }

    /// Add an angle interaction between the particles with names `names`
    pub fn add_angle_interaction<T>(&mut self, i: &str, j: &str, k: &str, pot: T)
    where T: AnglePotential + 'static {
        let ikind = self.get_kind(i);
        let jkind = self.get_kind(j);
        let kkind = self.get_kind(k);

        self.interactions.add_angle(ikind, jkind, kkind, pot);
    }

    /// Add an angle interaction between the particles with names `names`
    pub fn add_dihedral_interaction<T>(&mut self, i: &str, j: &str, k: &str, m: &str, pot: T)
    where T: DihedralPotential + 'static {
        let ikind = self.get_kind(i);
        let jkind = self.get_kind(j);
        let kkind = self.get_kind(k);
        let mkind = self.get_kind(m);

        self.interactions.add_dihedral(ikind, jkind, kkind, mkind, pot);
    }

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

    /// Get an iterator over the `Particle` in this universe
    #[inline] pub fn iter(&self) -> slice::Iter<Particle> {
        self.particles.iter()
    }

    /// Get a mutable iterator over the `Particle` in this universe
    #[inline] pub fn iter_mut(&mut self) -> slice::IterMut<Particle> {
        self.particles.iter_mut()
    }
}

/******************************************************************************/

use ::simulation::Compute;
use ::simulation::{PotentialEnergy, KineticEnergy, TotalEnergy};
use ::simulation::Temperature;
use ::simulation::Volume;
use ::simulation::{Virial, Stress, Pressure};

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
    fn index<'a>(&'a self, index: usize) -> &'a Particle {
        &self.particles[index]
    }
}

impl IndexMut<usize> for Universe {
    #[inline]
    fn index_mut<'a>(&'a mut self, index: usize) -> &'a mut Particle {
        &mut self.particles[index]
    }
}

#[cfg(test)]
mod tests {
    use ::universe::*;
    use ::types::*;
    use ::potentials::*;

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
    fn topology() {
        use std::collections::HashSet;

        let mut universe = Universe::new();
        universe.add_particle(Particle::new("He"));
        universe.add_particle(Particle::new("He"));
        universe.add_particle(Particle::new("He"));
        universe.add_particle(Particle::new("He"));

        {
            let topology = universe.topology_mut();
            topology.add_bond(0, 3);
            topology.add_bond(1, 2);
        }

        let topology = universe.topology();
        let mut bonds = HashSet::new();
        bonds.insert(Bond::new(1, 2));
        bonds.insert(Bond::new(0, 3));
        assert_eq!(topology.bonds(), &bonds);
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

        universe.add_pair_interaction("He", "He", LennardJones{sigma: 0.3, epsilon: 2.0});
        universe.add_pair_interaction("He", "He", Harmonic{k: 100.0, x0: 1.1});
        assert_eq!(universe.pair_potentials(0, 0).len(), 2);

        universe.add_bond_interaction("He", "He", Harmonic{k: 100.0, x0: 1.1});
        assert_eq!(universe.bond_potentials(0, 0).len(), 1);

        universe.add_angle_interaction("He", "He", "He", Harmonic{k: 100.0, x0: 1.1});
        assert_eq!(universe.angle_potentials(0, 0, 0).len(), 1);

        universe.add_dihedral_interaction("He", "He", "He", "He", CosineHarmonic::new(0.3, 2.0));
        assert_eq!(universe.dihedral_potentials(0, 0, 0, 0).len(), 1);
    }

    #[test]
    #[should_panic]
    fn pairs_errors() {
        let mut universe = Universe::new();
        universe.add_particle(Particle::new("He"));
        universe.pair_potentials(0, 0);
    }

    #[test]
    #[should_panic]
    fn angles_errors() {
        let mut universe = Universe::new();
        universe.add_particle(Particle::new("He"));
        universe.angle_potentials(0, 0, 0);
    }
    #[test]
    #[should_panic]
    fn dihedrals_errors() {
        let mut universe = Universe::new();
        universe.add_particle(Particle::new("He"));
        universe.dihedral_potentials(0, 0, 0, 0);
    }
}
