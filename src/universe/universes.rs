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

use ::potentials::PairPotential;
use ::types::Vector3D;

use super::Particle;
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
    /// Particles kinds, associating particles names and indexes
    kinds: HashMap<String, usize>,
    /// Interactions is a hash map associating particles kinds and potentials
    interactions: Interactions,
}

impl Universe {
    /// Create a new empty Universe
    pub fn new() -> Universe {
        Universe{
            particles: Vec::new(),
            kinds: HashMap::new(),
            interactions: Interactions::new(),
            cell: UnitCell::new(),
        }
    }

    /// Create an empty universe with a specific UnitCell
    pub fn from_cell(cell: UnitCell) -> Universe {
        let mut universe = Universe::new();
        universe.set_cell(cell);
        return universe;
    }

    /// Get the universe unit cell
    pub fn cell<'a>(&'a self) -> &'a UnitCell {&self.cell}
    /// Set the universe unit cell
    pub fn set_cell(&mut self, cell: UnitCell) {self.cell = cell;}

    /// Insert a particle at the end of the internal list
    pub fn add_particle(&mut self, p: Particle) {
        let mut part = p;
        if part.kind() == usize::max_value() {
            // If no value have been precised, set one from the internal list
            // of particles kinds.
            let kind = self.get_kind(part.name());
            part.set_kind(kind);
        }
        self.particles.push(part);
    }
    /// Get the number of particles in this universe
    pub fn size(&self) -> usize {self.particles.len()}

    /// Get the list of pair interaction between the atom i and the atom j
    pub fn pairs<'a>(&'a self, i: usize, j: usize) -> &'a Vec<Box<PairPotential>> {
        let ikind = self.particles[i].kind();
        let jkind = self.particles[j].kind();
        &self.interactions.pairs[&(ikind, jkind)]
    }

    /// Add an interaction between the particles with names `names`
    pub fn add_pair_interaction<T>(&mut self, i: &str, j: &str, pot: T)
    where T: PairPotential + 'static {
        let ikind = self.get_kind(i);
        let jkind = self.get_kind(j);

        if !self.interactions.pairs.contains_key(&(ikind, jkind)) {
            self.interactions.pairs.insert((ikind, jkind), Vec::new());
        }
        let pairs = self.interactions.pairs.get_mut(&(ikind, jkind)).unwrap();
        pairs.push(Box::new(pot));
    }

    /// Get the distance between the particles at indexes `i` and `j`
    pub fn distance(&self, i: usize, j:usize) -> f64 {
        self.cell.distance(self.particles[i].position(), self.particles[j].position())
    }

    /// Wrap the vector i->j in the cell.
    pub fn wrap_vector(&self, i: usize, j:usize) -> Vector3D {
        let mut res = *self.particles[i].position() - *self.particles[j].position();
        self.cell.wrap_vector(&mut res);
        return res;
    }

    /// Get or create the usize kind index for the name `name` of a particle
    fn get_kind(&mut self, name: &str) -> usize {
        if self.kinds.contains_key(name) {
            self.kinds[name]
        } else {
            let index = self.kinds.len();
            self.kinds.insert(name.to_string(), index);
            return index;
        }
    }

    /// Get an iterator over the `Particle` in this universe
    pub fn iter(&self) -> slice::Iter<Particle> {
        self.particles.iter()
    }

    /// Get a mutable iterator over the `Particle` in this universe
    pub fn iter_mut(&mut self) -> slice::IterMut<Particle> {
        self.particles.iter_mut()
    }
}

/******************************************************************************/

use ::simulation::Compute;
use ::simulation::{PotentialEnergy, KineticEnergy, TotalEnergy, Temperature};

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
    fn pairs() {
        let mut universe = Universe::new();
        universe.add_particle(Particle::new("He"));

        universe.add_pair_interaction("He", "He", LennardJones{sigma: 0.3, epsilon: 2.0});
        universe.add_pair_interaction("He", "He", Harmonic{k: 100.0, r0: 1.1});

        assert_eq!(universe.pairs(0, 0).len(), 2);
    }
}
