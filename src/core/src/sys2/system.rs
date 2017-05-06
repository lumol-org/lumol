// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 Lumol's contributors — BSD license

use std::ops::{Deref, DerefMut};
use std::slice;
use std::collections::BTreeMap;

use energy::{PairInteraction, BondPotential, AnglePotential, DihedralPotential};
use energy::{GlobalPotential, CoulombicPotential};

use sys::{Configuration, Particle, ParticleKind, UnitCell};
use sys::Composition;

use sys2::Interactions;

/// The `System` type hold all the data about a simulated system.
///
/// This data contains:
///
///   - an unit cell, containing the system;
///   - a list of particles in the system;
///   - a list of molecules in the system;
///   - a list of interactions, associating particles kinds and potentials
///
/// In the implementation, the particles contained in a molecule are guaranteed
/// to be contiguous in memory. This allow for faster access when iterating
/// over molecules, and easier molecule removal from the system.
#[derive(Clone)]
pub struct System {
    /// The system configuration
    configuration: Configuration,
    /// All the interactions in this system
    interactions: Interactions,
    /// Association particles names to particle kinds
    kinds: BTreeMap<String, ParticleKind>,
    /// The current simulation step
    step: u64,
    /// Externally managed temperature for the system
    external_temperature: Option<f64>,
}

impl System {
    /// Create a new empty `System`
    pub fn new() -> System {
        System::with_configuration(Configuration::new())
    }

    /// Create an empty system with a specific unit cell
    pub fn with_cell(cell: UnitCell) -> System {
        let mut configuration = Configuration::new();
        configuration.cell = cell;
        System::with_configuration(configuration)
    }

    /// Create a system with the specified `configuration`
    fn with_configuration(configuration: Configuration) -> System {
        System {
            configuration: configuration,
            kinds: BTreeMap::new(),
            interactions: Interactions::new(),
            step: 0,
            external_temperature: None,
        }
    }

    fn get_kind(&mut self, name: &str) -> ParticleKind {
        match self.kinds.get(name).cloned() {
            Some(kind) => kind,
            None => {
                let kind = ParticleKind(self.kinds.len() as u32);
                let _ = self.kinds.insert(String::from(name), kind);
                kind
            }
        }
    }

    /// Insert a particle at the end of the internal list.
    pub fn add_particle(&mut self, mut particle: Particle) {
        if particle.kind == ParticleKind::invalid() {
            particle.kind = self.get_kind(particle.name());
        }
        self.configuration.add_particle(particle);
    }

    /// Get the number of particles of each kind in the configuration
    pub fn composition(&self) -> Composition {
        let mut composition = Composition::new();
        composition.resize(self.kinds.len());
        for particle in self {
            composition[particle.kind] += 1;
        }
        return composition;
    }

    /// Get a list of all the particles kinds in the system.
    pub fn particle_kinds(&self) -> Vec<ParticleKind> {
        self.kinds.values().cloned().collect()
    }

    /// Get the current step of the system
    pub fn step(&self) -> u64 {
        self.step
    }

    /// Increment the system step
    pub fn increment_step(&mut self) {
        self.step += 1;
    }

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
}

/// Functions related to interactions
impl System {
    /// Add the `potential` pair interaction for the pair `(i, j)`
    pub fn add_pair_potential(&mut self, i: &str, j: &str, potential: PairInteraction) {
        let kind_i = self.get_kind(i);
        let kind_j = self.get_kind(j);
        self.interactions.add_pair(kind_i, kind_j, potential)
    }

    /// Add the `potential` bonded interaction for the pair `(i, j)`
    pub fn add_bond_potential(&mut self, i: &str, j: &str, potential: Box<BondPotential>) {
        let kind_i = self.get_kind(i);
        let kind_j = self.get_kind(j);
        self.interactions.add_bond(kind_i, kind_j, potential)
    }

    /// Add the `potential` angle interaction for the angle `(i, j, k)`
    pub fn add_angle_potential(&mut self, i: &str, j: &str, k: &str, potential: Box<AnglePotential>) {
        let kind_i = self.get_kind(i);
        let kind_j = self.get_kind(j);
        let kind_k = self.get_kind(k);
        self.interactions.add_angle(kind_i, kind_j, kind_k, potential)
    }

    /// Add the `potential` dihedral interaction for the dihedral angle `(i, j,
    /// k, m)`
    pub fn add_dihedral_potential(&mut self, i: &str, j: &str, k: &str, m: &str, potential: Box<DihedralPotential>) {
        let kind_i = self.get_kind(i);
        let kind_j = self.get_kind(j);
        let kind_k = self.get_kind(k);
        let kind_m = self.get_kind(m);
        self.interactions.add_dihedral(kind_i, kind_j, kind_k, kind_m, potential)
    }

    /// Set the coulombic interaction for all pairs to `potential`
    pub fn set_coulomb_potential(&mut self, potential: Box<CoulombicPotential>) {
        self.interactions.coulomb = Some(potential);
    }

    /// Add the `potential` global interaction
    pub fn add_global_potential(&mut self, potential: Box<GlobalPotential>) {
        self.interactions.globals.push(potential);
    }

    /// Get the list of pair potential acting between the particles at indexes
    /// `i` and `j`.
    pub fn pair_potentials(&self, i: usize, j: usize) -> &[PairInteraction] {
        let kind_i = self[i].kind;
        let kind_j = self[j].kind;
        let pairs = self.interactions.pairs(kind_i, kind_j);
        if pairs.is_empty() {
            warn_once!(
                "No potential defined for the pair ({}, {})",
                self[i].name(), self[j].name()
            );
        }
        return pairs;
    }

    /// Get the list of bonded potential acting between the particles at indexes
    /// `i` and `j`.
    pub fn bond_potentials(&self, i: usize, j: usize) -> &[Box<BondPotential>] {
        let kind_i = self[i].kind;
        let kind_j = self[j].kind;
        let bonds = self.interactions.bonds(kind_i, kind_j);
        if bonds.is_empty() {
            warn_once!(
                "No potential defined for the bond ({}, {})",
                self[i].name(), self[j].name()
            );
        }
        return bonds;
    }

    /// Get the list of angle interaction acting between the particles at
    /// indexes `i`, `j` and `k`.
    pub fn angle_potentials(&self, i: usize, j: usize, k: usize) -> &[Box<AnglePotential>] {
        let kind_i = self[i].kind;
        let kind_j = self[j].kind;
        let kind_k = self[k].kind;
        let angles = self.interactions.angles(kind_i, kind_j, kind_k);
        if angles.is_empty() {
            warn_once!(
                "No potential defined for the angle ({}, {}, {})",
                self[i].name(), self[j].name(), self[k].name()
            );
        }
        return angles;
    }

    /// Get the list of dihedral angles interaction acting between the particles
    /// at indexes `i`, `j`, `k` and `m`.
    pub fn dihedral_potentials(&self, i: usize, j: usize, k: usize, m: usize) -> &[Box<DihedralPotential>] {
        let kind_i = self[i].kind;
        let kind_j = self[j].kind;
        let kind_k = self[k].kind;
        let kind_m = self[m].kind;
        let dihedrals = self.interactions.dihedrals(kind_i, kind_j, kind_k, kind_m);
        if dihedrals.is_empty() {
            warn_once!(
                "No potential defined for the dihedral angle ({}, {}, {}, {})",
                self[i].name(), self[j].name(), self[k].name(), self[m].name()
            );
        }
        return dihedrals;
    }

    /// Get the coulombic interaction for the system
    pub fn coulomb_potential(&self) -> Option<&Box<CoulombicPotential>> {
        self.interactions.coulomb.as_ref()
    }

    /// Get all global interactions for the system
    pub fn global_potentials(&self) -> &[Box<GlobalPotential>] {
        &self.interactions.globals
    }

    /// Get maximum cutoff from `coulomb`, `pairs` and `global` interactions.
    pub fn maximum_cutoff(&self) -> Option<f64> {
        self.interactions.maximum_cutoff()
    }
}

impl Deref for System {
    type Target = Configuration;

    fn deref(&self) -> &Configuration {
        &self.configuration
    }
}

impl DerefMut for System {
    fn deref_mut(&mut self) -> &mut Configuration {
        &mut self.configuration
    }
}

impl<'a> IntoIterator for &'a System {
    type Item = &'a Particle;
    type IntoIter = slice::Iter<'a, Particle>;

    #[inline]
    fn into_iter(self) -> slice::Iter<'a, Particle> {
        self.configuration.iter()
    }
}

impl<'a> IntoIterator for &'a mut System {
    type Item = &'a mut Particle;
    type IntoIter = slice::IterMut<'a, Particle>;

    #[inline]
    fn into_iter(self) -> slice::IterMut<'a, Particle> {
        self.configuration.iter_mut()
    }
}

#[cfg(test)]
mod tests {
    use super::System;
    use sys::*;

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
    #[should_panic]
    fn negative_external_temperature() {
        let mut system = System::new();
        system.external_temperature(Some(-1.0));
    }

    #[test]
    fn deref() {
        let mut system = System::new();
        system.add_particle(Particle::new("H"));
        system.add_particle(Particle::new("O"));
        system.add_particle(Particle::new("H"));
        assert_eq!(system.molecules().len(), 3);

        // This uses deref_mut
        let _ = system.add_bond(0, 1);
        let _ = system.add_bond(2, 1);

        // This uses deref
        assert_eq!(system.molecules().len(), 1);
    }

    #[test]
    fn add_particle() {
        let mut system = System::new();
        system.add_particle(Particle::new("H"));
        system.add_particle(Particle::new("O"));
        system.add_particle(Particle::new("H"));

        assert_eq!(system[0].kind, ParticleKind(0));
        assert_eq!(system[1].kind, ParticleKind(1));
        assert_eq!(system[2].kind, ParticleKind(0));
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

        let composition = system.composition();
        assert_eq!(composition.len(), 4);
        assert_eq!(composition[ParticleKind(0)], 3);
        assert_eq!(composition[ParticleKind(1)], 2);
        assert_eq!(composition[ParticleKind(2)], 1);
        assert_eq!(composition[ParticleKind(3)], 1);
    }

    #[test]
    fn missing_interaction() {
        let mut system = System::new();
        system.add_particle(Particle::new("He"));
        system.add_particle(Particle::new("He"));
        system.add_particle(Particle::new("He"));
        system.add_particle(Particle::new("He"));
        assert_eq!(system.pair_potentials(0, 0).len(), 0);
        assert_eq!(system.bond_potentials(0, 0).len(), 0);
        assert_eq!(system.angle_potentials(0, 0, 0).len(), 0);
        assert_eq!(system.dihedral_potentials(0, 0, 0, 0).len(), 0);
    }
}
