// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 Lumol's contributors — BSD license

use std::ops::{Deref, DerefMut};
use std::collections::BTreeMap;

use types::{Vector3D, Matrix3};

use energy::{PairInteraction, BondPotential, AnglePotential, DihedralPotential};
use energy::{GlobalPotential, CoulombicPotential};

use sys::{Configuration, Particle, ParticleKind, UnitCell};
use sys::{Composition, Interactions, EnergyEvaluator};

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

    /// Create a system with the specified `configuration`, and no interactions.
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
        if let Some(&kind) = self.kinds.get(name) {
            return kind;
        } else {
            let kind = ParticleKind(self.kinds.len() as u32);
            let _ = self.kinds.insert(String::from(name), kind);
            kind
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
        for particle in self.particles() {
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



    /// Guess the bonds in the configuration using the chemfiles algorithm.
    ///
    /// This function removes any existing bond, and tries to guess bonds using
    /// a distance criteria. Because this function does a round trip to
    /// chemfiles, it might be costly when dealing with big systems.
    pub fn guess_bonds(&mut self) {
        use ::sys::chfl::{ToChemfiles, ToLumol};
        let mut frame = self.to_chemfiles().expect(
            "can not convert the configuration to a chemfiles frame"
        );
        frame.guess_topology().expect(
            "can not guess system topology in chemfiles"
        );
        *self = frame.to_lumol().expect(
            "can not convert the chemfiles frame to a configuration"
        );
    }
}

/// Functions related to interactions
impl System {
    /// Get an helper struct to evaluate the energy of this system.
    pub fn energy_evaluator(&self) -> EnergyEvaluator {
        EnergyEvaluator::new(self)
    }

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
        let kind_i = self.particle(i).kind;
        let kind_j = self.particle(j).kind;
        let pairs = self.interactions.pairs(kind_i, kind_j);
        if pairs.is_empty() {
            warn_once!(
                "No potential defined for the pair ({}, {})",
                self.particle(i).name(), self.particle(j).name()
            );
        }
        return pairs;
    }

    /// Get read-only access to the interactions for this system
    pub(crate) fn interactions(&self) -> &Interactions {
        &self.interactions
    }

    /// Get the list of bonded potential acting between the particles at indexes
    /// `i` and `j`.
    pub fn bond_potentials(&self, i: usize, j: usize) -> &[Box<BondPotential>] {
        let kind_i = self.particle(i).kind;
        let kind_j = self.particle(j).kind;
        let bonds = self.interactions.bonds(kind_i, kind_j);
        if bonds.is_empty() {
            warn_once!(
                "No potential defined for the bond ({}, {})",
                self.particle(i).name(), self.particle(j).name()
            );
        }
        return bonds;
    }

    /// Get the list of angle interaction acting between the particles at
    /// indexes `i`, `j` and `k`.
    pub fn angle_potentials(&self, i: usize, j: usize, k: usize) -> &[Box<AnglePotential>] {
        let kind_i = self.particle(i).kind;
        let kind_j = self.particle(j).kind;
        let kind_k = self.particle(k).kind;
        let angles = self.interactions.angles(kind_i, kind_j, kind_k);
        if angles.is_empty() {
            warn_once!(
                "No potential defined for the angle ({}, {}, {})",
                self.particle(i).name(), self.particle(j).name(),
                self.particle(k).name()
            );
        }
        return angles;
    }

    /// Get the list of dihedral angles interaction acting between the particles
    /// at indexes `i`, `j`, `k` and `m`.
    pub fn dihedral_potentials(&self, i: usize, j: usize, k: usize, m: usize) -> &[Box<DihedralPotential>] {
        let kind_i = self.particle(i).kind;
        let kind_j = self.particle(j).kind;
        let kind_k = self.particle(k).kind;
        let kind_m = self.particle(m).kind;
        let dihedrals = self.interactions.dihedrals(kind_i, kind_j, kind_k, kind_m);
        if dihedrals.is_empty() {
            warn_once!(
                "No potential defined for the dihedral angle ({}, {}, {}, {})",
                self.particle(i).name(), self.particle(j).name(),
                self.particle(k).name(), self.particle(m).name()
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

#[cfg(test)]
mod tests {
    use super::System;
    use sys::{Particle, ParticleKind};

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

        assert_eq!(system.particle(0).kind, ParticleKind(0));
        assert_eq!(system.particle(1).kind, ParticleKind(1));
        assert_eq!(system.particle(2).kind, ParticleKind(0));
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
