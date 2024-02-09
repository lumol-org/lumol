// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 Lumol's contributors — BSD license

use std::ops::{Deref, DerefMut};
use std::cmp::{max, min};

use soa_derive::soa_zip;
use log_once::warn_once;

use crate::{Matrix3, Vector3D};
use crate::{AnglePotential, BondPotential, DihedralPotential, PairInteraction};
use crate::{CoulombicPotential, GlobalPotential};
use crate::{Composition, EnergyEvaluator, Interactions};
use crate::{Configuration, Molecule, UnitCell};

/// The number of degrees of freedom simulated in a given system
#[derive(Clone, PartialEq, Debug)]
pub enum DegreesOfFreedom {
    /// All particles are explicitly simulated
    Particles,
    /// All molecules are simulated as rigid bodies
    Molecules,
    /// All particles are explicitly simulated, but some degrees of freedom
    /// are frozen. The usize value is the number of frozen degree of freedom.
    Frozen(usize),
}

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
    /// Externally managed temperature for the system
    external_temperature: Option<f64>,
    /// Number of degrees of freedom simulated in the system. This default to
    /// `DegreesOfFreedom::Particles`, and is set in the simulation setup.
    pub simulated_degrees_of_freedom: DegreesOfFreedom,
    /// The current simulation step
    pub step: u64,
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
            interactions: Interactions::new(),
            step: 0,
            external_temperature: None,
            simulated_degrees_of_freedom: DegreesOfFreedom::Particles,
        }
    }

    /// Add a molecule to the system
    pub fn add_molecule(&mut self, mut molecule: Molecule) {
        for (kind, name) in soa_zip!(molecule.particles_mut(), [mut kind, name]) {
            *kind = self.interactions.get_kind(name);
        }
        self.configuration.add_molecule(molecule);
    }

    /// Get the composition in particles and molecules of the configuration
    pub fn composition(&self) -> Composition {
        let mut composition = Composition::new();
        for &kind in self.particles().kind {
            composition.add_particle(kind);
        }
        for molecule in self.molecules() {
            composition.add_molecule(molecule.hash());
        }
        return composition;
    }

    /// Use an external temperature for all the system properties. Calling this
    /// with `Some(temperature)` will replace all the computation of the
    /// temperature from the velocities with the given values. Calling it with
    /// `None` will use the velocities.
    ///
    /// The default is to use the velocities unless this function is called.
    pub fn simulated_temperature(&mut self, temperature: Option<f64>) {
        if let Some(temperature) = temperature {
            assert!(temperature >= 0.0, "External temperature must be positive");
        }
        self.external_temperature = temperature;
    }
}

/// Functions related to interactions
impl System {
    /// Get an helper struct to evaluate the energy of this system.
    pub fn energy_evaluator(&self) -> EnergyEvaluator<'_> {
        EnergyEvaluator::new(self)
    }

    /// Set the pair interaction `potential` for atoms with types `i` and `j`
    #[allow(clippy::manual_assert)]
    pub fn set_pair_potential(&mut self, (i, j): (&str, &str), potential: PairInteraction) {
        if self.cell.lengths().iter().any(|&d| 0.5 * d < potential.cutoff()) {
            panic!(
                "Can not add a potential with a cutoff bigger than half of the \
                smallest cell length. Try increasing the cell size or decreasing \
                the cutoff."
            );
        }
        self.interactions.set_pair((i, j), potential);
    }

    /// Set the bond interaction `potential` for atoms with types `i` and `j`
    pub fn set_bond_potential(&mut self, (i, j): (&str, &str), potential: Box<dyn BondPotential>) {
        self.interactions.set_bond((i, j), potential);
    }

    /// Set the angle interaction `potential` for atoms with types `i`, `j`, and `k`
    pub fn set_angle_potential(
        &mut self,
        (i, j, k): (&str, &str, &str),
        potential: Box<dyn AnglePotential>,
    ) {
        self.interactions.set_angle((i, j, k), potential);
    }

    /// Set the dihedral angle interaction `potential` for atoms with types
    /// `i`, `j`, `k`, and `m`.
    pub fn set_dihedral_potential(
        &mut self,
        (i, j, k, m): (&str, &str, &str, &str),
        potential: Box<dyn DihedralPotential>,
    ) {
        self.interactions.set_dihedral((i, j, k, m), potential);
    }

    /// Set the coulombic interaction for all pairs to `potential`
    #[allow(clippy::manual_assert)]
    pub fn set_coulomb_potential(&mut self, potential: Box<dyn CoulombicPotential>) {
        if let Some(cutoff) = potential.cutoff() {
            if self.cell.lengths().iter().any(|&d| 0.5 * d < cutoff) {
                panic!(
                    "Can not add a potential with a cutoff bigger than half of the \
                    smallest cell length. Try increasing the cell size or decreasing \
                    the cutoff."
                );
            }
        }
        self.interactions.coulomb = Some(potential);
    }

    /// Add the `potential` global interaction
    pub fn add_global_potential(&mut self, potential: Box<dyn GlobalPotential>) {
        self.interactions.globals.push(potential);
    }

    /// Get the pair potential acting between the particles at indexes `i` and `j`.
    pub fn pair_potential(&self, i: usize, j: usize) -> Option<&PairInteraction> {
        let kind_i = self.particles().kind[i];
        let kind_j = self.particles().kind[j];
        return self.interactions.pair((kind_i, kind_j));
    }

    /// Get read-only access to the interactions for this system
    pub(crate) fn interactions(&self) -> &Interactions {
        &self.interactions
    }

    /// Get the bond potential acting between the particles at indexes `i` and
    /// `j`.
    pub fn bond_potential(&self, i: usize, j: usize) -> Option<&dyn BondPotential> {
        let kind_i = self.particles().kind[i];
        let kind_j = self.particles().kind[j];
        return self.interactions.bond((kind_i, kind_j));
    }

    /// Get the angle potential acting between the particles at indexes `i`, `j`
    /// and `k`.
    pub fn angle_potential(&self, i: usize, j: usize, k: usize) -> Option<&dyn AnglePotential> {
        let kind_i = self.particles().kind[i];
        let kind_j = self.particles().kind[j];
        let kind_k = self.particles().kind[k];
        return self.interactions.angle((kind_i, kind_j, kind_k));
    }

    /// Get the dihedral angles potential acting between the particles at
    /// indexes `i`, `j`, `k` and `m`.
    pub fn dihedral_potential(
        &self,
        i: usize,
        j: usize,
        k: usize,
        m: usize,
    ) -> Option<&dyn DihedralPotential> {
        let kind_i = self.particles().kind[i];
        let kind_j = self.particles().kind[j];
        let kind_k = self.particles().kind[k];
        let kind_m = self.particles().kind[m];
        return self.interactions.dihedral((kind_i, kind_j, kind_k, kind_m));
    }

    /// Get the coulombic interaction for the system
    pub fn coulomb_potential(&self) -> Option<&dyn CoulombicPotential> {
        self.interactions.coulomb.as_deref()
    }

    /// Get all global interactions for the system
    pub fn global_potentials(&self) -> &[Box<dyn GlobalPotential>] {
        &self.interactions.globals
    }

    /// Get maximum cutoff from `coulomb`, `pairs` and `global` interactions.
    pub fn maximum_cutoff(&self) -> Option<f64> {
        self.interactions.maximum_cutoff()
    }
}

use crate::compute::{KineticEnergy, PotentialEnergy, TotalEnergy};
use crate::compute::{Pressure, Stress, Virial};
use crate::compute::{PressureAtTemperature, StressAtTemperature};
use crate::compute::Compute;
use crate::compute::Forces;
use crate::compute::Temperature;
use crate::compute::Volume;

/// Functions to get physical properties of a system.
impl System {
    /// Get the number of degrees of freedom in the system
    pub fn degrees_of_freedom(&self) -> usize {
        match self.simulated_degrees_of_freedom {
            DegreesOfFreedom::Particles => 3 * self.size(),
            DegreesOfFreedom::Frozen(frozen) => 3 * self.size() - frozen,
            DegreesOfFreedom::Molecules => 3 * self.molecules().count(),
        }
    }

    /// Get the kinetic energy of the system.
    pub fn kinetic_energy(&self) -> f64 {
        KineticEnergy.compute(self)
    }

    /// Get the potential energy of the system.
    pub fn potential_energy(&self) -> f64 {
        PotentialEnergy.compute(self)
    }

    /// Get the total energy of the system.
    pub fn total_energy(&self) -> f64 {
        TotalEnergy.compute(self)
    }

    /// Get the temperature of the system.
    pub fn temperature(&self) -> f64 {
        match self.external_temperature {
            Some(value) => value,
            None => Temperature.compute(self),
        }
    }

    /// Get the volume of the system.
    pub fn volume(&self) -> f64 {
        Volume.compute(self)
    }

    /// Get the virial of the system as a tensor
    pub fn virial(&self) -> Matrix3 {
        Virial.compute(self)
    }

    /// Get the pressure of the system from the virial equation, at the system
    /// instantaneous temperature.
    pub fn pressure(&self) -> f64 {
        match self.external_temperature {
            Some(temperature) => {
                PressureAtTemperature {
                    temperature: temperature,
                }.compute(self)
            }
            None => Pressure.compute(self),
        }
    }

    /// Get the stress tensor of the system from the virial equation.
    pub fn stress(&self) -> Matrix3 {
        match self.external_temperature {
            Some(temperature) => {
                StressAtTemperature {
                    temperature: temperature,
                }.compute(self)
            }
            None => Stress.compute(self),
        }
    }

    /// Get the forces acting on all the particles in the system
    pub fn forces(&self) -> Vec<Vector3D> {
        Forces.compute(self)
    }
}

impl System {
    /// Check the system before running a simulation
    pub fn check(&self) {
        self.check_potentials();
    }

    fn check_potentials(&self) {
        // Check pairs
        for i in 0..self.size() {
            let kind_i = self.particles().kind[i];
            for j in (i + 1)..self.size() {
                let kind_j = self.particles().kind[j];
                if self.interactions.pair((kind_i, kind_j)).is_none() {
                    warn_once!(
                        "no potential defined for the pair {:?}",
                        self.sorted_names_pair(i, j)
                    );
                }
            }
        }

        // check molecular potentials
        for molecule in self.molecules() {
            for bond in molecule.bonds() {
                let kind_i = self.particles().kind[bond.i()];
                let kind_j = self.particles().kind[bond.j()];
                if self.interactions.bond((kind_i, kind_j)).is_none() {
                    warn_once!(
                        "no potential defined for the bond {:?}",
                        self.sorted_names_pair(bond.i(), bond.j())
                    );
                }
            }

            for angle in molecule.angles() {
                let kind_i = self.particles().kind[angle.i()];
                let kind_j = self.particles().kind[angle.j()];
                let kind_k = self.particles().kind[angle.k()];
                if self.interactions.angle((kind_i, kind_j, kind_k)).is_none() {
                    warn_once!(
                        "no potential defined for the angle {:?}",
                        self.sorted_names_angle(angle.i(), angle.j(), angle.k())
                    );
                }
            }

            for dihedral in molecule.dihedrals() {
                let kind_i = self.particles().kind[dihedral.i()];
                let kind_j = self.particles().kind[dihedral.j()];
                let kind_k = self.particles().kind[dihedral.k()];
                let kind_m = self.particles().kind[dihedral.m()];
                if self.interactions.dihedral((kind_i, kind_j, kind_k, kind_m)).is_none() {
                    warn_once!(
                        "no potential defined for the dihedral angle {:?}",
                        self.sorted_names_dihedral(
                            dihedral.i(), dihedral.j(), dihedral.k(), dihedral.m()
                        )
                    );
                }
            }
        }

        // check the need for a coulombic potential
        let charge2 = self.particles().charge.iter().map(|q| q * q).sum::<f64>();
        if charge2 > 1e-3 && self.interactions.coulomb.is_none() {
            warn_once!("no coulombic potential solver defined, but the system is charged");
        }
    }

    fn sorted_names_pair(&self, i: usize, j: usize) -> (&str, &str) {
        // Use the same sorting as interactions
        let name_i = &self.particles().name[i];
        let name_j = &self.particles().name[j];
        if name_i < name_j {
            (name_i, name_j)
        } else {
            (name_j, name_i)
        }
    }

    fn sorted_names_angle(&self, i: usize, j: usize, k: usize) -> (&str, &str, &str) {
        // Use the same sorting as interactions
        let name_i = &self.particles().name[i];
        let name_j = &self.particles().name[j];
        let name_k = &self.particles().name[k];
        if name_i < name_k {
            (name_i, name_j, name_k)
        } else {
            (name_k, name_j, name_i)
        }
    }

    fn sorted_names_dihedral(&self, i: usize, j: usize, k: usize, m: usize) -> (&str, &str, &str, &str) {
        // Use the same sorting as interactions
        let name_i = &self.particles().name[i];
        let name_j = &self.particles().name[j];
        let name_k = &self.particles().name[k];
        let name_m = &self.particles().name[m];

        match (max(name_i, name_j), max(name_k, name_m)) {
            (ij, km) if ij == km => {
                if min(name_i, name_j) < min(name_k, name_m) {
                    (name_i, name_j, name_k, name_m)
                } else {
                    (name_m, name_k, name_j, name_i)
                }
            },
            (ij, km) if ij < km => (name_i, name_j, name_k, name_m),
            (_, _) => (name_m, name_k, name_j, name_i),
        }
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
    use crate::{System, Molecule, Particle, ParticleKind};

    #[test]
    #[should_panic(expected="External temperature must be positive")]
    fn negative_simulated_temperature() {
        let mut system = System::new();
        system.simulated_temperature(Some(-1.0));
    }

    #[test]
    fn deref() {
        let mut system = System::new();
        system.add_molecule(Molecule::new(Particle::new("H")));
        system.add_molecule(Molecule::new(Particle::new("O")));
        system.add_molecule(Molecule::new(Particle::new("H")));
        assert_eq!(system.molecules().count(), 3);

        // This uses deref_mut
        let _ = system.add_bond(0, 1);
        let _ = system.add_bond(2, 1);

        // This uses deref
        assert_eq!(system.molecules().count(), 1);
    }

    #[test]
    fn add_molecule() {
        let mut system = System::new();
        system.add_molecule(Molecule::new(Particle::new("H")));
        system.add_molecule(Molecule::new(Particle::new("O")));
        system.add_molecule(Molecule::new(Particle::new("H")));

        assert_eq!(system.particles().kind[0], ParticleKind(0));
        assert_eq!(system.particles().kind[1], ParticleKind(1));
        assert_eq!(system.particles().kind[2], ParticleKind(0));
    }

    #[test]
    fn composition() {
        let mut system = System::new();
        system.add_molecule(Molecule::new(Particle::new("H")));
        system.add_molecule(Molecule::new(Particle::new("O")));
        system.add_molecule(Molecule::new(Particle::new("O")));
        system.add_molecule(Molecule::new(Particle::new("H")));
        system.add_molecule(Molecule::new(Particle::new("C")));
        system.add_molecule(Molecule::new(Particle::new("U")));
        system.add_molecule(Molecule::new(Particle::new("H")));

        let composition = system.composition();
        assert_eq!(composition.particles(ParticleKind(0)), 3);
        assert_eq!(composition.particles(ParticleKind(1)), 2);
        assert_eq!(composition.particles(ParticleKind(2)), 1);
        assert_eq!(composition.particles(ParticleKind(3)), 1);
    }

    #[test]
    fn missing_interaction() {
        let mut system = System::new();
        system.add_molecule(Molecule::new(Particle::new("He")));
        system.add_molecule(Molecule::new(Particle::new("He")));
        system.add_molecule(Molecule::new(Particle::new("He")));
        system.add_molecule(Molecule::new(Particle::new("He")));
        assert!(system.pair_potential(0, 0).is_none());
        assert!(system.bond_potential(0, 0).is_none());
        assert!(system.angle_potential(0, 0, 0).is_none());
        assert!(system.dihedral_potential(0, 0, 0, 0).is_none());
    }

    #[test]
    fn check_potentials() {
        use std::sync::{Arc, Mutex};
        use std::fmt::Write;
        struct TestLogger {
            message: Arc<Mutex<String>>,
        }

        impl log::Log for TestLogger {
            fn enabled(&self, _: &log::Metadata<'_>) -> bool {true}

            fn log(&self, record: &log::Record<'_>) {
                if record.level() == log::Level::Warn {
                    let mut message = self.message.lock().unwrap();
                    writeln!(&mut *message, "{}", record.args()).unwrap();
                }
            }

            fn flush(&self) {}
        }

        let message = Arc::new(Mutex::new(String::from("\n")));
        let logger = TestLogger {message: message.clone()};
        log::set_boxed_logger(Box::new(logger)).unwrap();
        log::set_max_level(log::LevelFilter::Info);


        let mut system = System::new();
        system.add_molecule(Molecule::new(Particle::new("He")));
        system.add_molecule(Molecule::new(Particle::new("He")));
        system.add_molecule(Molecule::new(Particle::new("Ar")));
        system.add_molecule(Molecule::new(Particle::new("Ar")));
        let _ = system.add_bond(0, 1);
        let _ = system.add_bond(1, 2);
        let _ = system.add_bond(2, 3);

        system.check();

        let expected_warnings = r#"
no potential defined for the pair ("He", "He")
no potential defined for the pair ("Ar", "He")
no potential defined for the pair ("Ar", "Ar")
no potential defined for the bond ("He", "He")
no potential defined for the bond ("Ar", "Ar")
no potential defined for the bond ("Ar", "He")
no potential defined for the angle ("Ar", "Ar", "He")
no potential defined for the angle ("Ar", "He", "He")
no potential defined for the dihedral angle ("Ar", "Ar", "He", "He")
"#;

        let messages = message.lock().unwrap();
        assert_eq!(expected_warnings.lines().count(), messages.lines().count());
        for line in messages.lines() {
            assert!(expected_warnings.contains(line));
        }
    }
}
