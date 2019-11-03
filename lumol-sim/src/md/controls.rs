// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

//! Control algorithms operate on the system state after the integration, and
//! can be used to adjust the simulated system in various ways.

use soa_derive::soa_zip;

use lumol_core::System;
use lumol_core::{Matrix3, Vector3D};

/// Trait for controlling some parameters in a system during a simulation.
pub trait Control {
    /// Function called once at the beginning of the simulation, which allow
    /// for some setup of the control algorithm if needed.
    fn setup(&mut self, _: &System) {}

    /// This will be called once at every step of the simulation, after the
    /// integration step.
    fn control(&mut self, system: &mut System);

    /// Function called once at the end of the simulation.
    fn finish(&mut self, _: &System) {}
}

/// Remove global translation from the system
pub struct RemoveTranslation;

impl Control for RemoveTranslation {
    fn control(&mut self, system: &mut System) {
        let total_mass = system.particles().mass.iter().sum();

        let mut com_velocity = Vector3D::zero();
        for (&mass, velocity) in soa_zip!(system.particles(), [mass, velocity]) {
            com_velocity += velocity * mass / total_mass;
        }

        for velocity in system.particles_mut().velocity {
            *velocity -= com_velocity;
        }
    }
}

/// Remove global rotation from the system
pub struct RemoveRotation;

impl Control for RemoveRotation {
    fn control(&mut self, system: &mut System) {
        // Center-of-mass
        let com = system.center_of_mass();

        // Angular momentum
        let mut moment = Vector3D::zero();
        let mut inertia = Matrix3::zero();
        for (&mass, position, velocity) in soa_zip!(system.particles(), [mass, position, velocity]) {
            let delta = position - com;
            moment += mass * (delta ^ velocity);
            inertia += -mass * delta.tensorial(&delta);
        }

        let trace = inertia.trace();
        inertia[0][0] += trace;
        inertia[1][1] += trace;
        inertia[2][2] += trace;

        // The angular velocity omega is defined by `L = I w` with L the angular
        // momentum, and I the inertia matrix.
        let angular = inertia.inverse() * moment;
        for (position, velocity) in soa_zip!(system.particles_mut(), [position, mut velocity]) {
            *velocity -= (position - com) ^ angular;
        }
    }
}

/// Rewrap all molecules' centers of mass to lie within the unit cell.
///
/// Individual atoms in a molecule may still lie outside of the cell.
pub struct Rewrap;

impl Control for Rewrap {
    fn control(&mut self, system: &mut System) {
        let cell = system.cell;
        for mut molecule in system.molecules_mut() {
            molecule.wrap(&cell);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use lumol_core::{Particle, Molecule, System, UnitCell};

    #[test]
    fn remove_translation() {
        let mut system = System::with_cell(UnitCell::cubic(10.0));
        system.add_molecule(Molecule::new(Particle::with_position("Ag", [0.0, 0.0, 0.0].into())));
        system.add_molecule(Molecule::new(Particle::with_position("Ag", [1.0, 1.0, 1.0].into())));

        system.particles_mut().velocity[0] = [1.0, 2.0, 0.0].into();
        system.particles_mut().velocity[1] = [1.0, 0.0, 0.0].into();

        RemoveTranslation.control(&mut system);
        assert_eq!(system.particles().velocity[0], Vector3D::new(0.0, 1.0, 0.0));
        assert_eq!(system.particles().velocity[1], Vector3D::new(0.0, -1.0, 0.0));
    }

    #[test]
    fn remove_rotation() {
        let mut system = System::with_cell(UnitCell::cubic(10.0));
        system.add_molecule(Molecule::new(Particle::with_position("Ag", [0.0, 0.0, 0.0].into())));
        system.add_molecule(Molecule::new(Particle::with_position("Ag", [1.0, 0.0, 0.0].into())));

        system.particles_mut().velocity[0] = [0.0, 1.0, 0.0].into();
        system.particles_mut().velocity[1] = [0.0, -1.0, 2.0].into();

        RemoveRotation.control(&mut system);
        assert_eq!(system.particles().velocity[0], Vector3D::new(0.0, 0.0, 1.0));
        assert_eq!(system.particles().velocity[1], Vector3D::new(0.0, 0.0, 1.0));
    }

    #[test]
    fn rewrap() {
        let mut system = System::with_cell(UnitCell::cubic(10.0));
        system.add_molecule(Molecule::new(Particle::with_position("Ag", [0.0, 0.0, 0.0].into())));
        system.add_molecule(Molecule::new(Particle::with_position("Ag", [15.0, 0.0, 0.0].into())));

        Rewrap.control(&mut system);
        assert_eq!(system.particles().position[0], Vector3D::new(0.0, 0.0, 0.0));
        assert_eq!(system.particles().position[1], Vector3D::new(5.0, 0.0, 0.0));
    }
}
