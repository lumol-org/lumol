// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

use crate::propagator::{Propagator, TemperatureStrategy};
use lumol_core::{System, DegreesOfFreedom};

use super::{Control, Integrator, Thermostat};
use super::VelocityVerlet;

/// Molecular Dynamics propagator for the simulation.
pub struct MolecularDynamics {
    /// The integrator we should use to propagate the equations of motion.
    integrator: Box<dyn Integrator>,
    /// Optional thermostat algorithm
    thermostat: Option<Box<dyn Thermostat>>,
    /// Control algorithms in the simulation.
    controls: Vec<Box<dyn Control>>,
}

impl MolecularDynamics {
    /// Create a new `MolecularDynamics` propagator using a `VelocityVerlet`
    /// integrator.
    pub fn new(dt: f64) -> MolecularDynamics {
        MolecularDynamics::from_integrator(Box::new(VelocityVerlet::new(dt)))
    }

    /// Create a new `MolecularDynamics` propagator using the specified
    /// `integrator`.
    pub fn from_integrator(integrator: Box<dyn Integrator>) -> MolecularDynamics {
        MolecularDynamics {
            integrator: integrator,
            thermostat: None,
            controls: Vec::new(),
        }
    }

    /// Add a control algorithm to the internal list of controls.
    pub fn add_control(&mut self, control: Box<dyn Control>) {
        self.controls.push(control);
    }

    /// Set the thermostat to use with this simulation
    pub fn set_thermostat(&mut self, thermostat: Box<dyn Thermostat>) {
        self.thermostat = Some(thermostat);
    }
}

impl Propagator for MolecularDynamics {
    fn temperature_strategy(&self) -> TemperatureStrategy {
        TemperatureStrategy::Velocities
    }

    fn degrees_of_freedom(&self, _: &System) -> DegreesOfFreedom {
        // default to particles for now. change this if/when constrains are
        // implemented
        DegreesOfFreedom::Particles
    }

    fn setup(&mut self, system: &System) {
        self.integrator.setup(system);
        for control in &mut self.controls {
            control.setup(system);
        }
    }

    fn propagate(&mut self, system: &mut System) {
        self.integrator.integrate(system);

        if let Some(ref mut thermostat) = self.thermostat {
            thermostat.apply(system);
        }

        for control in &mut self.controls {
            control.control(system);
        }
    }

    fn finish(&mut self, system: &System) {
        for control in &mut self.controls {
            control.finish(system);
        }
    }
}
