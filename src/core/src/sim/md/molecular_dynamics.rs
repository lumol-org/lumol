// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

use sim::{Propagator, TemperatureStrategy};
use sys::System;

use super::{Control, Integrator, Thermostat};
use super::VelocityVerlet;

/// Molecular Dynamics propagator for the simulation.
pub struct MolecularDynamics {
    /// The integrator we should use to propagate the equations of motion.
    integrator: Box<Integrator>,
    /// Optional thermostat algorithm
    thermostat: Option<Box<Thermostat>>,
    /// Control algorithms in the simulation.
    controls: Vec<Box<Control>>,
}

impl MolecularDynamics {
    /// Create a new `MolecularDynamics` propagator using a `VelocityVerlet`
    /// integrator.
    pub fn new(dt: f64) -> MolecularDynamics {
        MolecularDynamics::from_integrator(Box::new(VelocityVerlet::new(dt)))
    }

    /// Create a new `MolecularDynamics` propagator using the specified
    /// `integrator`.
    pub fn from_integrator(integrator: Box<Integrator>) -> MolecularDynamics {
        MolecularDynamics {
            integrator: integrator,
            thermostat: None,
            controls: Vec::new(),
        }
    }

    /// Add a control algorithm to the internal list of controls.
    pub fn add_control(&mut self, control: Box<Control>) {
        self.controls.push(control);
    }

    /// Set the thermostat to use with this simulation
    pub fn set_thermostat(&mut self, thermostat: Box<Thermostat>) {
        self.thermostat = Some(thermostat);
    }
}

impl Propagator for MolecularDynamics {
    fn temperature_strategy(&self) -> TemperatureStrategy {
        TemperatureStrategy::Velocities
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
            thermostat.control(system);
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
