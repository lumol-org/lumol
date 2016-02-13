// Cymbalum, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

use system::System;
use simulation::Propagator;

use super::Integrator;
use super::Control;
use super::VelocityVerlet;

/// Molecular Dynamics propagator for the simulation.
pub struct MolecularDynamics {
    /// The integrator we should use to propagate the equations of motion.
    integrator: Box<Integrator>,
    /// Control algorithms in the simulation.
    controls: Vec<Box<Control>>,
}

impl MolecularDynamics {
    /// Create a new `MolecularDynamics` propagator using a `VelocityVerlet`
    /// integrator.
    pub fn new(dt: f64) -> MolecularDynamics {
        MolecularDynamics::from_integrator(VelocityVerlet::new(dt))
    }

    /// Create a new `MolecularDynamics` propagator using the specified
    /// `integrator`.
    pub fn from_integrator<I>(integrator: I) -> MolecularDynamics where I: Integrator + 'static {
        MolecularDynamics{
            integrator: Box::new(integrator),
            controls: Vec::new(),
        }
    }

    /// Add a control algorithm to the internal list of controls.
    pub fn add_control<C>(&mut self, control: C) where C: Control + 'static {
        self.controls.push(Box::new(control));
    }
}

impl Propagator for MolecularDynamics {
    fn setup(&mut self, system: &System) {
        self.integrator.setup(system);
        for control in &mut self.controls {
            control.setup(system);
        }
    }

    fn propagate(&mut self, system: &mut System) {
        self.integrator.integrate(system);
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
