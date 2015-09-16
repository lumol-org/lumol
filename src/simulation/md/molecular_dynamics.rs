/*
 * Cymbalum, Molecular Simulation in Rust
 * Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
*/

use ::universe::Universe;

use ::simulation::Propagator;
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
    fn setup(&mut self, universe: &Universe) {
        self.integrator.setup(universe);
        for control in self.controls.iter_mut() {
            control.setup(universe);
        }
    }

    fn propagate(&mut self, universe: &mut Universe) {
        self.integrator.integrate(universe);
        for control in self.controls.iter_mut() {
            control.control(universe);
        }
    }

    fn finish(&mut self, universe: &Universe) {
        for control in self.controls.iter_mut() {
            control.finish(universe);
        }
    }
}
