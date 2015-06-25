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
use super::VelocityVerlet;

/// Molecular Dynamics propagator for the simulation.
pub struct MolecularDynamics {
    /// The integrator we should use to propagate the equations of motion.
    integrator: Box<Integrator>,
}

impl MolecularDynamics {
    pub fn new(dt: f64) -> MolecularDynamics {
        MolecularDynamics{
            integrator: Box::new(VelocityVerlet::new(dt)),
        }
    }

    pub fn from_integrator<I>(integrator: I) -> MolecularDynamics where I: Integrator + 'static {
        MolecularDynamics{
            integrator: Box::new(integrator),
        }
    }
}

impl Propagator for MolecularDynamics {
    fn setup(&mut self, universe: &Universe) {
        self.integrator.setup(universe);
    }

    fn propagate(&mut self, universe: &mut Universe) {
        self.integrator.integrate(universe);
    }
}
