/*
 * Cymbalum, Molecular Simulation in Rust
 * Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
*/

use ::universe::Universe;
use ::types::Vector3D;

/// The Integrator trait define integrators interface for Molecular Dynamics.
/// An integrator is an algorithm responsible for propagating the equations of
/// motion in the system.
pub trait Integrator {
    /// Setup the integrator. This function is called once by every simulation
    /// run.
    fn setup(&mut self, _: &Universe) {}
    /// Integrate the equations of motion. This is called at every step of the
    /// simulation.
    fn integrate(&mut self, universe: &mut Universe);
}

/// Velocity-Verlet integrator. This one is reversible and symplectic.
pub struct VelocityVerlet {
    /// Timestep for the integrator
    timestep: f64,
    /// Caching the accelerations
    accelerations: Vec<Vector3D>
}

impl VelocityVerlet {
    /// Create a new integrator with a timestep of `timestep`.
    pub fn new(timestep: f64) -> VelocityVerlet {
        VelocityVerlet{
            timestep: timestep,
            accelerations: Vec::new(),
        }
    }
}

impl Integrator for VelocityVerlet {
    fn setup(&mut self, universe: &Universe) {
        // TODO(unstable): use Vec::resize here
        self.accelerations.clear();
        self.accelerations.reserve_exact(universe.size());
        for _ in 0..universe.size() {
            self.accelerations.push(Vector3D::new(0.0, 0.0, 0.0));
        }
    }

    fn integrate(&mut self, universe: &mut Universe) {
        let dt = self.timestep;
        let natoms = universe.size();

        // Update velocities at t + ∆t/2 and positions at t + ∆t
        for i in 0..natoms {
            universe[i].add_velocity(0.5 * dt * self.accelerations[i]);
            let new_vel = *universe[i].velocity() * dt;
            universe[i].add_position(new_vel);
        }

        let forces = universe.forces();
        // Update accelerations at t + ∆t and velocities at t + ∆t
        for i in 0..natoms {
            self.accelerations[i] = forces[i] / universe[i].mass();
            universe[i].add_velocity(0.5 * dt * self.accelerations[i]);
        }
    }
}
