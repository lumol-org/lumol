/*
 * Cymbalum, Molecular Simulation in Rust
 * Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
*/

use ::types::{Vector3D, Matrix3};
use ::simulation::{Compute, Forces};
use ::universe::Universe;

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
    /// Storing the accelerations
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
        self.accelerations = vec![Vector3D::new(0.0, 0.0, 0.0); universe.size()];
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

        let forces = Forces.compute(&universe);
        // Update accelerations at t + ∆t and velocities at t + ∆t
        for i in 0..natoms {
            self.accelerations[i] = forces[i] / universe[i].mass();
            universe[i].add_velocity(0.5 * dt * self.accelerations[i]);
        }
    }
}

/******************************************************************************/
/// Verlet integrator. This one is reversible and symplectic.
pub struct Verlet {
    /// Timestep for the integrator
    timestep: f64,
    /// Previous positions
    prevpos: Vec<Vector3D>,
}

impl Verlet {
    /// Create a new integrator with a timestep of `timestep`.
    pub fn new(timestep: f64) -> Verlet {
        Verlet{
            timestep: timestep,
            prevpos: Vec::new(),
        }
    }
}

impl Integrator for Verlet {
    fn setup(&mut self, universe: &Universe) {
        self.prevpos = vec![Vector3D::new(0.0, 0.0, 0.0); universe.size()];

        let dt = self.timestep;
        // Approximate the positions at t - ∆t
        for (i, part) in universe.iter().enumerate() {
            self.prevpos[i] = *part.position() - *part.velocity() * dt;
        }
    }

    fn integrate(&mut self, universe: &mut Universe) {
        let dt = self.timestep;
        let dt2 = self.timestep * self.timestep;

        let forces = Forces.compute(&universe);
        for (i, p) in universe.iter_mut().enumerate() {
            // Save positions at t
            let tmp = p.position().clone();
            // Update positions at t + ∆t
            let position = 2.0 * tmp - self.prevpos[i] + dt2/p.mass() * forces[i];
            // Update velocities at t
            let velocity = (position - self.prevpos[i]) / (2.0 * dt);
            p.set_position(position);
            p.set_velocity(velocity);
            // Update saved position
            self.prevpos[i] = tmp;
        }
    }
}

/******************************************************************************/
/// Leap-frog integrator. This one is reversible and symplectic.
pub struct LeapFrog {
    /// Timestep for the integrator
    timestep: f64,
    /// Storing the accelerations
    accelerations: Vec<Vector3D>
}

impl LeapFrog {
    /// Create a new integrator with a timestep of `timestep`.
    pub fn new(timestep: f64) -> LeapFrog {
        LeapFrog{
            timestep: timestep,
            accelerations: Vec::new(),
        }
    }
}

impl Integrator for LeapFrog {
    fn setup(&mut self, universe: &Universe) {
        self.accelerations = vec![Vector3D::new(0.0, 0.0, 0.0); universe.size()];
    }

    fn integrate(&mut self, universe: &mut Universe) {
        let dt = self.timestep;
        let dt2 = self.timestep * self.timestep;

        for (i, p) in universe.iter_mut().enumerate() {
            let velocity = p.velocity().clone();
            p.add_position(velocity * dt + 0.5 * self.accelerations[i] * dt2);
        }

        let forces = Forces.compute(&universe);
        for (i, p) in universe.iter_mut().enumerate() {
            let mass = p.mass();
            let acceleration = forces[i]/mass;
            p.add_velocity(0.5 * (self.accelerations[i] + acceleration)* dt);
            self.accelerations[i] = acceleration;
        }
    }
}

/******************************************************************************/
/// This is needed for the BerendsenBarostat implentations. The value comes
/// from the DL_POLY source code.
const WATER_COMPRESSIBILITY: f64 = 7372.0;

/// Berendsen barostat integrator based on velocity-Verlet. This one neither
/// reversible nor symplectic.
pub struct BerendsenBarostat {
    /// Timestep for the integrator
    timestep: f64,
    /// Target pressure for the barostat
    pressure: f64,
    /// Barostat timestep, expressed in units of the main timestep.
    baro_timestep: f64,
    /// Storing the accelerations
    accelerations: Vec<Vector3D>,
    /// Storing the scaling factor
    eta: f64,
}

impl BerendsenBarostat {
    /// Create a new Berendsen barostat with a timestep of `timestep`, and a
    /// target pressure of `pressure`. The barostat timestep is 1000.
    pub fn new(timestep: f64, pressure: f64) -> BerendsenBarostat {
        BerendsenBarostat::with_timestep(timestep, pressure, 1000.0)
    }

    /// Create a new Berendsen barostat with a timestep of `timestep`, a
    /// target pressure of `pressure` and a barostat timestep is of
    /// `baro_timestep`. The barostat timestep is expressed in terms of the main
    /// timestep. With `baro_timestep == 1000` and `timestep == 0.8 fs`, then
    /// the effective barostat timestep will be `800 fs`.
    pub fn with_timestep(timestep: f64, pressure: f64, baro_timestep: f64) -> BerendsenBarostat {
        BerendsenBarostat{
            timestep: timestep,
            pressure: pressure,
            baro_timestep: baro_timestep,
            accelerations: Vec::new(),
            eta: 1.0,
        }
    }

}

impl Integrator for BerendsenBarostat {
    fn setup(&mut self, universe: &Universe) {
        self.accelerations = vec![Vector3D::new(0.0, 0.0, 0.0); universe.size()];
    }

    fn integrate(&mut self, universe: &mut Universe) {
        let dt = self.timestep;

        // Update velocities at t + ∆t/2 and positions at t + ∆t
        for (i, part) in universe.iter_mut().enumerate() {
            part.add_velocity(0.5 * dt * self.accelerations[i]);
            let position = *part.position();
            // Scale all positions
            let position = self.eta * position;
            let position = position + *part.velocity() * dt;
            part.set_position(position);
        }

        universe.cell_mut().scale_mut(self.eta*self.eta*self.eta * Matrix3::one());

        self.eta = f64::cbrt(1.0 - WATER_COMPRESSIBILITY / self.baro_timestep * (self.pressure - universe.pressure()));

        let forces = Forces.compute(&universe);
        // Update accelerations at t + ∆t and velocities at t + ∆t
        for (i, part) in universe.iter_mut().enumerate() {
            self.accelerations[i] = forces[i] / part.mass();
            part.add_velocity(0.5 * dt * self.accelerations[i]);
        }
    }
}
