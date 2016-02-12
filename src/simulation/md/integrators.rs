/* Cymbalum, Molecular Simulation in Rust - Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
 */

use types::{Vector3D, Matrix3};
use simulation::{Compute, Forces};
use system::System;

/// The Integrator trait define integrators interface for Molecular Dynamics.
/// An integrator is an algorithm responsible for propagating the equations of
/// motion in the system.
pub trait Integrator {
    /// Setup the integrator. This function is called once by every simulation
    /// run.
    fn setup(&mut self, _: &System) {}
    /// Integrate the equations of motion. This is called at every step of the
    /// simulation.
    fn integrate(&mut self, system: &mut System);
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
    fn setup(&mut self, system: &System) {
        self.accelerations = vec![Vector3D::new(0.0, 0.0, 0.0); system.size()];
    }

    fn integrate(&mut self, system: &mut System) {
        let dt = self.timestep;

        // Update velocities at t + ∆t/2 and positions at t + ∆t
        for (i, part) in system.iter_mut().enumerate() {
            part.velocity = part.velocity + 0.5 * dt * self.accelerations[i];
            let dr = part.velocity * dt;
            part.position = part.position + dr;
        }

        let forces = Forces.compute(&system);
        // Update accelerations at t + ∆t and velocities at t + ∆t
        for (i, part) in system.iter_mut().enumerate() {
            self.accelerations[i] = forces[i] / part.mass;
            part.velocity = part.velocity + 0.5 * dt * self.accelerations[i];
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
    fn setup(&mut self, system: &System) {
        self.prevpos = vec![Vector3D::new(0.0, 0.0, 0.0); system.size()];

        let dt = self.timestep;
        // Approximate the positions at t - ∆t
        for (i, part) in system.iter().enumerate() {
            self.prevpos[i] = part.position - part.velocity * dt;
        }
    }

    fn integrate(&mut self, system: &mut System) {
        let dt = self.timestep;
        let dt2 = self.timestep * self.timestep;

        let forces = Forces.compute(&system);
        for (i, part) in system.iter_mut().enumerate() {
            // Save positions at t
            let tmp = part.position.clone();
            // Update positions at t + ∆t
            let position = 2.0 * tmp - self.prevpos[i] + dt2/part.mass * forces[i];
            // Update velocities at t
            let velocity = (position - self.prevpos[i]) / (2.0 * dt);

            part.position = position;
            part.velocity = velocity;
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
    fn setup(&mut self, system: &System) {
        self.accelerations = vec![Vector3D::new(0.0, 0.0, 0.0); system.size()];
    }

    fn integrate(&mut self, system: &mut System) {
        let dt = self.timestep;
        let dt2 = self.timestep * self.timestep;

        for (i, part) in system.iter_mut().enumerate() {
            part.position = part.position + part.velocity * dt + 0.5 * self.accelerations[i] * dt2;
        }

        let forces = Forces.compute(&system);
        for (i, part) in system.iter_mut().enumerate() {
            let mass = part.mass;
            let acceleration = forces[i]/mass;
            part.velocity = part.velocity + 0.5 * (self.accelerations[i] + acceleration)* dt;
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
    /// Create a new Berendsen barostat with an integration timestep of
    /// `timestep`, and a target pressure of `pressure`. The barostat timestep
    /// is 1000.
    pub fn new(timestep: f64, pressure: f64) -> BerendsenBarostat {
        BerendsenBarostat{
            timestep: timestep,
            pressure: pressure,
            baro_timestep: 1000.0,
            accelerations: Vec::new(),
            eta: 1.0,
        }
    }

    /// Set the barostat timestep. It should be expressed in terms of the
    /// integration timestep. With a barostat timestep of 1000 and an
    /// integration timestep of `0.8 fs`, the effective barostat timestep will
    /// be `800 fs`.
    pub fn timestep(mut self, timestep: f64) -> BerendsenBarostat {
        self.baro_timestep =  timestep;
        return self;
    }
}

impl Integrator for BerendsenBarostat {
    fn setup(&mut self, system: &System) {
        self.accelerations = vec![Vector3D::new(0.0, 0.0, 0.0); system.size()];
    }

    fn integrate(&mut self, system: &mut System) {
        let dt = self.timestep;

        // Update velocities at t + ∆t/2 and positions at t + ∆t
        for (i, part) in system.iter_mut().enumerate() {
            part.velocity = part.velocity + 0.5 * dt * self.accelerations[i];
            // Scale all positions
            part.position = self.eta * part.position;
            part.position = part.position + part.velocity * dt;
        }

        system.cell_mut().scale_mut(self.eta*self.eta*self.eta * Matrix3::one());

        self.eta = f64::cbrt(1.0 - WATER_COMPRESSIBILITY / self.baro_timestep * (self.pressure - system.pressure()));

        let forces = Forces.compute(&system);
        // Update accelerations at t + ∆t and velocities at t + ∆t
        for (i, part) in system.iter_mut().enumerate() {
            self.accelerations[i] = forces[i] / part.mass;
            part.velocity = part.velocity + 0.5 * dt * self.accelerations[i];
        }
    }
}

/******************************************************************************/
/// Anisotropic Berendsen barostat integrator based on velocity-Verlet. This one
/// neither reversible nor symplectic.
pub struct AnisoBerendsenBarostat {
    /// Timestep for the integrator
    timestep: f64,
    /// Target stress matrix for the barostat
    stress: Matrix3,
    /// Barostat timestep, expressed in units of the main timestep.
    baro_timestep: f64,
    /// Storing the accelerations
    accelerations: Vec<Vector3D>,
    /// Storing the scaling factor
    eta: Matrix3,
}

impl AnisoBerendsenBarostat {
    /// Create a new anisotropic Berendsen barostat with an integration timestep
    /// of `timestep`, and a target stress matrix of `stress`. The barostat
    /// timestep is 1000.
    pub fn new(timestep: f64, stress: Matrix3) -> AnisoBerendsenBarostat {
        AnisoBerendsenBarostat{
            timestep: timestep,
            stress: stress,
            baro_timestep: 1000.0,
            accelerations: Vec::new(),
            eta: Matrix3::one(),
        }
    }

    /// Create a new anisotropic Berendsen barostat with an integration timestep
    /// of `timestep`, using an hydrostatic stress matrix corresponding to the
    /// pressure `pressure`. The barostat timestep is 1000.
    pub fn hydrostatic(timestep: f64, pressure: f64) -> AnisoBerendsenBarostat {
        AnisoBerendsenBarostat::new(timestep, pressure * Matrix3::one())
    }

    /// Set the barostat timestep. It should be expressed in terms of the
    /// integration timestep. With a barostat timestep of 1000 and an
    /// integration timestep of `0.8 fs`, the effective barostat timestep will
    /// be `800 fs`.
    pub fn timestep(mut self, timestep: f64) -> AnisoBerendsenBarostat {
        self.baro_timestep =  timestep;
        return self;
    }
}

impl Integrator for AnisoBerendsenBarostat {
    fn setup(&mut self, system: &System) {
        self.accelerations = vec![Vector3D::new(0.0, 0.0, 0.0); system.size()];
    }

    fn integrate(&mut self, system: &mut System) {
        let dt = self.timestep;

        // Update velocities at t + ∆t/2 and positions at t + ∆t
        for (i, part) in system.iter_mut().enumerate() {
            part.velocity = part.velocity + 0.5 * dt * self.accelerations[i];
            // Scale all positions
            part.position = self.eta * part.position;
            part.position = part.position + part.velocity * dt;
        }

        system.cell_mut().scale_mut(self.eta);

        let factor = self.timestep * WATER_COMPRESSIBILITY / self.baro_timestep;
        self.eta = Matrix3::one() - factor * (self.stress - system.stress());

        // Make the eta matrix symetric here
        for i in 0..3 {
            for j in 0..i {
                self.eta[(i, j)] = 0.5 * (self.eta[(i, j)] + self.eta[(j, i)]);
                self.eta[(j, i)] = self.eta[(i, j)];
            }
        }

        let forces = Forces.compute(&system);
        // Update accelerations at t + ∆t and velocities at t + ∆t
        for (i, part) in system.iter_mut().enumerate() {
            self.accelerations[i] = forces[i] / part.mass;
            part.velocity = part.velocity + 0.5 * dt * self.accelerations[i];
        }
    }
}
