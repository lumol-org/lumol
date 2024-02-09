// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license
use soa_derive::soa_zip;

use lumol_core::{System, Matrix3, Vector3D};

/// The `Integrator` trait define integrator interface for molecular dynamics.
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

/// Velocity-Verlet integrator.
///
/// This integrator is time-reversible and symplectic (volume preserving).
pub struct VelocityVerlet {
    /// Timestep for the integrator
    timestep: f64,
    /// Storing the accelerations
    accelerations: Vec<Vector3D>,
}

impl VelocityVerlet {
    /// Create a new integrator with a timestep of `timestep`.
    pub fn new(timestep: f64) -> VelocityVerlet {
        VelocityVerlet {
            timestep: timestep,
            accelerations: Vec::new(),
        }
    }
}

impl Integrator for VelocityVerlet {
    fn setup(&mut self, system: &System) {
        self.accelerations = vec![Vector3D::zero(); system.size()];
    }

    fn integrate(&mut self, system: &mut System) {
        let dt = self.timestep;

        // Update velocities at t + ∆t/2 and positions at t + ∆t
        for (position, velocity, acceleration) in soa_zip!(
            system.particles_mut(), [mut position, mut velocity], &self.accelerations
        ) {
            *velocity += 0.5 * dt * acceleration;
            *position += velocity * dt;
        }

        let forces = system.forces();
        // Update accelerations at t + ∆t
        for (&mass, acceleration, force) in soa_zip!(
            system.particles(), [mass], &mut self.accelerations, forces
        ) {
            *acceleration = force / mass;
        }

        // Update velocities at t + ∆t
        for (velocity, acceleration) in soa_zip!(
            system.particles_mut(), [mut velocity], &self.accelerations
        ) {
            *velocity += 0.5 * dt * acceleration;
        }
    }
}

/// Verlet integrator.
///
/// This integrator is time-reversible and symplectic (volume preserving).
pub struct Verlet {
    /// Timestep for the integrator
    timestep: f64,
    /// Previous positions
    prevpos: Vec<Vector3D>,
}

impl Verlet {
    /// Create a new integrator with a timestep of `timestep`.
    pub fn new(timestep: f64) -> Verlet {
        Verlet {
            timestep: timestep,
            prevpos: Vec::new(),
        }
    }
}

impl Integrator for Verlet {
    fn setup(&mut self, system: &System) {
        self.prevpos = vec![Vector3D::zero(); system.size()];

        let dt = self.timestep;
        // Approximate the positions at t - ∆t
        for (position, velocity, prevpos) in soa_zip!(
            system.particles(), [position, velocity], &mut self.prevpos
        ) {
            *prevpos = position - velocity * dt;
        }
    }

    fn integrate(&mut self, system: &mut System) {
        let forces = system.forces();
        let dt = self.timestep;
        let dt2 = dt * dt;

        for (position, velocity, mass, prevpos, force) in soa_zip!(
            system.particles_mut(), [mut position, mut velocity, mass], &mut self.prevpos, forces
        ) {
            // Save positions at t
            let tmp = *position;
            // Update positions at t + ∆t
            *position = 2.0 * (*position) - (*prevpos) + dt2 / mass * force;
            // Update velocities at t
            *velocity = ((*position) - (*prevpos)) / (2.0 * dt);
            // Update saved position
            *prevpos = tmp;
        }
    }
}

/// Leap-frog integrator.
///
/// This integrator is time-reversible and symplectic (volume preserving).
pub struct LeapFrog {
    /// Timestep for the integrator
    timestep: f64,
    /// Storing the accelerations
    accelerations: Vec<Vector3D>,
}

impl LeapFrog {
    /// Create a new integrator with a timestep of `timestep`.
    pub fn new(timestep: f64) -> LeapFrog {
        LeapFrog {
            timestep: timestep,
            accelerations: Vec::new(),
        }
    }
}

impl Integrator for LeapFrog {
    fn setup(&mut self, system: &System) {
        self.accelerations = vec![Vector3D::zero(); system.size()];
    }

    fn integrate(&mut self, system: &mut System) {
        let dt = self.timestep;
        let dt2 = dt * dt;

        for (position, velocity, acceleration) in soa_zip!(
            system.particles_mut(), [mut position, velocity], &self.accelerations
        ) {
            *position += velocity * dt + 0.5 * acceleration * dt2;
        }

        let forces = system.forces();
        for (velocity, &mass, acceleration, force) in soa_zip!(
            system.particles_mut(), [mut velocity, mass], &mut self.accelerations, &forces
        ) {
            let new_acceleration = force / mass;
            *velocity += 0.5 * ((*acceleration) + new_acceleration) * dt;
            *acceleration = new_acceleration;
        }
    }
}

/// This is needed for the `BerendsenBarostat` implementation. The value comes
/// from the DL_POLY source code.
const WATER_COMPRESSIBILITY: f64 = 7372.0;

/// Berendsen barostat integrator based on velocity-Verlet.
///
/// This integrator is **neither** time-reversible nor symplectic.
pub struct BerendsenBarostat {
    /// Timestep for the integrator
    timestep: f64,
    /// Target pressure for the barostat
    pressure: f64,
    /// Barostat time scale, expressed in units of the timestep.
    tau: f64,
    /// Storing the accelerations
    accelerations: Vec<Vector3D>,
    /// Storing the scaling factor
    eta: f64,
}

impl BerendsenBarostat {
    /// Create a new Berendsen barostat with an integration timestep of
    /// `timestep`, and a target pressure of `pressure` and the barostat time
    /// scale `tau`.
    pub fn new(timestep: f64, pressure: f64, tau: f64) -> BerendsenBarostat {
        BerendsenBarostat {
            timestep: timestep,
            pressure: pressure,
            tau: tau,
            accelerations: Vec::new(),
            eta: 1.0,
        }
    }
}

impl Integrator for BerendsenBarostat {
    fn setup(&mut self, system: &System) {
        self.accelerations = vec![Vector3D::zero(); system.size()];
    }

    #[allow(clippy::manual_assert)]
    fn integrate(&mut self, system: &mut System) {
        let dt = self.timestep;

        // Update velocities at t + ∆t/2 and positions at t + ∆t
        for (position, velocity, acceleration) in soa_zip!(
            system.particles_mut(), [mut position, mut velocity], &self.accelerations
        ) {
            *velocity += 0.5 * dt * acceleration;
            // Scale all positions
            *position *= self.eta;
            *position += velocity * dt;
        }

        system.cell.scale_mut(self.eta * self.eta * self.eta * Matrix3::one());

        if let Some(maximum_cutoff) = system.maximum_cutoff() {
            if system.cell.lengths().iter().any(|&d| 0.5 * d <= maximum_cutoff) {
                panic!(
                    "Tried to decrease the cell size in Berendesen barostat \
                     but the new size is smaller than the interactions cut off \
                     radius. You can try to increase the cell size or the number \
                     of particles."
                );
            }
        };

        let eta3 = 1.0 - WATER_COMPRESSIBILITY / self.tau * (self.pressure - system.pressure());
        self.eta = f64::cbrt(eta3);

        let forces = system.forces();
        // Update accelerations at t + ∆t and velocities at t + ∆t
        for (velocity, &mass, acceleration, force) in soa_zip!(
            system.particles_mut(), [mut velocity, mass], &mut self.accelerations, &forces
        ) {
            *acceleration = force / mass;
            *velocity += 0.5 * dt * acceleration;
        }
    }
}

/// Anisotropic Berendsen barostat integrator based on velocity-Verlet.
///
/// This integrator is **neither** time-reversible nor symplectic.
pub struct AnisoBerendsenBarostat {
    /// Timestep for the integrator
    timestep: f64,
    /// Target stress matrix for the barostat
    stress: Matrix3,
    /// Barostat time scale, expressed in units of the timestep
    tau: f64,
    /// Storing the accelerations
    accelerations: Vec<Vector3D>,
    /// Storing the scaling factor
    eta: Matrix3,
}

impl AnisoBerendsenBarostat {
    /// Create a new anisotropic Berendsen barostat with an integration timestep
    /// of `timestep`, and a target stress matrix of `stress` and the barostat
    /// time scale `tau`.
    pub fn new(timestep: f64, stress: Matrix3, tau: f64) -> AnisoBerendsenBarostat {
        AnisoBerendsenBarostat {
            timestep: timestep,
            stress: stress,
            tau: tau,
            accelerations: Vec::new(),
            eta: Matrix3::one(),
        }
    }

    /// Create a new anisotropic Berendsen barostat with an integration timestep
    /// of `timestep`, using an hydrostatic stress matrix corresponding to the
    /// pressure `pressure` and the barostat time scale `tau`.
    pub fn hydrostatic(timestep: f64, pressure: f64, tau: f64) -> AnisoBerendsenBarostat {
        AnisoBerendsenBarostat::new(timestep, pressure * Matrix3::one(), tau)
    }
}

impl Integrator for AnisoBerendsenBarostat {
    fn setup(&mut self, system: &System) {
        self.accelerations = vec![Vector3D::zero(); system.size()];
    }

    #[allow(clippy::manual_assert)]
    fn integrate(&mut self, system: &mut System) {
        let dt = self.timestep;

        // Update velocities at t + ∆t/2 and positions at t + ∆t
        for (position, velocity, acceleration) in soa_zip!(
            system.particles_mut(), [mut position, mut velocity], &self.accelerations
        ) {
            *velocity += 0.5 * dt * acceleration;
            // Scale all positions
            *position = self.eta * (*position);
            *position += velocity * dt;
        }

        system.cell.scale_mut(self.eta);

        if let Some(maximum_cutoff) = system.maximum_cutoff() {
            if system.cell.lengths().iter().any(|&d| 0.5 * d <= maximum_cutoff) {
                panic!(
                    "Tried to decrease the cell size in anisotropic Berendesen \
                     barostat but the new size is smaller than the interactions \
                     cut off radius. You can try to increase the cell size or \
                     the number of particles."
                );
            }
        };

        let factor = self.timestep * WATER_COMPRESSIBILITY / self.tau;
        self.eta = Matrix3::one() - factor * (self.stress - system.stress());

        // Make the eta matrix symmetric here
        for i in 0..3 {
            for j in 0..i {
                self.eta[i][j] = 0.5 * (self.eta[i][j] + self.eta[j][i]);
                self.eta[j][i] = self.eta[i][j];
            }
        }

        let forces = system.forces();
        // Update accelerations at t + ∆t and velocities at t + ∆t
        for (velocity, &mass, acceleration, force) in soa_zip!(
            system.particles_mut(), [mut velocity, mass], &mut self.accelerations, &forces
        ) {
            *acceleration = force / mass;
            *velocity += 0.5 * dt * acceleration;
        }
    }
}
