// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

use lumol_core::System;

use crate::velocities;

/// Trait for thermostat algorithms some parameters in a system during a simulation.
pub trait Thermostat {
    /// Function called once at the beginning of the simulation, which allow
    /// for some setup of the thermostat algorithm if needed.
    fn setup(&mut self, _: &System) {}

    /// Main thermostating function. THis should update the system velocities
    /// in some way to produce constant temperature
    fn apply(&mut self, system: &mut System);

    /// Function called once at the end of the simulation.
    fn finish(&mut self, _: &System) {}
}


/// Velocity rescaling thermostat.
///
/// This algorithm controls the temperature by rescaling all the velocities when
/// the temperature differs exceedingly from the desired temperature. A
/// tolerance parameter prevent this algorithm from running too often: if
/// tolerance is 10K and the target temperature is 300K, the algorithm will only
/// run if the instant temperature is below 290K or above 310K.
///
/// **WARNING**: This thermostat does NOT produces a NVT or NPT ensemble. It
/// will not even produce correct average temperature, except if the rescaling
/// is done at every step. It can It can still be usefull in the equilibration
/// of a system at a given temperature before an actual simulation. A good
/// alternative is the CSVR thermostat, which produces correct ensemble.
pub struct RescaleThermostat {
    /// Target temperature
    temperature: f64,
    /// Tolerance in temperature
    tol: f64,
}

impl RescaleThermostat {
    /// Create a new `RescaleThermostat` acting at temperature `temperature`, with a
    /// tolerance of `5% * temperature`.
    pub fn new(temperature: f64) -> RescaleThermostat {
        assert!(temperature >= 0.0, "The temperature must be positive in thermostats.");
        let tol = 0.05 * temperature;
        RescaleThermostat::with_tolerance(temperature, tol)
    }

    /// Create a new `RescaleThermostat` acting at temperature `T`, with a
    /// tolerance of `tol`. For rescaling all the steps, use `tol = 0`.
    pub fn with_tolerance(temperature: f64, tol: f64) -> RescaleThermostat {
        RescaleThermostat {
            temperature: temperature,
            tol: tol,
        }
    }
}

impl Thermostat for RescaleThermostat {
    fn apply(&mut self, system: &mut System) {
        let instant_temperature = system.temperature();
        if f64::abs(instant_temperature - self.temperature) > self.tol {
            velocities::scale(system, self.temperature);
        }
    }
}

/// Berendsen (or weak coupling) thermostat.
///
/// The Berendsen thermostat sets the simulation temperature by exponentially
/// relaxing to a desired temperature. A more complete description of this
/// algorithm can be found in the original article [1].
///
/// **WARNING**: This thermostat does NOT produces a reliable NVT or NPT
/// ensemble (See [2]). While it produces correct average temperature, it does
/// not reproduce the fluctuations of said temperature. It can still be usefull,
/// especialy for the equilibration part of a simulation. Good alternatives
/// include the CSVR or Nosé-Hoover thermostats (not yet implemented in lumol),
/// which produce correct ensembles.
///
/// [1] Berendsen et al. J. Chem Phys 81, 3684 (1984); doi: 10.1063/1.448118
///
/// [2] Braun et al. J. Chem. Theo. Comp. 14, 10 (2018) doi: 10.1021/acs.jctc
pub struct BerendsenThermostat {
    /// Target temperature
    temperature: f64,
    /// Timestep of the thermostat, expressed as a multiplicative factor of the
    /// integrator timestep.
    tau: f64,
}

impl BerendsenThermostat {
    /// Create a new `BerendsenThermostat` acting at temperature `T`, with a
    /// timestep of `tau` times the integrator timestep.
    pub fn new(temperature: f64, tau: f64) -> BerendsenThermostat {
        assert!(temperature >= 0.0, "The temperature must be positive in thermostats.");
        assert!(tau >= 1.0, "The timestep must be larger than 1 in berendsen thermostat.");
        BerendsenThermostat {
            temperature: temperature,
            tau: tau,
        }
    }
}

impl Thermostat for BerendsenThermostat {
    fn apply(&mut self, system: &mut System) {
        let instant_temperature = system.temperature();
        let factor =
            f64::sqrt(1.0 + 1.0 / self.tau * (self.temperature / instant_temperature - 1.0));
        for velocity in system.particles_mut().velocity {
            *velocity *= factor;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use lumol_core::{Vector3D, Particle, Molecule, System, UnitCell};
    use crate::velocities::{BoltzmannVelocities, InitVelocities};

    use approx::assert_ulps_eq;

    fn testing_system() -> System {
        let mut system = System::with_cell(UnitCell::cubic(20.0));

        for i in 0..10 {
            for j in 0..10 {
                for k in 0..10 {
                    let mut particle = Particle::new("Cl");
                    particle.position = Vector3D::new(i as f64 * 2.0, j as f64 * 2.0, k as f64 * 2.0);
                    system.add_molecule(Molecule::new(particle));
                }
            }
        }

        let mut velocities = BoltzmannVelocities::new(300.0);
        velocities.init(&mut system);
        return system;
    }

    #[test]
    fn rescale_thermostat() {
        let mut system = testing_system();
        let temperature = system.temperature();
        assert_ulps_eq!(temperature, 300.0, epsilon = 1e-12);

        let mut thermostat = RescaleThermostat::with_tolerance(250.0, 100.0);
        thermostat.apply(&mut system);
        let temperature = system.temperature();
        assert_ulps_eq!(temperature, 300.0, epsilon = 1e-12);

        let mut thermostat = RescaleThermostat::with_tolerance(250.0, 10.0);
        thermostat.apply(&mut system);
        let temperature = system.temperature();
        assert_ulps_eq!(temperature, 250.0, epsilon = 1e-12);
    }

    #[test]
    #[should_panic]
    fn negative_temperature_rescale() {
        let _ = RescaleThermostat::new(-56.0);
    }

    #[test]
    fn berendsen_thermostat() {
        let mut system = testing_system();
        let temperature = system.temperature();
        assert_ulps_eq!(temperature, 300.0, epsilon = 1e-9);

        let mut thermostat = BerendsenThermostat::new(250.0, 100.0);
        for _ in 0..3000 {
            thermostat.apply(&mut system);
        }
        let temperature = system.temperature();
        assert_ulps_eq!(temperature, 250.0, epsilon = 1e-9);
    }

    #[test]
    #[should_panic]
    fn negative_temperature_berendsen() {
        let _ = BerendsenThermostat::new(-56.0, 1000.0);
    }

    #[test]
    #[should_panic]
    fn negative_timestep_berendsen() {
        let _ = BerendsenThermostat::new(56.0, -2.0);
    }

    #[test]
    #[should_panic]
    fn too_small_timestep_berendsen() {
        let _ = BerendsenThermostat::new(56.0, 0.3);
    }
}
