// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

use lumol_core::System;
use lumol_core::consts::K_BOLTZMANN;

use rand::{self, SeedableRng};
use rand_distr::{Distribution, Normal, Gamma};

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
    /// Create a new `BerendsenThermostat` acting at the given `temperature`,
    /// with a timestep of `tau` times the integrator timestep.
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
        let factor = f64::sqrt(1.0 + (self.temperature / instant_temperature - 1.0) / self.tau);
        for velocity in system.particles_mut().velocity {
            *velocity *= factor;
        }
    }
}

/// A thermostat using the Canonical Sampling Through Velocities Rescaling
/// algorithm.
///
/// This algorithm works by combining a Berendsen-type thermostat that ensure a
/// fast relaxation to the equilibrium temperature, and a stochastic thermostat
/// that ensure that the NVT distribution of kinetic enegies is actually
/// sampled.
///
/// For a more in-depth description of the algorithm, see [1].
///
/// [1] Bussi et al. J. Chem. Phys. 126, 014101 (2007) doi: 10.1063/1.2408420
pub struct CSVRThermostat {
    /// Target kinetic energy for the system, per degree of freedom
    target_kinetic_per_dof: f64,
    /// Timestep of the thermostat, expressed as a multiplicative factor of the
    /// integrator timestep.
    tau: f64,
    /// Random number generator for the stochatsic propagation of kinetic energy
    rng: Box<dyn rand::RngCore>,
    /// normal (i.e. gaussian) distribution
    normal: Normal<f64>,
}

impl CSVRThermostat {
    /// Create a new `CSVRThermostat` enforcing the given `temperature`, with a
    /// timestep of `tau` times the integrator timestep.
    pub fn new(temperature: f64, tau: f64) -> CSVRThermostat {
        let rng = Box::new(rand_xorshift::XorShiftRng::from_seed([
            0xeb, 0xa8, 0xe4, 0x29, 0xca, 0x60, 0x44, 0xb0,
            0xd3, 0x77, 0xc6, 0xa0, 0x21, 0x71, 0x37, 0xf7,
        ]));
        return CSVRThermostat::from_rng(temperature, tau, rng)
    }

    /// Create a new `CSVRThermostat` enforcing the given `temperature`, with a
    /// timestep of `tau` times the integrator timestep, using the given `rng`
    /// when generating random noise.
    pub fn from_rng(temperature: f64, tau: f64, rng: Box<dyn rand::RngCore>) -> CSVRThermostat {
        assert!(temperature >= 0.0, "The temperature must be positive in thermostats.");
        assert!(tau >= 1.0, "The timestep must be larger than 1 in CSVR thermostat.");
        CSVRThermostat {
            target_kinetic_per_dof: K_BOLTZMANN * temperature / 2.0,
            tau: tau,
            rng: rng,
            normal: Normal::new(0.0, 1.0).expect("bad normal distribution"),
        }
    }

    /// Get the sum of n independent gaussian noises squared, i.e. the Wiener
    /// noise in equation 4 of Bussi2007.
    ///
    /// This functions returns two values: a single gaussian noise and the sum
    /// of n squared gaussian noises.
    fn sum_noises(&mut self, n: usize) -> (f64, f64) {
        if n == 0 {
            return (0.0, 0.0);
        } else if n == 1 {
            let rr = self.normal.sample(&mut self.rng);
            return (rr, rr * rr);
        }

        let rr = self.normal.sample(&mut self.rng);
        if n % 2 == 0 {
            let gamma = Gamma::new((n / 2) as f64, 1.0).expect("bad gamma distribution");
            return (rr, 2.0 * gamma.sample(&mut self.rng));
        } else {
            let gamma = Gamma::new(((n - 1) / 2) as f64, 1.0).expect("bad gamma distribution");
            let wiener = 2.0 * gamma.sample(&mut self.rng) + rr * rr;
            return (rr, wiener);
        }
    }
}

impl Thermostat for CSVRThermostat {
    fn apply(&mut self, system: &mut System) {
        let kinetic = system.kinetic_energy();
        let kinetic_factor = self.target_kinetic_per_dof / kinetic;
        let exp_1 = f64::exp(-1.0/self.tau);
        let exp_2 = (1.0 - exp_1) * kinetic_factor;

        let dof = system.degrees_of_freedom();
        let (gauss, wiener) = self.sum_noises(dof - 1);

        let scale = exp_1 + exp_2 * (gauss * gauss + wiener) + 2.0 * gauss * f64::sqrt(exp_1 * exp_2);
        let alpha = f64::sqrt(scale);
        for velocity in system.particles_mut().velocity {
            *velocity *= alpha;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // The actual thermostat code is tested in lumol-sim/tests/thermostats.rs

    #[test]
    #[should_panic(expected="The temperature must be positive in thermostats.")]
    fn negative_temperature_rescale() {
        let _ = RescaleThermostat::new(-56.0);
    }

    #[test]
    #[should_panic(expected="The temperature must be positive in thermostats.")]
    fn negative_temperature_berendsen() {
        let _ = BerendsenThermostat::new(-56.0, 1000.0);
    }

    #[test]
    #[should_panic(expected="The timestep must be larger than 1 in berendsen thermostat.")]
    fn negative_timestep_berendsen() {
        let _ = BerendsenThermostat::new(56.0, -2.0);
    }

    #[test]
    #[should_panic(expected="The timestep must be larger than 1 in berendsen thermostat.")]
    fn too_small_timestep_berendsen() {
        let _ = BerendsenThermostat::new(56.0, 0.3);
    }

    #[test]
    #[should_panic(expected="The temperature must be positive in thermostats.")]
    fn negative_temperature_csvr() {
        let _ = CSVRThermostat::new(-56.0, 1000.0);
    }

    #[test]
    #[should_panic(expected="The timestep must be larger than 1 in CSVR thermostat.")]
    fn negative_timestep_csvr() {
        let _ = CSVRThermostat::new(56.0, -2.0);
    }

    #[test]
    #[should_panic(expected="The timestep must be larger than 1 in CSVR thermostat.")]
    fn too_small_timestep_csvr() {
        let _ = CSVRThermostat::new(56.0, 0.3);
    }
}
