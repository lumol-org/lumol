/* Cymbalum, Molecular Simulation in Rust - Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
 */

//! While running a simulation, we often want to have control over some
//! simulation parameters: the temperature, the pressure, etc. This is the goal
//! of the control algorithms, all implmenting of the `Control` trait.

use universe::Universe;

/// Trait for controling some parameters in an universe during a simulation.
pub trait Control {
    /// Function called once at the beggining of the simulation, which allow
    /// for some setup of the control algorithm if needed.
    fn setup(&mut self, _: &Universe) {}

    /// Do your job, control algorithm!
    fn control(&mut self, universe: &mut Universe);

    /// Function called once at the end of the simulation.
    fn finish(&mut self, _: &Universe) {}
}


/******************************************************************************/
/// Velocity rescaling thermostat.
///
/// The velocity rescale algorithm controls the temperature by rescaling all
/// the velocities when the temperature differs exceedingly from the desired
/// temperature. A tolerance parameter prevent this algorithm from running too
/// often: if tolerance is 10K and the target temperature is 300K, the algorithm
/// will only run if the instant temperature is below 290K or above 310K.
pub struct RescaleThermostat {
    /// Target temperature
    T: f64,
    /// Tolerance in temperature
    tol: f64,
}

impl RescaleThermostat {
    /// Create a new `RescaleThermostat` acting at temperature `T`, with a
    /// tolerance of `5% T`.
    pub fn new(T: f64) -> RescaleThermostat {
        assert!(T >= 0.0, "The temperature must be positive in thermostats.");
        let tol = 0.05 * T;
        RescaleThermostat::with_tolerance(T, tol)
    }

    /// Create a new `RescaleThermostat` acting at temperature `T`, with a
    /// tolerance of `tol`. For rescaling all the steps, use `tol = 0`.
    pub fn with_tolerance(T: f64, tol: f64) -> RescaleThermostat {
        RescaleThermostat{T: T, tol: tol}
    }
}

impl Control for RescaleThermostat {
    fn control(&mut self, universe: &mut Universe) {
        let T_inst = universe.temperature();
        if f64::abs(T_inst - self.T) > self.tol {
            let factor = f64::sqrt(self.T / T_inst);
            for particle in universe {
                particle.velocity = factor * particle.velocity;
            }
        }
    }
}

/******************************************************************************/
/// Berendsen thermostat.
///
/// The berendsen thermostat sets the simulation temperature by exponentially
/// relaxing to a desired temperature. A more complete description of this
/// algorithm can be found in the original article [1].
///
/// [1] H.J.C. Berendsen, et al. J. Chem Phys 81, 3684 (1984); doi: 10.1063/1.448118
pub struct BerendsenThermostat {
    /// Target temperature
    T: f64,
    /// Timestep of the thermostat, expressed as a multiplicative factor of the
    /// integrator timestep.
    tau: f64,
}

impl BerendsenThermostat {
    /// Create a new `BerendsenThermostat` acting at temperature `T`, with a
    /// timestep of `tau` times the integrator timestep.
    pub fn new(T: f64, tau: f64) -> BerendsenThermostat {
        assert!(T >= 0.0, "The temperature must be positive in thermostats.");
        assert!(tau >= 0.0, "The timestep must be positive in berendsen thermostat.");
        BerendsenThermostat{T: T, tau: tau}
    }
}

impl Control for BerendsenThermostat {
    fn control(&mut self, universe: &mut Universe) {
        let T_inst = universe.temperature();
        let factor = f64::sqrt(1.0 + 1.0 / self.tau * (self.T / T_inst - 1.0));
        for particle in universe {
            particle.velocity = factor * particle.velocity;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use universe::*;
    use types::*;

    fn testing_universe() -> Universe {
        let mut universe = Universe::from_cell(UnitCell::cubic(20.0));;

        for i in 0..10 {
            for j in 0..10 {
                for k in 0..10 {
                    let mut p = Particle::new("Cl");
                    p.position = Vector3D::new(i as f64*2.0, j as f64*2.0, k as f64*2.0);
                    universe.add_particle(p);
                }
            }
        }

        let mut velocities = BoltzmanVelocities::new(300.0);
        velocities.init(&mut universe);
        return universe;
    }

    #[test]
    fn rescale_thermostat() {
        let mut universe = testing_universe();
        let T0 = universe.temperature();
        assert_approx_eq!(T0, 300.0, 1e-12);

        let mut thermostat = RescaleThermostat::with_tolerance(250.0, 100.0);
        thermostat.control(&mut universe);
        let T1 = universe.temperature();
        assert_approx_eq!(T1, 300.0, 1e-12);

        let mut thermostat = RescaleThermostat::with_tolerance(250.0, 10.0);
        thermostat.control(&mut universe);
        let T2 = universe.temperature();
        assert_approx_eq!(T2, 250.0, 1e-12);
    }

    #[test]
    fn berendsen_thermostat() {
        let mut universe = testing_universe();
        let T0 = universe.temperature();
        assert_approx_eq!(T0, 300.0, 1e-12);

        let mut thermostat = BerendsenThermostat::new(250.0, 100.0);
        for _ in 0..1000 {
            thermostat.control(&mut universe);
        }
        let T1 = universe.temperature();
        assert_approx_eq!(T1, 250.0, 1e-2);
    }

    #[test]
    #[should_panic]
    fn negative_temperature_rescale() {
        RescaleThermostat::new(-56.0);
    }

    #[test]
    #[should_panic]
    fn negative_temperature_berendsen() {
        BerendsenThermostat::new(-56.0, 1000.0);
    }
}
