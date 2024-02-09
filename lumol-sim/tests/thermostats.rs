// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

use lumol_core::{Vector3D, Particle, Molecule, System, UnitCell};
use lumol_core::consts::K_BOLTZMANN;

use lumol_sim::{BoltzmannVelocities, InitVelocities};
use lumol_sim::md::{Integrator, VelocityVerlet};
use lumol_sim::md::{Thermostat, RescaleThermostat, BerendsenThermostat, CSVRThermostat};

use approx::{assert_ulps_eq, assert_relative_eq};

// An ideal gas system
fn testing_system() -> System {
    let mut system = System::with_cell(UnitCell::cubic(20.0));

    for i in 0..10 {
        for j in 0..10 {
            for k in 0..10 {
                let mut particle = Particle::new("He");
                particle.position = Vector3D::new(i as f64 * 2.0, j as f64 * 2.0, k as f64 * 2.0);
                system.add_molecule(Molecule::new(particle));
            }
        }
    }

    let mut velocities = BoltzmannVelocities::new(300.0);
    velocities.init(&mut system);

    assert_ulps_eq!(system.temperature(), 300.0, epsilon = 1e-9);

    system
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
fn berendsen_thermostat() {
    let mut system = testing_system();

    let mut thermostat = BerendsenThermostat::new(250.0, 10.0);
    let mut integrator = VelocityVerlet::new(1.0);
    integrator.setup(&system);
    thermostat.setup(&system);
    // equilibrate
    for _ in 0..100 {
        integrator.integrate(&mut system);
        thermostat.apply(&mut system);
    }

    // accumulate
    let mut temperatures = Vec::new();
    for _ in 0..100 {
        integrator.integrate(&mut system);
        thermostat.apply(&mut system);
        temperatures.push(system.temperature());
    }

    let mean = temperatures.iter().sum::<f64>() / temperatures.len() as f64;
    assert_relative_eq!(mean, 250.0, epsilon=1e-3);
}


#[test]
fn csvr_thermostat() {
    let mut system = testing_system();
    let temperature = system.temperature();
    assert_ulps_eq!(temperature, 300.0, epsilon = 1e-9);


    let mut thermostat = CSVRThermostat::new(250.0, 10.0);
    let mut integrator = VelocityVerlet::new(1.0);
    integrator.setup(&system);
    thermostat.setup(&system);

    // equilibrate
    for _ in 0..100 {
        integrator.integrate(&mut system);
        thermostat.apply(&mut system);
    }

    // accumulate
    let mut kinetic = Vec::new();
    for _ in 0..100 {
        integrator.integrate(&mut system);
        thermostat.apply(&mut system);
        kinetic.push(system.kinetic_energy());
    }

    let temperature = 250.0;
    let dof = system.degrees_of_freedom() as f64;
    let mean = kinetic.iter().sum::<f64>() / kinetic.len() as f64;
    let expected = dof * K_BOLTZMANN * temperature / 2.0;
    assert_relative_eq!(mean, expected, epsilon=1e-3);

    let variance = kinetic.iter()
        .map(|t| (mean - t) * (mean - t))
        .sum::<f64>() / kinetic.len() as f64;

    // expected variance
    let expected = dof * (K_BOLTZMANN * temperature) * (K_BOLTZMANN * temperature) / 2.0;
    assert_relative_eq!(variance, expected, epsilon=1e-3);
}
