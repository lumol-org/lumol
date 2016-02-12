/*
 * Cymbalum, Molecular Simulation in Rust
 * Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
*/
//! Testing physical properties of a Lennard-Jones gaz of Helium using Molecular
//! dynamics
extern crate cymbalum;
use self::cymbalum::*;

use std::sync::{Once, ONCE_INIT};
static START: Once = ONCE_INIT;

use std::path::Path;

fn get_system() -> System {
    let data_dir = Path::new(file!()).parent().unwrap();
    let configuration = data_dir.join("data").join("helium.xyz");
    let mut system = System::from_file(configuration.to_str().unwrap()).unwrap();
    system.set_cell(UnitCell::cubic(10.0));

    let mut velocities = BoltzmanVelocities::new(units::from(300.0, "K").unwrap());
    velocities.init(&mut system);
    return system;
}

fn get_system_with_interaction() -> System {
    let mut system = get_system();
    system.add_pair_interaction("He", "He",
        Box::new(LennardJones{
            sigma: units::from(2.0, "A").unwrap(),
            epsilon: units::from(0.2, "kJ/mol").unwrap()
        })
    );

    return system;
}

#[test]
fn constant_energy_velocity_verlet() {
    START.call_once(|| {Logger::stdout();});
    let mut system = get_system_with_interaction();
    let mut simulation = Simulation::new(
        MolecularDynamics::from_integrator(
            VelocityVerlet::new(units::from(1.0, "fs").unwrap())
        )
    );
    let e_initial = system.total_energy();
    simulation.run(&mut system, 1000);
    let e_final = system.total_energy();
    assert!(f64::abs((e_initial - e_final)/e_final) < 1e-3);
}


#[test]
fn constant_energy_verlet() {
    START.call_once(|| {Logger::stdout();});
    let mut system = get_system_with_interaction();
    let mut simulation = Simulation::new(
        MolecularDynamics::from_integrator(
            Verlet::new(units::from(1.0, "fs").unwrap())
        )
    );
    let e_initial = system.total_energy();
    simulation.run(&mut system, 1000);
    let e_final = system.total_energy();
    assert!(f64::abs((e_initial - e_final)/e_final) < 1e-2);
}


#[test]
fn constant_energy_leap_frog() {
    START.call_once(|| {Logger::stdout();});
    let mut system = get_system_with_interaction();
    let mut simulation = Simulation::new(
        MolecularDynamics::from_integrator(
            LeapFrog::new(units::from(1.0, "fs").unwrap())
        )
    );
    let e_initial = system.total_energy();
    simulation.run(&mut system, 1000);
    let e_final = system.total_energy();
    assert!(f64::abs((e_initial - e_final)/e_final) < 1e-3);
}

#[test]
fn perfect_gaz() {
    START.call_once(|| {Logger::stdout();});
    let mut system = get_system_with_interaction();
    let mut simulation = Simulation::new(
        MolecularDynamics::from_integrator(
            VelocityVerlet::new(units::from(1.0, "fs").unwrap())
        )
    );

    // dilating the system!
    for particle in system.iter_mut() {
        particle.position = 10.0 * particle.position;
    }
    system.set_cell(UnitCell::cubic(100.0));

    simulation.run(&mut system, 1000);
    let pressure = system.pressure();
    let volume = system.volume();
    let temperature = system.temperature();
    let n = system.size() as f64;

    assert!(f64::abs(pressure * volume - n * constants::K_BOLTZMANN * temperature) < 1e-3);
}

#[test]
fn berendsen_barostat() {
    START.call_once(|| {Logger::stdout();});
    let mut system = get_system_with_interaction();
    let mut simulation = Simulation::new(
        MolecularDynamics::from_integrator(
            BerendsenBarostat::new(
                units::from(1.0, "fs").unwrap(),
                units::from(5000.0, "bar").unwrap()
            )
        )
    );

    simulation.run(&mut system, 1000);
    let pressure = units::from(5000.0, "bar").unwrap();
    assert!(f64::abs((system.pressure() - pressure)/pressure) < 5e-2);
}

#[test]
fn cutoff_computation() {
    START.call_once(|| {Logger::stdout();});
    let mut system = get_system();

    system.add_pair_interaction("He", "He",
        Box::new(CutoffComputation::new(
            Box::new(LennardJones{
                sigma: units::from(2.0, "A").unwrap(),
                epsilon: units::from(0.2, "kJ/mol").unwrap()
            }),
            units::from(7.0, "A").unwrap()
        ))
    );

    let mut simulation = Simulation::new(
        MolecularDynamics::from_integrator(
            VelocityVerlet::new(units::from(1.0, "fs").unwrap())
        )
    );
    simulation.run(&mut system, 100);

    let e_initial = system.total_energy();
    simulation.run(&mut system, 1000);
    let e_final = system.total_energy();
    assert!(f64::abs((e_initial - e_final)/e_final) < 2e-3);
}


#[test]
fn table_computation() {
    START.call_once(|| {Logger::stdout();});
    let mut system = get_system();

    system.add_pair_interaction("He", "He",
        Box::new(TableComputation::new(
            Box::new(LennardJones{
                sigma: units::from(2.0, "A").unwrap(),
                epsilon: units::from(0.2, "kJ/mol").unwrap()
            }),
            1000,
            units::from(7.0, "A").unwrap()
        ))
    );

    let mut simulation = Simulation::new(
        MolecularDynamics::from_integrator(
            VelocityVerlet::new(units::from(1.0, "fs").unwrap())
        )
    );
    simulation.run(&mut system, 100);

    let e_initial = system.total_energy();
    simulation.run(&mut system, 1000);
    let e_final = system.total_energy();
    assert!(f64::abs((e_initial - e_final)/e_final) < 2e-3);
}
