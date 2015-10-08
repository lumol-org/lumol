/*
 * Cymbalum, Molecular Simulation in Rust
 * Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
*/
#![allow(non_snake_case)]
//! Testing physical properties of a Lennard-Jones gaz of Helium.
extern crate cymbalum;
use self::cymbalum::*;

use std::sync::{Once, ONCE_INIT};
static START: Once = ONCE_INIT;

use std::path::Path;

fn get_universe() -> Universe {
    let data_dir = Path::new(file!()).parent().unwrap();
    let configuration = data_dir.join("data").join("helium.xyz");
    let mut universe = Universe::from_file(configuration.to_str().unwrap()).unwrap();
    universe.set_cell(UnitCell::cubic(10.0));

    let mut velocities = BoltzmanVelocities::new(units::from(300.0, "K").unwrap());
    velocities.init(&mut universe);
    return universe;
}

fn get_universe_with_interaction() -> Universe {
    let mut universe = get_universe();
    universe.add_pair_interaction("He", "He",
        Box::new(LennardJones{
            sigma: units::from(2.0, "A").unwrap(),
            epsilon: units::from(0.2, "kJ/mol").unwrap()
        })
    );

    return universe;
}

#[test]
fn constant_energy_velocity_verlet() {
    START.call_once(|| {Logger::stdout();});
    let mut universe = get_universe_with_interaction();
    let mut simulation = Simulation::new(
        MolecularDynamics::from_integrator(
            VelocityVerlet::new(units::from(1.0, "fs").unwrap())
        )
    );
    let E_initial = universe.total_energy();
    simulation.run(&mut universe, 1000);
    let E_final = universe.total_energy();
    assert!(f64::abs(E_initial - E_final)/E_final < 1e-3);
}


#[test]
fn constant_energy_verlet() {
    START.call_once(|| {Logger::stdout();});
    let mut universe = get_universe_with_interaction();
    let mut simulation = Simulation::new(
        MolecularDynamics::from_integrator(
            Verlet::new(units::from(1.0, "fs").unwrap())
        )
    );
    let E_initial = universe.total_energy();
    simulation.run(&mut universe, 1000);
    let E_final = universe.total_energy();
    assert!(f64::abs(E_initial - E_final)/E_final < 1e-2);
}


#[test]
fn constant_energy_leap_frog() {
    START.call_once(|| {Logger::stdout();});
    let mut universe = get_universe_with_interaction();
    let mut simulation = Simulation::new(
        MolecularDynamics::from_integrator(
            LeapFrog::new(units::from(1.0, "fs").unwrap())
        )
    );
    let E_initial = universe.total_energy();
    simulation.run(&mut universe, 1000);
    let E_final = universe.total_energy();
    assert!(f64::abs(E_initial - E_final)/E_final < 1e-3);
}

#[test]
fn perfect_gaz() {
    START.call_once(|| {Logger::stdout();});
    let mut universe = get_universe_with_interaction();
    let mut simulation = Simulation::new(
        MolecularDynamics::from_integrator(
            VelocityVerlet::new(units::from(1.0, "fs").unwrap())
        )
    );

    // dilating the universe!
    for particle in universe.iter_mut() {
        let pos = particle.position().clone();
        particle.set_position(10.0 * pos);
    }
    universe.set_cell(UnitCell::cubic(100.0));

    simulation.run(&mut universe, 1000);
    let P = universe.pressure();
    let V = universe.volume();
    let T = universe.temperature();
    let N = universe.size() as f64;

    assert!(f64::abs(P * V - N * constants::K_BOLTZMANN * T) < 1e-3);
}

#[test]
fn berendsen_barostat() {
    START.call_once(|| {Logger::stdout();});
    let mut universe = get_universe_with_interaction();
    let mut md = MolecularDynamics::from_integrator(
        BerendsenBarostat::with_timestep(
            units::from(1.0, "fs").unwrap(),
            units::from(5000.0, "bar").unwrap(),
            100.0
        )
    );
    md.add_control(BerendsenThermostat::new(units::from(300.0, "K").unwrap(), 10.0));
    let mut simulation = Simulation::new(md);

    simulation.run(&mut universe, 5000);
    let P = universe.pressure();
    let T = universe.temperature();

    assert!(f64::abs(P - units::from(5000.0, "bar").unwrap())/P < 8e-2);
    assert!(f64::abs(T - units::from(300.0, "K").unwrap())/T < 1e-2);
}

#[test]
fn cutoff_computation() {
    START.call_once(|| {Logger::stdout();});
    let mut universe = get_universe();

    universe.add_pair_interaction("He", "He",
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
    simulation.run(&mut universe, 100);

    let E_initial = universe.total_energy();
    simulation.run(&mut universe, 1000);
    let E_final = universe.total_energy();
    assert!(f64::abs(E_initial - E_final)/E_final < 1e-3);
}


#[test]
fn table_computation() {
    START.call_once(|| {Logger::stdout();});
    let mut universe = get_universe();

    universe.add_pair_interaction("He", "He",
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
    simulation.run(&mut universe, 100);

    let E_initial = universe.total_energy();
    simulation.run(&mut universe, 1000);
    let E_final = universe.total_energy();
    assert!(f64::abs(E_initial - E_final)/E_final < 1e-3);
}
