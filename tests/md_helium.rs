// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 Lumol's contributors — BSD license

//! Testing physical properties of a Lennard-Jones gas of Helium using Molecular
//! dynamics
extern crate lumol;
extern crate lumol_input as input;

use lumol::Logger;
use lumol::sys::{System, Trajectory, UnitCell};
use lumol::sys::veloc::{BoltzmannVelocities, InitVelocities};
use lumol::energy::{LennardJones, PairInteraction};
use lumol::energy::TableComputation;
use lumol::sim::{Simulation, MolecularDynamics};
use lumol::sim::md::{LeapFrog, VelocityVerlet, Verlet, BerendsenBarostat};
use lumol::consts::K_BOLTZMANN;
use lumol::units;

use input::Input;

use std::sync::{Once, ONCE_INIT};
static START: Once = ONCE_INIT;

use std::path::Path;

fn get_system() -> System {
    let data_dir = Path::new(file!()).parent().unwrap();
    let configuration = data_dir.join("data").join("helium.xyz");
    let mut system = Trajectory::open(configuration)
                                .and_then(|mut traj| traj.read())
                                .unwrap();
    system.set_cell(UnitCell::cubic(10.0));

    let mut velocities = BoltzmannVelocities::new(units::from(300.0, "K").unwrap());
    velocities.init(&mut system);
    return system;
}

fn get_system_with_interaction() -> System {
    let mut system = get_system();
    let lj = Box::new(LennardJones{
        sigma: units::from(2.0, "A").unwrap(),
        epsilon: units::from(0.2, "kJ/mol").unwrap()
    });

    system.interactions_mut().add_pair("He", "He", PairInteraction::new(lj, 10.0));
    return system;
}

#[test]
fn constant_energy_velocity_verlet() {
    START.call_once(|| {Logger::stdout();});

    let path = Path::new(file!()).parent().unwrap().join("data")
                                 .join("md_helium/velocity_verlet.toml");
    let mut config = Input::new(path).unwrap().read().unwrap();

    let e_initial = config.system.total_energy();
    config.simulation.run(&mut config.system, config.nsteps);
    let e_final = config.system.total_energy();
    assert!(f64::abs((e_initial - e_final)/e_final) < 1e-3);
}


#[test]
fn constant_energy_verlet() {
    START.call_once(|| {Logger::stdout();});

    let path = Path::new(file!()).parent().unwrap().join("data")
                                 .join("md_helium/verlet.toml");
    let mut config = Input::new(path).unwrap().read().unwrap();

    let e_initial = config.system.total_energy();
    config.simulation.run(&mut config.system, config.nsteps);
    let e_final = config.system.total_energy();
    assert!(f64::abs((e_initial - e_final)/e_final) < 1e-2);
}


#[test]
fn constant_energy_leap_frog() {
    START.call_once(|| {Logger::stdout();});

    let path = Path::new(file!()).parent().unwrap().join("data")
                                 .join("md_helium/leap_frog.toml");
    let mut config = Input::new(path).unwrap().read().unwrap();

    let e_initial = config.system.total_energy();
    config.simulation.run(&mut config.system, config.nsteps);
    let e_final = config.system.total_energy();
    assert!(f64::abs((e_initial - e_final)/e_final) < 1e-3);
}

#[test]
fn perfect_gas() {
     START.call_once(|| {Logger::stdout();});

    let path = Path::new(file!()).parent().unwrap().join("data")
                                 .join("md_helium/perfect_gas.toml");
    let mut config = Input::new(path).unwrap().read().unwrap();

    config.simulation.run(&mut config.system, config.nsteps);

    let pressure = config.system.pressure();
    let volume = config.system.volume();
    let temperature = config.system.temperature();
    let n = config.system.size() as f64;

    assert!(f64::abs(pressure * volume - n * K_BOLTZMANN * temperature) < 1e-3);
}

#[test]
fn berendsen_barostat() {
    START.call_once(|| {Logger::stdout();});
      let path = Path::new(file!()).parent().unwrap().join("data")
                                 .join("md_helium/berendsen_barostat.toml");
    let mut config = Input::new(path).unwrap().read().unwrap();

    config.simulation.run(&mut config.system, config.nsteps);

    let pressure = units::from(5000.0, "bar").unwrap();
    assert!(f64::abs((config.system.pressure() - pressure)/pressure) < 5e-2);
}

#[test]
fn shifted() {
    START.call_once(|| {Logger::stdout();});

    let path = Path::new(file!()).parent().unwrap().join("data")
                                 .join("md_helium/shifted.toml");
    let mut config = Input::new(path).unwrap().read().unwrap();

    let e_initial = config.system.total_energy();
    config.simulation.run(&mut config.system, config.nsteps);
    let e_final = config.system.total_energy();
    assert!(f64::abs((e_initial - e_final)/e_final) < 2e-3);
}


#[test]
fn table_computation() {
     START.call_once(|| {Logger::stdout();});

    let path = Path::new(file!()).parent().unwrap().join("data")
                                 .join("md_helium/table_computation.toml");
    let mut config = Input::new(path).unwrap().read().unwrap();

    let e_initial = config.system.total_energy();
    config.simulation.run(&mut config.system, config.nsteps);
    let e_final = config.system.total_energy();
    assert!(f64::abs((e_initial - e_final)/e_final) < 2e-3);
}
