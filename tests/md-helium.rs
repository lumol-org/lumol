// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

//! Testing physical properties of a Lennard-Jones gas of Helium using Molecular
//! dynamics
use lumol::input::Input;
use lumol::consts::K_BOLTZMANN;
use lumol::units;

use std::path::Path;
use std::sync::Once;
static START: Once = Once::new();

mod utils;


#[test]
fn constant_energy_velocity_verlet() {
    START.call_once(::env_logger::init);
    let path = Path::new(file!()).parent()
                                 .unwrap()
                                 .join("data")
                                 .join("md-helium")
                                 .join("nve-velocity-verlet.toml");
    let mut config = Input::new(path).unwrap().read().unwrap();

    let e_initial = config.system.total_energy();
    config.simulation.run(&mut config.system, config.nsteps);
    let e_final = config.system.total_energy();
    assert!(f64::abs((e_initial - e_final) / e_final) < 5e-3);
}


#[test]
fn constant_energy_verlet() {
    START.call_once(::env_logger::init);
    let path = Path::new(file!()).parent()
                                 .unwrap()
                                 .join("data")
                                 .join("md-helium")
                                 .join("nve-verlet.toml");
    let mut config = Input::new(path).unwrap().read().unwrap();

    let e_initial = config.system.total_energy();
    config.simulation.run(&mut config.system, config.nsteps);
    let e_final = config.system.total_energy();
    assert!(f64::abs((e_initial - e_final) / e_final) < 1e-2);
}


#[test]
fn constant_energy_leap_frog() {
    START.call_once(::env_logger::init);
    let path = Path::new(file!()).parent()
                                 .unwrap()
                                 .join("data")
                                 .join("md-helium")
                                 .join("nve-leap-frog.toml");
    let mut config = Input::new(path).unwrap().read().unwrap();

    let e_initial = config.system.total_energy();
    config.simulation.run(&mut config.system, config.nsteps);
    let e_final = config.system.total_energy();
    assert!(f64::abs((e_initial - e_final) / e_final) < 5e-3);
}

#[test]
fn perfect_gas() {
    START.call_once(::env_logger::init);
    let path = Path::new(file!()).parent()
                                 .unwrap()
                                 .join("data")
                                 .join("md-helium")
                                 .join("nve-perfect-gas.toml");
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
    START.call_once(::env_logger::init);
    let path = Path::new(file!()).parent()
                                 .unwrap()
                                 .join("data")
                                 .join("md-helium")
                                 .join("npt-berendsen-barostat.toml");
    let mut config = Input::new(path).unwrap().read().unwrap();

    let collector = utils::Collector::starting_at(4000);
    let temperatures = collector.temperatures();
    let pressures = collector.pressures();

    config.simulation.add_output(Box::new(collector));
    config.simulation.run(&mut config.system, config.nsteps);

    let expected = units::from(5000.0, "bar").unwrap();
    let pressure = crate::utils::mean(pressures);
    assert!(f64::abs(pressure - expected) / expected < 5e-2);

    let expected = units::from(273.0, "K").unwrap();
    let temperature = crate::utils::mean(temperatures);
    assert!(f64::abs(temperature - expected) / expected < 1e-2);
}

#[test]
fn shifted() {
    START.call_once(::env_logger::init);
    let path = Path::new(file!()).parent()
                                 .unwrap()
                                 .join("data")
                                 .join("md-helium")
                                 .join("nve-shifted.toml");
    let mut config = Input::new(path).unwrap().read().unwrap();

    let e_initial = config.system.total_energy();
    config.simulation.run(&mut config.system, config.nsteps);
    let e_final = config.system.total_energy();
    assert!(f64::abs((e_initial - e_final) / e_final) < 2e-3);
}


#[test]
fn table_computation() {
    START.call_once(::env_logger::init);
    let path = Path::new(file!()).parent()
                                 .unwrap()
                                 .join("data")
                                 .join("md-helium")
                                 .join("nve-table.toml");
    let mut config = Input::new(path).unwrap().read().unwrap();

    let e_initial = config.system.total_energy();
    config.simulation.run(&mut config.system, config.nsteps);
    let e_final = config.system.total_energy();
    assert!(f64::abs((e_initial - e_final) / e_final) < 5e-3);
}
