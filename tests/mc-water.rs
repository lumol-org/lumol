// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license
use lumol::input::Input;
use lumol::units;

use std::path::Path;
use std::sync::Once;
static START: Once = Once::new();

mod utils;

// This test only run a Monte Carlo simulation of water, but do not test
// anything for now. It should test the g(r) function someday.
#[test]
fn wolf_nvt() {
    START.call_once(::env_logger::init);
    let path = Path::new(file!()).parent()
                                 .unwrap()
                                 .join("data")
                                 .join("mc-water")
                                 .join("nvt-wolf.toml");

    let mut config = Input::new(path).unwrap().read().unwrap();
    config.simulation.run(&mut config.system, config.nsteps);
}

#[test]
fn wolf_npt() {
    START.call_once(::env_logger::init);
    let path = Path::new(file!()).parent()
                                 .unwrap()
                                 .join("data")
                                 .join("mc-water")
                                 .join("npt-wolf.toml");

    let mut config = Input::new(path).unwrap().read().unwrap();

    let collector = utils::Collector::starting_at((config.nsteps - 5_000) as u64);
    let pressures = collector.pressures();

    config.simulation.add_output(Box::new(collector));
    config.simulation.run(&mut config.system, config.nsteps);

    let pressure = utils::mean(pressures);
    let expected = units::from(1000.0, "bar").unwrap();
    let tolerance = units::from(800.0, "bar").unwrap();

    assert!(f64::abs(pressure - expected) < tolerance);
}

// This test only run a Monte Carlo simulation of water, but do not test
// anything for now. It should test the g(r) function someday.
#[test]
fn ewald_nvt() {
    START.call_once(::env_logger::init);
    let path = Path::new(file!()).parent()
                                 .unwrap()
                                 .join("data")
                                 .join("mc-water")
                                 .join("nvt-ewald.toml");

    let mut config = Input::new(path).unwrap().read().unwrap();
    config.simulation.run(&mut config.system, config.nsteps);
}


#[test]
fn ewald_npt() {
    START.call_once(::env_logger::init);
    let path = Path::new(file!()).parent()
                                 .unwrap()
                                 .join("data")
                                 .join("mc-water")
                                 .join("npt-ewald.toml");

    let mut config = Input::new(path).unwrap().read().unwrap();

    let collector = utils::Collector::starting_at((config.nsteps - 5_000) as u64);
    let pressures = collector.pressures();

    config.simulation.add_output(Box::new(collector));
    config.simulation.run(&mut config.system, config.nsteps);

    let pressure = utils::mean(pressures);
    let expected = units::from(1000.0, "bar").unwrap();
    let tolerance = units::from(800.0, "bar").unwrap();

    assert!(f64::abs(pressure - expected) < tolerance);
}
