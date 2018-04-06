// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

//! Testing physical properties of a Lennard-Jones gaz of Helium using
//! Monte Carlo simulation
extern crate env_logger;
extern crate lumol;
extern crate lumol_input as input;

use lumol::consts::K_BOLTZMANN;
use lumol::units;

use input::Input;

use std::path::Path;
use std::sync::{Once, ONCE_INIT};
static START: Once = ONCE_INIT;


#[test]
fn perfect_gas() {
    START.call_once(::env_logger::init);
    let path = Path::new(file!()).parent()
                                 .unwrap()
                                 .join("data")
                                 .join("mc-helium")
                                 .join("nvt.toml");

    let mut config = Input::new(path).unwrap().read().unwrap();
    config.simulation.run(&mut config.system, config.nsteps);
    let pressure = config.system.pressure();
    let volume = config.system.volume();
    let temperature = units::from(300.0, "K").unwrap();

    let pv = pressure * volume;
    let nkt = config.system.size() as f64 * K_BOLTZMANN * temperature;
    let msg = format!("{} {}", f64::abs(pv - nkt), f64::abs(pv - nkt) / pv);
    assert!(f64::abs(pv - nkt) / pv < 2e-2, msg);
}
