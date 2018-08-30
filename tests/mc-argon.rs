// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license
extern crate env_logger;
extern crate lumol;
extern crate lumol_input as input;

use input::Input;
use lumol::units;

use std::path::Path;
use std::sync::{Once, ONCE_INIT};
static START: Once = ONCE_INIT;

mod utils;

#[test]
fn npt() {
    START.call_once(::env_logger::init);
    let path = Path::new(file!()).parent()
                                 .unwrap()
                                 .join("data")
                                 .join("mc-argon")
                                 .join("npt.toml");

    let mut config = Input::new(path).unwrap().read().unwrap();

    let collecter = utils::Collecter::new(500);
    let pressures = collecter.pressures();

    config.simulation.add_output(Box::new(collecter));
    config.simulation.run(&mut config.system, config.nsteps);

    let expected = units::from(200.0, "bar").unwrap();
    let pressure = ::utils::mean(pressures.clone());
    assert!(f64::abs(pressure - expected) / expected < 1e-2);
}
