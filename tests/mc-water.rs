// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 Lumol's contributors — BSD license
extern crate lumol;
extern crate lumol_input as input;
extern crate env_logger;

use input::Input;

use std::path::Path;
use std::sync::{Once, ONCE_INIT};
static START: Once = ONCE_INIT;


// This test only run a Monte-Carlo simulation of water, but do not test
// anything for now. It should test the g(r) function someday.
#[test]
fn wolf() {
    START.call_once(|| {env_logger::init().unwrap();});
    let path = Path::new(file!()).parent().unwrap()
                                 .join("data")
                                 .join("mc-water")
                                 .join("nvt-wolf.toml");

    let mut config = Input::new(path).unwrap().read().unwrap();
    config.simulation.run(&mut config.system, config.nsteps);
}

#[test]
fn ewald() {
    START.call_once(|| {env_logger::init().unwrap();});
    let path = Path::new(file!()).parent().unwrap()
                                 .join("data")
                                 .join("mc-water")
                                 .join("nvt-ewald.toml");

    let mut config = Input::new(path).unwrap().read().unwrap();
    config.simulation.run(&mut config.system, config.nsteps);
}
