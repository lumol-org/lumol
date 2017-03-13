// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 Lumol's contributors — BSD license
extern crate lumol;
extern crate lumol_input as input;

use lumol::Logger;
use input::Input;

use std::path::Path;
use std::sync::{Once, ONCE_INIT};
static START: Once = ONCE_INIT;


// This test only run a Monte-Carlo simulation of water, but do not test
// anything for now. It should test the g(r) function someday.
#[test]
fn wolf() {
    START.call_once(|| {Logger::stdout();});
    let path = Path::new(file!()).parent().unwrap().join("data")
                                 .join("mc_water_wolf.toml");

    let mut config = Input::new(path).unwrap().read().unwrap();
    config.simulation.run(&mut config.system, config.nsteps);
}

#[test]
fn ewald() {
    START.call_once(|| {Logger::stdout();});
    let path = Path::new(file!()).parent().unwrap().join("data")
                                 .join("mc_water_ewald.toml");

    let mut config = Input::new(path).unwrap().read().unwrap();
    config.simulation.run(&mut config.system, config.nsteps);
}
