// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

//! Testing physical properties of f-SPC water
use lumol::input::Input;

use std::path::Path;
use std::sync::Once;
static START: Once = Once::new();


#[test]
fn constant_energy_ewald() {
    START.call_once(::env_logger::init);
    let path = Path::new(file!()).parent()
                                 .unwrap()
                                 .join("data")
                                 .join("md-water")
                                 .join("nve-ewald.toml");
    let mut config = Input::new(path).unwrap().read().unwrap();

    let e_initial = config.system.total_energy();
    config.simulation.run(&mut config.system, config.nsteps);
    let e_final = config.system.total_energy();

    // FIXME? This thresold is really bad
    assert!(f64::abs((e_initial - e_final) / e_final) < 1e-1);
}

#[test]
fn constant_energy_wolf() {
    START.call_once(::env_logger::init);
    let path = Path::new(file!()).parent()
                                 .unwrap()
                                 .join("data")
                                 .join("md-water")
                                 .join("nve-wolf.toml");
    let mut config = Input::new(path).unwrap().read().unwrap();

    let e_initial = config.system.total_energy();
    config.simulation.run(&mut config.system, config.nsteps);
    let e_final = config.system.total_energy();
    assert!(f64::abs((e_initial - e_final) / e_final) < 3e-2);
}
