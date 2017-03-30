// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 Lumol's contributors — BSD license

//! Testing molecular dynamics of butane
extern crate lumol;
extern crate lumol_input as input;
extern crate env_logger;

use input::Input;

use std::path::Path;
use std::sync::{Once, ONCE_INIT};
static START: Once = ONCE_INIT;


#[test]
fn bonds_detection() {
    START.call_once(|| {env_logger::init().unwrap();});
    let path = Path::new(file!()).parent().unwrap()
                                 .join("data")
                                 .join("md-butane")
                                 .join("nve.toml");
    let system = Input::new(path).unwrap().read_system().unwrap();
    assert_eq!(system.molecules().len(), 50);

    for molecule in system.molecules() {
        assert_eq!(molecule.bonds().len(), 3);
        assert_eq!(molecule.angles().len(), 2);
        assert_eq!(molecule.dihedrals().len(), 1);
    }
}

#[test]
fn constant_energy() {
    START.call_once(|| {env_logger::init().unwrap();});
    let path = Path::new(file!()).parent().unwrap()
                                 .join("data")
                                 .join("md-butane")
                                 .join("nve.toml");
    let mut config = Input::new(path).unwrap().read().unwrap();


    let e_initial = config.system.total_energy();
    config.simulation.run(&mut config.system, config.nsteps);
    let e_final = config.system.total_energy();
    assert!(f64::abs((e_initial - e_final)/e_final) < 1e-3);
}
