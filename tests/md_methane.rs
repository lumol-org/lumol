// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 Lumol's contributors — BSD license

//! Testing molecular dynamics of methane
extern crate lumol;
extern crate lumol_input as input;

use lumol::Logger;

use input::Input;

use std::path::Path;
use std::sync::{Once, ONCE_INIT};
static START: Once = ONCE_INIT;



#[test]
fn bonds_detection() {
    START.call_once(|| {Logger::stdout();});
    let path = Path::new(file!()).parent().unwrap().join("data")
                                     .join("md_methane.toml");
    let system = Input::new(path).unwrap().read_system().unwrap();
    assert_eq!(system.molecules().len(), 150);

    for molecule in system.molecules() {
        assert_eq!(molecule.bonds().len(), 4);
        assert_eq!(molecule.angles().len(), 6);
        assert_eq!(molecule.dihedrals().len(), 0);
    }
}

#[test]
fn constant_energy() {
    START.call_once(|| {Logger::stdout();});
    let path = Path::new(file!()).parent().unwrap().join("data")
                                     .join("md_methane.toml");
    let mut config = Input::new(path).unwrap().read().unwrap();


    let e_initial = config.system.total_energy();
    config.simulation.run(&mut config.system, config.nsteps);
    let e_final = config.system.total_energy();
    assert!(f64::abs((e_initial - e_final)/e_final) < 1e-2);
}