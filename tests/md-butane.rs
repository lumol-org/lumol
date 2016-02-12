/*
 * Cymbalum, Molecular Simulation in Rust
 * Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
*/
//! Testing molecular dynamics of butane
extern crate cymbalum;
use self::cymbalum::*;

use std::path::Path;
use std::sync::{Once, ONCE_INIT};
static START: Once = ONCE_INIT;

fn setup_system() -> System {
    let data_dir = Path::new(file!()).parent().unwrap().join("data");
    let configuration = data_dir.join("butane.xyz");
    let mut system = System::from_file_auto_bonds(configuration.to_str().unwrap()).unwrap();
    system.set_cell(UnitCell::cubic(20.0));

    let interactions = data_dir.join("butane.yml");
    input::read_interactions(&mut system, interactions).unwrap();

    let mut velocities = BoltzmanVelocities::new(units::from(300.0, "K").unwrap());
    velocities.init(&mut system);

    return system;
}

#[test]
fn bonds_detection() {
    START.call_once(|| {Logger::stdout();});
    let system = setup_system();
    assert_eq!(system.molecules().len(), 50);

    for molecule in system.molecules() {
        assert_eq!(molecule.bonds().len(), 3);
        assert_eq!(molecule.angles().len(), 2);
        assert_eq!(molecule.dihedrals().len(), 1);
    }
}

#[test]
fn constant_energy() {
    START.call_once(|| {Logger::stdout();});
    let mut system = setup_system();

    let mut simulation = Simulation::new(
        MolecularDynamics::new(units::from(1.0, "fs").unwrap())
    );

    let e_initial = system.total_energy();
    simulation.run(&mut system, 1000);
    let e_final = system.total_energy();
    assert!(f64::abs((e_initial - e_final)/e_final) < 1e-3);
}
