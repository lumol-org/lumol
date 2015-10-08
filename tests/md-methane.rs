/*
 * Cymbalum, Molecular Simulation in Rust
 * Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
*/
#![allow(non_snake_case)]
//! Testing physical properties of 4 molecules of water in a box
extern crate cymbalum;
use self::cymbalum::*;

use std::path::Path;
use std::sync::{Once, ONCE_INIT};
static START: Once = ONCE_INIT;

fn setup_universe() -> Universe {
    let data_dir = Path::new(file!()).parent().unwrap();
    let configuration = data_dir.join("data").join("methane.xyz");
    let mut universe = Universe::from_file_auto_bonds(configuration.to_str().unwrap()).unwrap();
    universe.set_cell(UnitCell::cubic(20.0));

    let potentials = data_dir.join("data").join("methane.yml");
    input::read_potentials(&mut universe, potentials).unwrap();

    let mut velocities = BoltzmanVelocities::new(units::from(300.0, "K").unwrap());
    velocities.init(&mut universe);

    return universe;
}

#[test]
fn bonds_detection() {
    START.call_once(|| {Logger::stdout();});
    let universe = setup_universe();
    assert_eq!(universe.topology().molecules().len(), 150);
    assert_eq!(universe.bonds().len(), 600);
    assert_eq!(universe.angles().len(), 900);
    assert_eq!(universe.dihedrals().len(), 0);
}

#[test]
fn constant_energy() {
    START.call_once(|| {Logger::stdout();});
    let mut universe = setup_universe();

    let mut simulation = Simulation::new(
        MolecularDynamics::new(units::from(1.0, "fs").unwrap())
    );

    let E_initial = universe.total_energy();
    simulation.run(&mut universe, 500);
    let E_final = universe.total_energy();
    assert!(f64::abs(E_initial - E_final)/E_final < 1e-2);
}
