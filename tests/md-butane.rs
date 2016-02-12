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

fn setup_universe() -> Universe {
    let data_dir = Path::new(file!()).parent().unwrap().join("data");
    let configuration = data_dir.join("butane.xyz");
    let mut universe = Universe::from_file_auto_bonds(configuration.to_str().unwrap()).unwrap();
    universe.set_cell(UnitCell::cubic(20.0));

    let interactions = data_dir.join("butane.yml");
    input::read_interactions(&mut universe, interactions).unwrap();

    let mut velocities = BoltzmanVelocities::new(units::from(300.0, "K").unwrap());
    velocities.init(&mut universe);

    return universe;
}

#[test]
fn bonds_detection() {
    START.call_once(|| {Logger::stdout();});
    let universe = setup_universe();
    assert_eq!(universe.molecules().len(), 50);

    for molecule in universe.molecules() {
        assert_eq!(molecule.bonds().len(), 3);
        assert_eq!(molecule.angles().len(), 2);
        assert_eq!(molecule.dihedrals().len(), 1);
    }
}

#[test]
fn constant_energy() {
    START.call_once(|| {Logger::stdout();});
    let mut universe = setup_universe();

    let mut simulation = Simulation::new(
        MolecularDynamics::new(units::from(1.0, "fs").unwrap())
    );

    let e_initial = universe.total_energy();
    simulation.run(&mut universe, 1000);
    let e_final = universe.total_energy();
    assert!(f64::abs((e_initial - e_final)/e_final) < 1e-3);
}
