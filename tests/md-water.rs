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

fn setup() -> (Simulation, Universe) {
    let mut universe = Universe::from_cell(UnitCell::cubic(7.0));

    let origins = vec![Vector3D::new(1.0, 1.0, 1.0),
                       Vector3D::new(5.0, 5.0, 1.0),
                       Vector3D::new(1.0, 5.0, 5.0),
                       Vector3D::new(5.0, 1.0, 5.0)];

    let h_1 = Vector3D::new(0.634859709040957, 0.8983057106778469, 0.0);
    let h_2 = Vector3D::new(-0.634859709040957, 0.8983057106778469, 0.0);

    for (i, origin) in origins.iter().enumerate() {
        universe.add_particle(Particle::new("O"));
        universe.add_particle(Particle::new("H"));
        universe.add_particle(Particle::new("H"));

        universe[3*i + 0].position = origin.clone();
        universe[3*i + 1].position = origin.clone() + h_1.clone();
        universe[3*i + 2].position = origin.clone() + h_2.clone();

        universe.add_bond(3*i, 3*i + 1);
        universe.add_bond(3*i, 3*i + 2);
    }

    let data_dir = Path::new(file!()).parent().unwrap();
    let potentials = data_dir.join("data").join("water.yml");
    input::read_interactions(&mut universe, potentials).unwrap();

    let mut velocities = BoltzmanVelocities::new(units::from(300.0, "K").unwrap());
    velocities.init(&mut universe);

    let simulation = Simulation::new(
        MolecularDynamics::new(units::from(1.0, "fs").unwrap())
    );
    return (simulation, universe);
}

#[test]
fn constant_energy() {
    Logger::stdout();
    let (mut simulation, mut universe) = setup();

    let E_initial = universe.total_energy();
    simulation.run(&mut universe, 1000);
    let E_final = universe.total_energy();
    assert!(f64::abs(E_initial - E_final)/E_final < 1.1e-2);
}
