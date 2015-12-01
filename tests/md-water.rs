/*
 * Cymbalum, Molecular Simulation in Rust
 * Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
*/
//! Testing physical properties of f-SPC water
extern crate cymbalum;
use self::cymbalum::*;

use std::path::Path;
use std::sync::{Once, ONCE_INIT};
static START: Once = ONCE_INIT;

fn setup(potential: &str, n: usize) -> (Simulation, Universe) {
    let mut universe = Universe::from_cell(UnitCell::cubic(28.0));

    let mut origins = Vec::new();
    let delta = 28.0 / n as f64;
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                origins.push(Vector3D::new(i as f64 * delta, j as f64 * delta, k as f64 * delta));
            }
        }
    }

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
    let potentials = data_dir.join("data").join(potential);
    input::read_interactions(&mut universe, potentials).unwrap();

    let mut velocities = BoltzmanVelocities::new(units::from(300.0, "K").unwrap());
    velocities.init(&mut universe);

    let simulation = Simulation::new(
        MolecularDynamics::new(units::from(1.0, "fs").unwrap())
    );
    return (simulation, universe);
}

#[test]
fn constant_energy_ewald() {
    START.call_once(|| {Logger::stdout();});
    let (mut simulation, mut universe) = setup("water-ewald.yml", 2);

    let e_initial = universe.total_energy();
    simulation.run(&mut universe, 1000);
    let e_final = universe.total_energy();
    assert!(f64::abs((e_initial - e_final)/e_final) < 3e-2);
}

#[test]
fn constant_energy_wolf() {
    START.call_once(|| {Logger::stdout();});
    let (mut simulation, mut universe) = setup("water-wolf.yml", 5);

    let e_initial = universe.total_energy();
    simulation.run(&mut universe, 1000);
    let e_final = universe.total_energy();
    assert!(f64::abs((e_initial - e_final)/e_final) < 3e-2);
}
