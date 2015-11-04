/*
 * Cymbalum, Molecular Simulation in Rust
 * Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
*/
//! Testing physical properties of a Lennard-Jones gaz of Helium using
//! Monte-Carlo simulation
extern crate cymbalum;
use self::cymbalum::*;

use std::sync::{Once, ONCE_INIT};
static START: Once = ONCE_INIT;

use std::path::Path;

fn get_universe() -> Universe {
    let data_dir = Path::new(file!()).parent().unwrap();
    let configuration = data_dir.join("data").join("helium.xyz");
    let mut universe = Universe::from_file(configuration.to_str().unwrap()).unwrap();
    universe.set_cell(UnitCell::cubic(10.0));

    let mut velocities = BoltzmanVelocities::new(units::from(300.0, "K").unwrap());
    velocities.init(&mut universe);

    universe.add_pair_interaction("He", "He",
        Box::new(LennardJones{
            sigma: units::from(2.0, "A").unwrap(),
            epsilon: units::from(0.2, "kJ/mol").unwrap()
        })
    );

    return universe;
}

#[test]
fn perfect_gaz() {
    START.call_once(|| {Logger::stdout();});
    let mut universe = get_universe();
    let mut mc = MonteCarlo::new(units::from(300.0, "K").unwrap());
    mc.add(Box::new(Translate::new(units::from(3.0, "A").unwrap())), 1.0);
    let mut simulation = Simulation::new(mc);

    // dilating the universe!
    for particle in universe.iter_mut() {
        particle.position = 10.0 * particle.position;
    }
    universe.set_cell(UnitCell::cubic(100.0));

    simulation.run(&mut universe, 5000);
    let pressure = universe.pressure();
    let volume = universe.volume();
    let temperature = universe.temperature();
    let n = universe.size() as f64;

    assert!(f64::abs(pressure * volume - n * constants::K_BOLTZMANN * temperature) < 1e-3);
}
