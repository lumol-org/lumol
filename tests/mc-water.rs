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

fn get_universe(potential: &str) -> Universe {
    let data_dir = Path::new(file!()).parent().unwrap().join("data");
    let configuration = data_dir.join("water.xyz");
    let mut universe = Universe::from_file(configuration.to_str().unwrap()).unwrap();
    universe.set_cell(UnitCell::cubic(18.0));

    let potentials = data_dir.join(potential);
    input::read_interactions(&mut universe, potentials).unwrap();
    return universe;
}

// This test only run a Monte-Carlo simulation of water, but do not test
// anything for now. It should test the g(r) function someday.
#[test]
fn run() {
    START.call_once(|| {Logger::stdout();});
    let mut universe = get_universe("water-wolf.yml");

    let mut mc = MonteCarlo::new(units::from(300.0, "K").unwrap());
    mc.add(Box::new(Translate::new(units::from(3.0, "A").unwrap())), 1.0);
    mc.add(Box::new(Rotate::new(units::from(20.0, "deg").unwrap())), 1.0);
    let mut simulation = Simulation::new(mc);

    // TODO: increase this value when the cache is implemented
    simulation.run(&mut universe, 500);
}
