// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

extern crate lumol;
use lumol::*;

use std::path::Path;
use std::sync::{Once, ONCE_INIT};
static START: Once = ONCE_INIT;

fn get_system(potential: &str) -> System {
    let data_dir = Path::new(file!()).parent().unwrap().join("data");
    let configuration = data_dir.join("water.xyz");
    let mut system = input::Trajectory::open(configuration)
                                        .and_then(|mut traj| traj.read_guess_bonds())
                                        .unwrap();
    system.set_cell(UnitCell::cubic(18.0));

    let potentials = data_dir.join(potential);
    input::read_interactions(&mut system, potentials).unwrap();
    return system;
}

// This test only run a Monte-Carlo simulation of water, but do not test
// anything for now. It should test the g(r) function someday.
#[test]
fn wolf() {
    START.call_once(|| {Logger::stdout();});
    let mut system = get_system("water-wolf.toml");

    let mut mc = MonteCarlo::new(units::from(300.0, "K").unwrap());
    mc.add(Box::new(Translate::new(units::from(3.0, "A").unwrap())), 1.0);
    mc.add(Box::new(Rotate::new(units::from(20.0, "deg").unwrap())), 1.0);
    let mut simulation = Simulation::new(Box::new(mc));

    simulation.run(&mut system, 5000);
}

#[test]
fn ewald() {
    START.call_once(|| {Logger::stdout();});
    let mut system = get_system("water-ewald.toml");

    let mut mc = MonteCarlo::new(units::from(300.0, "K").unwrap());
    mc.add(Box::new(Translate::new(units::from(3.0, "A").unwrap())), 1.0);
    mc.add(Box::new(Rotate::new(units::from(20.0, "deg").unwrap())), 1.0);
    let mut simulation = Simulation::new(Box::new(mc));

    simulation.run(&mut system, 1000);
}
