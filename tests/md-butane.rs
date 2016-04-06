// Cymbalum, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

//! Testing molecular dynamics of butane
extern crate cymbalum;
use cymbalum::*;

use std::path::Path;
use std::sync::{Once, ONCE_INIT};
static START: Once = ONCE_INIT;

fn setup_system() -> System {
    let data_dir = Path::new(file!()).parent().unwrap().join("data");
    let configuration = data_dir.join("butane.xyz");
    let mut system = input::Trajectory::open(configuration)
                                        .and_then(|mut traj| traj.read_guess_bonds())
                                        .unwrap();
    system.set_cell(UnitCell::cubic(20.0));

    let interactions = data_dir.join("butane.toml");
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
