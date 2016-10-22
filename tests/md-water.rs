// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

//! Testing physical properties of f-SPC water
extern crate lumol;
extern crate lumol_input as input;

use lumol::Logger;
use lumol::types::Vector3D;
use lumol::sys::{System, UnitCell, Particle};
use lumol::sys::veloc::{BoltzmannVelocities, InitVelocities};
use lumol::sim::{Simulation, MolecularDynamics};
use lumol::units;

use input::InteractionsInput;

use std::path::Path;
use std::sync::{Once, ONCE_INIT};
static START: Once = ONCE_INIT;

fn setup(potential: &str, n: usize) -> (Simulation, System) {
    let mut system = System::from_cell(UnitCell::cubic(28.0));

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
        system.add_particle(Particle::new("O"));
        system.add_particle(Particle::new("H"));
        system.add_particle(Particle::new("H"));

        system[3*i + 0].position = origin.clone();
        system[3*i + 1].position = origin.clone() + h_1.clone();
        system[3*i + 2].position = origin.clone() + h_2.clone();

        system.add_bond(3*i, 3*i + 1);
        system.add_bond(3*i, 3*i + 2);
    }

    let data_dir = Path::new(file!()).parent().unwrap();
    let potentials = data_dir.join("data").join(potential);
    let input = InteractionsInput::new(potentials).unwrap();
    input.read(&mut system).unwrap();

    let mut velocities = BoltzmannVelocities::new(units::from(300.0, "K").unwrap());
    velocities.init(&mut system);

    let simulation = Simulation::new(Box::new(
        MolecularDynamics::new(units::from(1.0, "fs").unwrap())
    ));
    return (simulation, system);
}

#[test]
fn constant_energy_ewald() {
    START.call_once(|| {Logger::stdout();});
    let (mut simulation, mut system) = setup("water-ewald.toml", 2);

    let e_initial = system.total_energy();
    simulation.run(&mut system, 1000);
    let e_final = system.total_energy();

    // TODO: use a better thresold when updating Ewald to work with triclinic
    // cells.
    assert!(f64::abs((e_initial - e_final)/e_final) < 1e-1);
}

#[test]
fn constant_energy_wolf() {
    START.call_once(|| {Logger::stdout();});
    let (mut simulation, mut system) = setup("water-wolf.toml", 5);

    let e_initial = system.total_energy();
    simulation.run(&mut system, 1000);
    let e_final = system.total_energy();
    assert!(f64::abs((e_initial - e_final)/e_final) < 3e-2);
}
