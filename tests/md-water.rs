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

extern crate env_logger;

extern crate cymbalum;
use self::cymbalum::*;

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

        universe[3*i + 0].set_position(origin.clone());
        universe[3*i + 1].set_position(origin.clone() + h_1.clone());
        universe[3*i + 2].set_position(origin.clone() + h_2.clone());

        let topology = universe.topology_mut();
        topology.add_bond(3*i, 3*i + 1);
        topology.add_bond(3*i, 3*i + 2);
    }

    universe.add_pair_interaction("O", "O",
        LennardJones{
            sigma: units::from(3.2, "A").unwrap(),
            epsilon: units::from(0.2583, "kcal/mol").unwrap()
        }
    );

    universe.add_pair_interaction("O", "H", NullPotential);
    universe.add_pair_interaction("H", "H", NullPotential);

    universe.add_bond_interaction("O", "H",
        Harmonic{
            x0: units::from(1.1, "A").unwrap(),
            k: units::from(390.0, "kcal/mol/A^2").unwrap()
        }
    );

    universe.add_angle_interaction("H", "O", "H",
        Harmonic{
            x0: units::from(109.5, "deg").unwrap(),
            k: units::from(70.0, "kcal/mol/rad^2").unwrap()
        }
    );

    let mut velocities = BoltzmanVelocities::new(units::from(300.0, "K").unwrap());
    velocities.init(&mut universe);

    let simulation = Simulation::new(
        MolecularDynamics::new(units::from(1.0, "fs").unwrap())
    );
    return (simulation, universe);
}

#[test]
fn constant_energy() {
    env_logger::init().unwrap();
    let (mut simulation, mut universe) = setup();

    let E_initial = universe.total_energy();
    simulation.run(&mut universe, 1000);
    let E_final = universe.total_energy();
    assert!(f64::abs(E_initial - E_final) < 1e-2);
}
