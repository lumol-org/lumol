// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 Lumol's contributors — BSD license

//! Geometry minization of a molecule of water
extern crate lumol;

use lumol::sys::{System, Particle};
use lumol::types::{Vector3D, Zero};
use lumol::energy::{PairInteraction, NullPotential, Harmonic};
use lumol::sim::{Simulation, Minimization};
use lumol::sim::min::SteepestDescent;
use lumol::out::EnergyOutput;
use lumol::units;

fn main() {
    let mut system = System::new();

    let alpha = units::from(50.0, "deg").unwrap();

    system.add_particle(Particle::new("O"));
    system[0].position = Vector3D::zero();
    system.add_particle(Particle::new("H"));
    system[1].position = Vector3D::new(1.2*f64::cos(alpha), 1.2*f64::sin(alpha), 0.0);
    system.add_particle(Particle::new("H"));
    system[2].position = Vector3D::new(1.2*f64::cos(-alpha), 1.2*f64::sin(-alpha), 0.0);

    system.add_bond(0, 1);
    system.add_bond(0, 2);

    {
        let interactions = system.interactions_mut();
        let null_interaction = PairInteraction::new(Box::new(NullPotential), 10.0);
        interactions.add_pair("O", "H", null_interaction.clone());
        interactions.add_pair("O", "O", null_interaction.clone());
        interactions.add_pair("H", "H", null_interaction.clone());

        interactions.add_bond("O", "H", Box::new(Harmonic{
            x0: units::from(1.1, "A").unwrap(),
            k: units::from(100.0, "kJ/mol/A^2").unwrap(),
        }));
        interactions.add_angle("H", "O", "H", Box::new(Harmonic{
            x0: units::from(109.0, "deg").unwrap(),
            k: units::from(30.0, "kJ/mol/deg").unwrap(),
        }));
    }


    let mut simulation = Simulation::new(
        Box::new(Minimization::new(
            Box::new(SteepestDescent::new()))
        )
    );
    simulation.add_output(Box::new(EnergyOutput::new("energy.dat").unwrap()));
    simulation.run(&mut system, 500);
}
