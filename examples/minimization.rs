// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

//! Geometry minization of a molecule of water
extern crate lumol;

use lumol::energy::{Harmonic, NullPotential, PairInteraction};
use lumol::out::EnergyOutput;
use lumol::sim::{Minimization, Simulation};
use lumol::sim::min::SteepestDescent;
use lumol::sys::{Particle, System};
use lumol::types::Vector3D;
use lumol::units;

fn main() {
    let mut system = System::new();

    let alpha = units::from(50.0, "deg").unwrap();
    let a_cos = 1.2 * f64::cos(alpha);
    let a_sin = 1.2 * f64::sin(alpha);
    system.add_particle(Particle::with_position("O", Vector3D::new(0.0, 0.0, 0.0)));
    system.add_particle(Particle::with_position("H", Vector3D::new(a_cos, a_sin, 0.0)));
    system.add_particle(Particle::with_position("H", Vector3D::new(a_cos, -a_sin, 0.0)));

    system.add_bond(0, 1);
    system.add_bond(0, 2);

    let null_interaction = PairInteraction::new(Box::new(NullPotential), 10.0);
    system.add_pair_potential(("O", "H"), null_interaction.clone());
    system.add_pair_potential(("O", "O"), null_interaction.clone());
    system.add_pair_potential(("H", "H"), null_interaction.clone());

    system.add_bond_potential(
        ("O", "H"),
        Box::new(Harmonic {
            x0: units::from(1.1, "A").unwrap(),
            k: units::from(100.0, "kJ/mol/A^2").unwrap(),
        }),
    );
    system.add_angle_potential(
        ("H", "O", "H"),
        Box::new(Harmonic {
            x0: units::from(109.0, "deg").unwrap(),
            k: units::from(30.0, "kJ/mol/deg").unwrap(),
        }),
    );

    let minimization = Minimization::new(Box::new(SteepestDescent::new()));
    let mut simulation = Simulation::new(Box::new(minimization));
    simulation.add_output(Box::new(EnergyOutput::new("energy.dat").unwrap()));
    simulation.run(&mut system, 500);
}
