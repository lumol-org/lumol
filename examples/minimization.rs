// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

//! Geometry minization of a molecule of water
extern crate lumol;

use lumol::types::Vector3D;
use lumol::units;

use lumol::energy::Harmonic;
use lumol::sys::{Particle, Molecule, System};

use lumol::sim::{Minimization, Simulation};
use lumol::sim::min::SteepestDescent;

fn main() -> Result<(), Box<std::error::Error>> {
    let mut system = System::new();

    let alpha = units::from(50.0, "deg")?;
    let a_cos = 1.2 * f64::cos(alpha);
    let a_sin = 1.2 * f64::sin(alpha);
    let mut molecule = Molecule::new(Particle::with_position("O", Vector3D::new(0.0, 0.0, 0.0)));
    molecule.add_particle_bonded_to(0, Particle::with_position("H", Vector3D::new(a_cos, a_sin, 0.0)));
    molecule.add_particle_bonded_to(0, Particle::with_position("H", Vector3D::new(a_cos, -a_sin, 0.0)));
    system.add_molecule(molecule);

    system.add_bond_potential(
        ("O", "H"),
        Box::new(Harmonic {
            x0: units::from(1.1, "A")?,
            k: units::from(100.0, "kJ/mol/A^2")?,
        }),
    );
    system.add_angle_potential(
        ("H", "O", "H"),
        Box::new(Harmonic {
            x0: units::from(109.0, "deg")?,
            k: units::from(30.0, "kJ/mol/deg")?,
        }),
    );

    let minimization = Minimization::new(Box::new(SteepestDescent::new()));
    let mut simulation = Simulation::new(Box::new(minimization));
    simulation.run(&mut system, 500);

    return Ok(());
}
