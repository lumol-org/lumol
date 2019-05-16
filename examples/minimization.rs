// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

//! Geometry minization of a molecule of water
use lumol::energy::Harmonic;
use lumol::{Particle, Molecule, System, Vector3D};
use lumol::units;

use lumol::sim::{Minimization, Simulation};
use lumol::sim::min::{SteepestDescent, Tolerance};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut system = System::new();

    let alpha = units::from(50.0, "deg")?;
    let a_cos = 1.2 * f64::cos(alpha);
    let a_sin = 1.2 * f64::sin(alpha);
    let mut molecule = Molecule::new(Particle::with_position("O", Vector3D::new(0.0, 0.0, 0.0)));
    molecule.add_particle_bonded_to(0, Particle::with_position("H", Vector3D::new(a_cos, a_sin, 0.0)));
    molecule.add_particle_bonded_to(0, Particle::with_position("H", Vector3D::new(a_cos, -a_sin, 0.0)));
    system.add_molecule(molecule);

    system.set_bond_potential(
        ("O", "H"),
        Box::new(Harmonic {
            x0: units::from(1.1, "A")?,
            k: units::from(100.0, "kJ/mol/A^2")?,
        }),
    );
    system.set_angle_potential(
        ("H", "O", "H"),
        Box::new(Harmonic {
            x0: units::from(109.0, "deg")?,
            k: units::from(30.0, "kJ/mol/deg")?,
        }),
    );

    let tolerance = Tolerance {
        energy: units::from(1e-5, "kJ/mol")?,
        force2: units::from(1e-5, "(kJ/mol/A)^2")?,
    };
    let minimization = Minimization::new(Box::new(SteepestDescent::new()), tolerance);
    let mut simulation = Simulation::new(Box::new(minimization));
    simulation.run(&mut system, 500);

    Ok(())
}
