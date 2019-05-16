// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

//! Using a custom potential in simulations
use lumol::energy::{PairInteraction, PairPotential, Potential};
use lumol::{Particle, Molecule, System, Vector3D};
use lumol::units;

use lumol::sim::{MolecularDynamics, Simulation};

/// Let's define a new version of the Lennard-Jones potential, using the
/// alternative form:
///
///         A         B
///  V =  -----  -  -----
///       r^12       r^6
///
#[derive(Clone)]
struct LJ {
    a: f64,
    b: f64,
}

// All we need to do is to implement the Potential trait
impl Potential for LJ {
    // The energy function give the energy at distance `r`
    fn energy(&self, r: f64) -> f64 {
        self.a / r.powi(12) - self.b / r.powi(6)
    }

    // The force function give the norm of the force at distance `r`
    fn force(&self, r: f64) -> f64 {
        12.0 * self.a / r.powi(13) - 6.0 * self.b / r.powi(7)
    }
}

// We want to use our LJ potential as a pair potential.
impl PairPotential for LJ {
    // The long-range correction to the energy at the given cutoff
    fn tail_energy(&self, cutoff: f64) -> f64 {
        -(1.0 / 9.0 * self.a / cutoff.powi(9) - 1.0 / 3.0 * self.b / cutoff.powi(3))
    }

    // The long-range correction to the virial at the given cutoff
    fn tail_virial(&self, cutoff: f64) -> f64 {
        -(12.0 / 9.0 * self.a / cutoff.powi(9) - 6.0 / 3.0 * self.b / cutoff.powi(3))
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut system = System::new();
    system.add_molecule(Molecule::new(Particle::with_position("F", Vector3D::new(0.0, 0.0, 0.0))));
    system.add_molecule(Molecule::new(Particle::with_position("F", Vector3D::new(1.5, 0.0, 0.0))));

    // We can now use our new potential in the system
    let lj = Box::new(LJ {
        a: units::from(675.5, "kJ/mol/A^12")?,
        b: units::from(40.26, "kJ/mol/A^6")?,
    });
    system.set_pair_potential(("F", "F"), PairInteraction::new(lj, 10.0));

    let md = MolecularDynamics::new(units::from(1.0, "fs")?);
    let mut simulation = Simulation::new(Box::new(md));
    simulation.run(&mut system, 1000);

    Ok(())
}
