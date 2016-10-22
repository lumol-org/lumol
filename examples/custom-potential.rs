//! Using a custom potential in simulations
extern crate lumol;

use lumol::Logger;
use lumol::types::{Vector3D, Zero};
use lumol::sys::{System, Particle};
use lumol::energy::{Potential, PairPotential, PairInteraction};
use lumol::sim::{Simulation, MolecularDynamics};
use lumol::units;

/// Let's define a new version of the LennardJones potential, using the
/// alternative form:
///
///         A         B
///  V =  -----  -  -----
///       r^12       r^6
///
#[derive(Clone)]
struct LJ {
    a: f64,
    b: f64
}

/// All we need to do is to implement the Potential trait
impl Potential for LJ {
    /// The energy function give the energy at distance `r`
    fn energy(&self, r: f64) -> f64 {
        self.a / r.powi(12) - self.b / r.powi(6)
    }

    /// The force function give the norm of the force at distance `r`
    fn force(&self, r: f64) -> f64 {
        12.0 * self.a / r.powi(13) - 6.0 * self.b / r.powi(7)
    }
}

// We want to use our LJ potential as a pair potential.
impl PairPotential for LJ {}

fn main() {
    Logger::stdout();
    let mut system = System::new();

    system.add_particle(Particle::new("F"));
    system[0].position = Vector3D::zero();
    system.add_particle(Particle::new("F"));
    system[1].position = Vector3D::new(1.5, 0.0, 0.0);

    // We can now use our new potential in the system
    let lj = Box::new(LJ {
        a: units::from(675.5, "kJ/mol/A^12").unwrap(),
        b: units::from(40.26, "kJ/mol/A^6").unwrap()
    });
    system.interactions_mut().add_pair("F", "F", PairInteraction::new(lj, 10.0));

    let mut simulation = Simulation::new(Box::new(
        MolecularDynamics::new(units::from(1.0, "fs").unwrap())
    ));
    simulation.run(&mut system, 1000);
}
