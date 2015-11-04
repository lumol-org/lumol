//! Using a custom potential in simulations
extern crate cymbalum;
use cymbalum::*;

/// Let's define a new version of the LennardJones potential, using the
/// alternative form:
///
///         A         B
///  V =  -----  -  -----
///       r^12       r^6
///
struct LJ {
    a: f64,
    b: f64
}

/// All we need to do is to implement the PotentialFunction trait
impl PotentialFunction for LJ {
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
    let mut universe = Universe::new();

    universe.add_particle(Particle::new("F"));
    universe[0].position = Vector3D::new(0.0, 0.0, 0.0);
    universe.add_particle(Particle::new("F"));
    universe[1].position = Vector3D::new(1.5, 0.0, 0.0);

    // We can now use our new potential in the universe
    universe.add_pair_interaction("F", "F",
        Box::new(LJ{
            a: units::from(675.5, "kJ/mol/A^12").unwrap(),
            b: units::from(40.26, "kJ/mol/A^6").unwrap()
        }
    ));

    let mut simulation = Simulation::new(MolecularDynamics::new(units::from(1.0, "fs").unwrap()));
    simulation.run(&mut universe, 1000);
}
