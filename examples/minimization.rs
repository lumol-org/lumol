//! Geometry minization of a molecule of water
extern crate cymbalum;
use cymbalum::*;

fn main() {
    Logger::stdout();
    let mut universe = Universe::new();

    let alpha = units::from(50.0, "deg").unwrap();

    universe.add_particle(Particle::new("O"));
    universe[0].position = Vector3D::new(0.0, 0.0, 0.0);
    universe.add_particle(Particle::new("H"));
    universe[1].position = Vector3D::new(1.2*f64::cos(alpha), 1.2*f64::sin(alpha), 0.0);
    universe.add_particle(Particle::new("H"));
    universe[2].position = Vector3D::new(1.2*f64::cos(-alpha), 1.2*f64::sin(-alpha), 0.0);

    universe.add_bond(0, 1);
    universe.add_bond(0, 2);

    universe.add_pair_interaction("O", "H", Box::new(NullPotential));
    universe.add_pair_interaction("O", "O", Box::new(NullPotential));
    universe.add_pair_interaction("H", "H", Box::new(NullPotential));

    universe.add_bond_interaction("O", "H", Box::new(Harmonic{
        x0: units::from(1.1, "A").unwrap(),
        k: units::from(100.0, "kJ/mol/A^2").unwrap(),
    }));
    universe.add_angle_interaction("H", "O", "H", Box::new(Harmonic{
        x0: units::from(109.0, "deg").unwrap(),
        k: units::from(30.0, "kJ/mol/deg").unwrap(),
    }));


    let mut simulation = Simulation::new(SteepestDescent::new());
    simulation.add_output(EnergyOutput::new("energy.dat").unwrap());
    simulation.run(&mut universe, 500);
}
