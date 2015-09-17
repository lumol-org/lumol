/*
 * Geometry minization of a molecule of water
 */

extern crate env_logger;
extern crate cymbalum;
use cymbalum::*;

fn main() {
    env_logger::init().unwrap();
    let mut universe = Universe::new();

    let alpha = units::from(50.0, "deg").unwrap();

    universe.add_particle(Particle::new("O"));
    universe[0].set_position(Vector3D::new(0.0, 0.0, 0.0));
    universe.add_particle(Particle::new("H"));
    universe[1].set_position(Vector3D::new(1.2*f64::cos(alpha), 1.2*f64::sin(alpha), 0.0));
    universe.add_particle(Particle::new("H"));
    universe[2].set_position(Vector3D::new(1.2*f64::cos(-alpha), 1.2*f64::sin(-alpha), 0.0));

    {
        let topology = universe.topology_mut();
        topology.add_bond(0, 1);
        topology.add_bond(0, 2);
    }

    universe.add_pair_interaction("O", "H", NullPotential);
    universe.add_pair_interaction("O", "O", NullPotential);
    universe.add_pair_interaction("H", "H", NullPotential);

    universe.add_bond_interaction("O", "H", Harmonic{
        x0: units::from(1.1, "A").unwrap(),
        k: units::from(100.0, "kJ/mol/A^2").unwrap(),
    });
    universe.add_angle_interaction("H", "O", "H", Harmonic{
        x0: units::from(109.0, "deg").unwrap(),
        k: units::from(30.0, "kJ/mol/deg").unwrap(),
    });


    let mut simulation = Simulation::new(GradientDescent::new());
    simulation.add_output(EnergyOutput::new("energy.dat"));
    simulation.run(&mut universe, 500);
}
