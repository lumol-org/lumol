extern crate env_logger;
extern crate cymbalum;
use cymbalum::*;

fn main() {
    env_logger::init().unwrap();
    let mut universe = Universe::from_cell(UnitCell::cubic(17.0));

    // Create a cubic crystal of Argon by hand.
    for i in 0..5 {
        for j in 0..5 {
            for k in 0..5 {
                let mut part = Particle::new("Ar");
                part.set_position(Vector3D::new(
                        i as f64 * 3.4,
                        j as f64 * 3.4,
                        k as f64 * 3.4
                ));
                universe.add_particle(part);
            }
        }
    }
    universe.add_pair_interaction("Ar", "Ar", LennardJones{
                                                sigma: units::from(3.4, "A").unwrap(),
                                                epsilon: units::from(1.0, "kJ/mol").unwrap()});

    let mut velocities = BoltzmanVelocities::new(units::from(300.0, "K").unwrap());
    velocities.seed(129);
    velocities.init(&mut universe);

    let mut simulation = Simulation::new(MolecularDynamics::new(units::from(1.0, "fs").unwrap()));
    simulation.add_output(TrajectoryOutput::new("trajectory.xyz"));
    simulation.add_output(EnergyOutput::new("energy.dat"));

    simulation.run(&mut universe, 5000);

    println!("All done!")
}
