//! Molecular dynamics simulation of a crystal NaCl, reading universe and
//! potentials from files.
extern crate cymbalum;
use cymbalum::*;

fn main() {
    Logger::stdout();

    // Read the universe fromt the `data/nacl.xyz` file
    let mut universe = Universe::from_file("data/NaCl.xyz").unwrap();
    // Set the unit cell, as there is no unit cell data in XYZ files
    universe.set_cell(UnitCell::cubic(units::from(22.5608, "A").unwrap()));
    // Read the interactions from the `data/NaCl.yml` YAML file
    input::read_interactions(&mut universe, "data/NaCl.yml").unwrap();

    let mut velocities = BoltzmanVelocities::new(units::from(300.0, "K").unwrap());
    velocities.init(&mut universe);

    let mut simulation = Simulation::new(MolecularDynamics::new(units::from(1.0, "fs").unwrap()));
    simulation.run(&mut universe, 1000);
}
