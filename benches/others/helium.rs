extern crate cymbalum;
use cymbalum::*;
use std::path::Path;

#[test]
fn main() {
    Logger::stdout();
    let data_dir = Path::new(file!()).parent().unwrap();
    let configuration = data_dir.join("data").join("helium.xyz");
    let mut system = System::from_file(configuration.to_str().unwrap()).unwrap();
    system.set_cell(UnitCell::cubic(10.0));

    system.add_pair_interaction("He", "He",
        LennardJones{
            sigma: units::from(2.0, "A").unwrap(),
            epsilon: units::from(0.3, "kJ/mol").unwrap()
        }
    );

    let mut velocities = BoltzmannVelocities::new(units::from(300.0, "K").unwrap());
    velocities.seed(42);
    velocities.init(&mut system);

    let mut simulation = Simulation::new(
        MolecularDynamics::new(units::from(1.0, "fs").unwrap())
    );
    simulation.run(&mut system, 5000);
}
