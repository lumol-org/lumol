use lumol::sys::{System, Trajectory};

use std::path::Path;

pub fn get_system(name: &str) -> System {
    let configuration = Path::new(file!())
                             .parent().unwrap()
                             .join("..")
                             .join("data")
                             .join(String::from(name) + ".pdb");
    let mut system = Trajectory::open(configuration)
                                .and_then(|mut trajectory| trajectory.read())
                                .unwrap();

    for particle in &mut system {
        let charge = match particle.name() {
            "Na" => 1.0,
            "Cl" => -1.0,
            "O" => -0.82,
            "H" => 0.41,
            other => panic!("Missing charge value for {}", other)
        };
        particle.charge = charge;
    }

    system
}
