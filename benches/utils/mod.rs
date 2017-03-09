use lumol::sys::{System, Trajectory};
use lumol_input::InteractionsInput;
use std::path::Path;

pub fn get_system(name: &str) -> System {
    let data = Path::new(file!()).parent().unwrap().join("..").join("data");

    let mut system = Trajectory::open(data.join(String::from(name) + ".pdb"))
                                .and_then(|mut trajectory| trajectory.read())
                                .unwrap();

    InteractionsInput::new(data.join(String::from(name) + ".toml"))
                      .and_then(|input| input.read(&mut system))
                      .unwrap();

    system
}
