use lumol::sys::{System, Trajectory};
use lumol_input::InteractionsInput;
use std::path::Path;

use rand::{XorShiftRng, SeedableRng};

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

pub fn get_rng(seed: u32) -> XorShiftRng {
    XorShiftRng::from_seed([seed, 784, 71255487, 5824])
}
