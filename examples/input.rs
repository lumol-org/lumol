// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

//! Example of a run using input files for the simulation and the system
//! This is the exact same simulation as the one in `binary.rs`
fn main() {
    let input = lumol::input::Input::new("data/simulation.toml").unwrap();
    match input.read() {
        Err(error) => println!("Error in input: {}", error),
        Ok(mut config) => {
            config.simulation.run(&mut config.system, config.nsteps);
        }
    }
}
