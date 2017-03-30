// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 Lumol's contributors — BSD license

//! Example of a run using input files for the simulation and the system
//! This is the exact same simulation as the one in `binary.rs`
extern crate lumol_input as input;

fn main() {
    match input::Input::new("data/simulation.toml").and_then(|input| input.read()) {
        Err(error) => println!("Error in input: {}", error),
        Ok(mut config) => {
            config.simulation.run(&mut config.system, config.nsteps);
        }
    }
}
