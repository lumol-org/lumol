//! Example of a run using input files for the simulation and the system
//! This is the exact same simulation as the one in `binary.rs`
extern crate lumol;
extern crate lumol_input as input;

use lumol::Logger;
use input::Input;

fn main() {
    Logger::stdout();
    match Input::new("data/simulation.toml").and_then(|input| input.read()) {
        Err(error) => println!("Error in input: {}", error),
        Ok(mut config) => {
            config.simulation.run(&mut config.system, config.nsteps);
        }
    }
}
