extern crate lumol;
extern crate lumol_input;

use lumol_input::read_config;

use std::env;
use std::process::exit;

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.contains(&"-h".into()) || args.contains(&"--help".into()) || args.len() != 2 {
        return usage();
    }

    let input = &args[1];
    let mut config = match read_config(input) {
        Ok(config) => config,
        Err(err) => {
            println!("Error in input file: {}", err);
            exit(2);
        }
    };

    config.simulation.run(&mut config.system, config.nsteps);
}

fn usage() {
    let name = env::args().next().unwrap_or("lumol".into());
    println!("Usage: {} input.toml", name);
}
