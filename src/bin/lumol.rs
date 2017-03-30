// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 Lumol's contributors — BSD license

extern crate lumol;
extern crate lumol_input;

#[macro_use]
extern crate log;
extern crate clap;

use lumol_input::Input;
use clap::{App, ArgMatches};

use std::process::exit;

fn parse_args<'a>() -> ArgMatches<'a> {
    App::new("lumol")
        .version(env!("CARGO_PKG_VERSION"))
        .about("An extensible molecular simulation engine")
        .args_from_usage("<input.toml>      'Simulation input file'")
        .get_matches()
}

fn main() {
    let args = parse_args();
    let input = args.value_of("input.toml").unwrap();
    let mut config = match Input::new(input).and_then(|input| input.read()) {
        Ok(config) => config,
        Err(err) => {
            error!("bad input file: {}", err);
            exit(2);
        }
    };

    config.simulation.run(&mut config.system, config.nsteps);
}
