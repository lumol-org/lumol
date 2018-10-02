// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license
extern crate lumol;
extern crate lumol_input;

extern crate chrono;
extern crate clap;
#[macro_use]
extern crate log;

use chrono::Duration;
use chrono::offset::Local;
use clap::{App, ArgMatches};
use lumol_input::Input;

fn parse_args<'a>() -> ArgMatches<'a> {
    App::new("lumol").version(lumol::VERSION)
                     .about("An extensible molecular simulation engine")
                     .args_from_usage("<input.toml>      'Simulation input file'")
                     .get_matches()
}

fn main() {
    let args = parse_args();

    let input = args.value_of("input.toml").unwrap();
    let input = match Input::new(input) {
        Ok(input) => input,
        Err(err) => {
            lumol_input::setup_default_logger();
            error!("invalid input file: {}", err);
            std::process::exit(2)
        }
    };

    let mut config = match input.read() {
        Ok(config) => config,
        Err(err) => {
            error!("bad input file: {}", err);
            std::process::exit(2)
        }
    };

    info!("Running lumol version {}", lumol::VERSION);

    let start = Local::now();
    info!(
        "Simulation started the {} at {}",
        start.format("%Y-%m-%d"),
        start.format("%H:%M:%S")
    );
    info!(" "); // Skip a line

    config.simulation.run(&mut config.system, config.nsteps);

    let end = Local::now();
    info!(
        "\nSimulation ended the {} at {}",
        end.format("%Y-%m-%d"),
        end.format("%H:%M:%S")
    );
    let elapsed = end.signed_duration_since(start);
    info!("Simulation ran for {}", format_elapsed(elapsed));
}

fn format_elapsed(elapsed: Duration) -> String {
    if elapsed.num_weeks() > 0 {
        let h = elapsed.num_hours() % 24;
        let d = elapsed.num_days() % 7;
        let w = elapsed.num_weeks();
        format!("{} weeks {} days {}h", w, d, h)
    } else if elapsed.num_days() > 0 {
        let m = elapsed.num_minutes() % 60;
        let h = elapsed.num_hours() % 24;
        let d = elapsed.num_days();
        format!("{} days {}h {}min", d, h, m)
    } else if elapsed.num_hours() > 0 {
        let s = elapsed.num_seconds() % 60;
        let m = elapsed.num_minutes() % 60;
        let h = elapsed.num_hours();
        format!("{}h {}min {}s", h, m, s)
    } else {
        let s = elapsed.num_seconds() % 60;
        let m = elapsed.num_minutes();
        format!("{}min {}s", m, s)
    }
}
