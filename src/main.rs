// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license
use std::fmt;

use backtrace::Backtrace;
use chrono::Duration;
use chrono::offset::Local;
use clap::{Command, Arg, ArgMatches};
use log::{info, error};

use lumol::input::Input;

fn parse_args() -> ArgMatches {
    Command::new("lumol").version(lumol::VERSION)
                     .about("An extensible molecular simulation engine")
                     .arg(
                        Arg::new("input.toml")
                            .required(true)
                            .help("Simulation input file")
                     )
                     .get_matches()
}

fn main() {
    std::panic::set_hook(Box::new(|info| {
        // Just in case no logger was created yet (very early panic), let's
        // create one!
        lumol_input::setup_default_logger();

        let playload = info.payload();
        let message = if let Some(message) = playload.downcast_ref::<&str>() {
            message
        } else if let Some(message) = playload.downcast_ref::<String>() {
            message
        } else {
            "<no message>"
        };

        error!("fatal error: {}", message);
        error!("{:?}", CleanedBacktrace::new());
    }));

    let args = parse_args();

    let input = args.get_one::<String>("input.toml").unwrap();
    let input = match Input::new(input) {
        Ok(input) => input,
        Err(err) => {
            panic!("invalid input file: {}", err);
        }
    };

    let mut config = match input.read() {
        Ok(config) => config,
        Err(err) => {
            panic!("bad input file: {}", err);
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

struct CleanedBacktrace {
    inner: Backtrace
}

impl CleanedBacktrace {
    fn new() -> CleanedBacktrace {
        CleanedBacktrace {
            inner: Backtrace::new()
        }
    }
}

impl fmt::Debug for CleanedBacktrace {
    fn fmt(&self, fmt: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(fmt, "stack backtrace:")?;

        let mut skiped_backtrace_generation = false;
        'outer: for frame in self.inner.frames() {
            for symbol in frame.symbols() {
                let name = if let Some(name) = symbol.name() {
                    // Use Display to get demangled name
                    format!("{}", name)
                } else {
                    "<unknown>".to_string()
                };

                if !skiped_backtrace_generation {
                    if name.starts_with("std::panicking::begin_panic") {
                        skiped_backtrace_generation = true;
                    }
                    continue 'outer;
                }

                if name.starts_with("std::rt::lang_start::") {
                    break 'outer;
                }

                write!(fmt, "\n - {}", name)?;
                if let (Some(file), Some(line)) = (symbol.filename(), symbol.lineno()) {
                    write!(fmt, "\n    at {}:{}", file.display(), line)?;
                }
            }
        }

        Ok(())
    }
}
