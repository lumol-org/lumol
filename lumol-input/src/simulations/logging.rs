// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license
use toml::value::Table;

use log::{self, Record, info};

use log4rs::append::Append;
use log4rs::append::console;
use log4rs::append::console::ConsoleAppender;
use log4rs::append::file::FileAppender;
use log4rs::config::{Appender, Config, Root};
use log4rs::encode::{Color, Encode, Style, Write};
use log4rs::filter::threshold::ThresholdFilter;

use crate::Input;
use crate::error::Error;
use crate::extract;

/// An encoder to configure the output style of logging messages
#[derive(Debug)]
struct LogEncoder;

impl Encode for LogEncoder {
    fn encode(&self, out: &mut dyn Write, record: &Record<'_>) -> anyhow::Result<()> {
        match record.level() {
            log::Level::Trace => write!(out, "[trace] ")?,
            log::Level::Debug => write!(out, "[debug] ")?,
            log::Level::Info => {}
            log::Level::Warn => {
                out.set_style(Style::new().text(Color::Red))?;
                write!(out, "[warning] ")?;
                out.set_style(&Style::new())?;
            }
            log::Level::Error => {
                out.set_style(Style::new().text(Color::Red).intense(true))?;
                write!(out, "[error] ")?;
            }
        }

        write!(out, "{}", record.args())?;

        match record.level() {
            log::Level::Trace | log::Level::Debug => {
                if let Some(file) = record.file() {
                    match record.line() {
                        Some(line) => write!(out, " [from {file}:{line}]")?,
                        None => write!(out, " [from {file}]")?,
                    }
                }
            }
            log::Level::Error => out.set_style(&Style::new())?,
            _ => {}
        }

        writeln!(out)?;
        Ok(())
    }
}


/// Ensure that a logger will be initialized, even in case of error in the log
/// section of the input file
struct EnsureLogger;
impl Drop for EnsureLogger {
    fn drop(&mut self) {
        // Setup default logger unconditionally, it won't do anything if another
        // logger is already initialized.
        setup_default_logger();
    }
}

impl Input {
    /// Setup the logger from the input file or to stdout as a default.
    pub(crate) fn setup_logging(&self) -> Result<(), Error> {
        let _guard = EnsureLogger;
        if let Some(loggers) = self.config.get("log") {
            let loggers = loggers.as_table().ok_or(
                Error::from("'log' section must be a table")
            )?;

            if loggers.get("target").is_some() {
                if loggers.get("targets").is_some() {
                    return Err(
                        Error::from("can not have both 'target' and 'targets' in the log section"),
                    );
                }

                let appender = read_appender(loggers, "main")?;
                let config = Config::builder().appender(appender)
                                              .build(
                    Root::builder().appender("main").build(log::LevelFilter::Trace),
                )
                                              .expect("Error in logging initialization");
                let _ = log4rs::init_config(config);
            } else if let Some(targets) = loggers.get("targets") {
                let targets = targets.as_array().ok_or(
                    Error::from("'targets' must be an array in 'log' section")
                )?;

                let mut appenders = Vec::new();
                for (i, target) in targets.iter().enumerate() {
                    let target = target.as_table().ok_or(
                        Error::from("'targets' must be an array of tables in 'log' section")
                    )?;
                    let appender = read_appender(target, &i.to_string())?;
                    appenders.push(appender);
                }

                let mut root = Root::builder();
                let mut config = Config::builder();
                for appender in appenders {
                    root = root.appender(appender.name());
                    config = config.appender(appender);
                }
                let config = config.build(root.build(log::LevelFilter::Trace))
                                   .expect("Error in logging initialization");
                let _ = log4rs::init_config(config);
            } else {
                return Err(Error::from("missing 'target' or 'targets' in log section"));
            }
        } else {
            setup_default_logger();
            info!("Printing logs to the console as default.");
        }
        Ok(())
    }
}

/// Setup a default logger to be able to print error messages
pub fn setup_default_logger() {
    // We just log everything to stdout
    let stdout = ConsoleAppender::builder()
        .target(console::Target::Stdout)
        .encoder(Box::new(LogEncoder))
        .build();

    let appender = Appender::builder()
        .filter(Box::new(ThresholdFilter::new(log::LevelFilter::Info)))
        .build("main", Box::new(stdout));

    let config = Config::builder()
        .appender(appender)
        .build(Root::builder().appender("main").build(log::LevelFilter::Info))
        .expect("Error in logging initialization");

    // We ignore the result of this call, because it can only fail if a
    // logger has already been initialized.
    let _ = log4rs::init_config(config);
}

fn read_appender(config: &Table, name: &str) -> Result<Appender, Error> {
    let allowed_keys = ["target", "targets", "level", "append"];
    for key in config.keys() {
        if !allowed_keys.contains(&&**key) {
            return Err(Error::from(format!("unknown '{key}' key in log section")));
        }
    }

    let level = config.get("level")
        .map_or(Some("info"), |level| level.as_str())
        .ok_or(Error::from("'level' must be a string in log target"))?;

    let level = match level {
        "trace" => log::LevelFilter::Trace,
        "debug" => log::LevelFilter::Debug,
        "info" => log::LevelFilter::Info,
        "warning" => log::LevelFilter::Warn,
        "error" => log::LevelFilter::Error,
        other => return Err(Error::from(format!("unknown logging level '{other}'"))),
    };

    let target = extract::str("target", config, "log target")?;
    let appender: Box<dyn Append> = match target {
        "<stdout>" => {
            Box::new(
                ConsoleAppender::builder()
                    .target(console::Target::Stdout)
                    .encoder(Box::new(LogEncoder))
                    .build(),
            )
        }
        "<stderr>" => {
            Box::new(
                ConsoleAppender::builder()
                    .target(console::Target::Stderr)
                    .encoder(Box::new(LogEncoder))
                    .build(),
            )
        }
        "" => return Err(Error::from("'target' can not be an empty string in log target")),
        filename => {
            let append_mode = config.get("append")
                .map_or(Some(false), |x| x.as_bool())
                .ok_or(Error::from("'append' must be a boolean in log file target"))?;

            let appender = FileAppender::builder()
                .append(append_mode)
                .encoder(Box::new(LogEncoder));

            let appender = try_io!(appender.build(filename), filename.into());
            Box::new(appender)
        }
    };

    let appender = Appender::builder()
        .filter(Box::new(ThresholdFilter::new(level)))
        .build(name, appender);

    Ok(appender)
}
