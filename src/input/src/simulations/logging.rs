// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 Lumol's contributors — BSD license
use error::Error;
use extract;
use Input;

use toml::value::Table;

use log::{LogLevel, LogLevelFilter};
use log::LogRecord;

use log4rs;
use log4rs::encode::{Encode, Write, Color, Style};
use log4rs::append::Append;
use log4rs::append::console;
use log4rs::append::console::ConsoleAppender;
use log4rs::append::file::FileAppender;
use log4rs::config::{Config, Appender, Root};
use log4rs::filter::threshold::ThresholdFilter;

type LogError = Box<::std::error::Error + Sync + Send>;

/// An encoder to configure the output style of logging messages
#[derive(Debug)]
struct LogEncoder;

impl Encode for LogEncoder {
    fn encode(&self, out: &mut Write, record: &LogRecord) -> Result<(), LogError> {
        match record.level() {
            LogLevel::Trace => try!(write!(out, "[trace] ")),
            LogLevel::Debug => try!(write!(out, "[debug] ")),
            LogLevel::Info => {},
            LogLevel::Warn => {
                try!(out.set_style(&Style::new().text(Color::Red)));
                try!(write!(out, "[warning] "));
                try!(out.set_style(&Style::new()));
            }
            LogLevel::Error => {
                try!(out.set_style(&Style::new().text(Color::Red).intense(true)));
                try!(write!(out, "[error] "));
            },
        }

        try!(write!(out, "{}", record.args()));

        match record.level() {
            LogLevel::Trace | LogLevel::Debug =>{
                let location = record.location();
                try!(write!(out, " [from {}:{}]", location.file(), location.line()))
            }
            LogLevel::Error => {
                try!(out.set_style(&Style::new()));
            },
            _ => {}
        }

        try!(write!(out, "\n"));
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
    // TODO: use restricted privacy here
    #[doc(hidden)]
    pub fn setup_logging(&self) -> Result<(), Error> {
        let _guard = EnsureLogger;
        if let Some(loggers) = self.config.get("log") {
            let loggers = try!(loggers.as_table().ok_or(
                Error::from("'log' section must be a table")
            ));

            if loggers.get("target").is_some() {
                if loggers.get("targets").is_some() {
                    return Err(Error::from(
                        "Can not have both 'target' and 'targets' in the log section"
                    ));
                }

                let appender = try!(read_appender(loggers, "main"));
                let config = Config::builder()
                    .appender(appender)
                    .build(Root::builder().appender("main").build(LogLevelFilter::Trace))
                    .expect("Error in logging initialization");
                let _ = log4rs::init_config(config);
            } else if let Some(targets) = loggers.get("targets") {
                let targets = try!(targets.as_array().ok_or(
                    Error::from("'targets' must be an array in 'log' section")
                ));

                let mut appenders = Vec::new();
                for (i, target) in targets.into_iter().enumerate() {
                    let target = try!(target.as_table().ok_or(Error::from(
                        "'targets' must be an array of tables in 'log' section"
                    )));
                    let appender = try!(read_appender(target, &i.to_string()));
                    appenders.push(appender);
                }

                let mut root = Root::builder();
                let mut config = Config::builder();
                for appender in appenders {
                    root = root.appender(appender.name());
                    config = config.appender(appender);
                }
                let config = config.build(root.build(LogLevelFilter::Trace))
                                   .expect("Error in logging initialization");
                let _ = log4rs::init_config(config);
            } else {
                return Err(Error::from(
                    "Missing 'target' or 'targets' in log section"
                ));
            }
        } else {
            setup_default_logger();
            info!("Logging to console as default.");
            info!("You can configure logging in the [log] section of input file.");
        }
        Ok(())
    }
}

fn setup_default_logger() {
    // We just log everything to stdout
    let stdout = ConsoleAppender::builder()
        .target(console::Target::Stdout)
        .encoder(Box::new(LogEncoder))
        .build();

    let appender = Appender::builder()
        .filter(Box::new(ThresholdFilter::new(LogLevelFilter::Info)))
        .build("main", Box::new(stdout));

    let config = Config::builder()
        .appender(appender)
        .build(Root::builder().appender("main").build(LogLevelFilter::Info))
        .expect("Error in logging initialization");

    // We ignore the result of this call, because it can only fail if a
    // logger has already been initialized.
    let _ = log4rs::init_config(config);
}

fn read_appender(config: &Table, name: &str) -> Result<Appender, Error> {
    let allowed_keys = ["target", "targets", "level", "append"];
    for key in config.keys() {
        if !allowed_keys.contains(&&**key) {
            return Err(Error::from(format!(
                "Unknown '{}' key in log section", key
            )))
        }
    }

    let level = try!(config.get("level")
                           .map(|level| level.as_str())
                           .unwrap_or(Some("info"))
                           .ok_or(Error::from(
                               "'level' must be a string in log target"
                           )));
    let level = match level {
        "trace" => LogLevelFilter::Trace,
        "debug" => LogLevelFilter::Debug,
        "info" => LogLevelFilter::Info,
        "warning" => LogLevelFilter::Warn,
        "error" => LogLevelFilter::Error,
        other => return Err(Error::from(format!(
            "Unknown logging level '{}'", other
        )))
    };

    let target = try!(extract::str("target", config, "log target"));
    let appender: Box<Append> = match target {
        "<stdout>" => {
            Box::new(ConsoleAppender::builder()
                                     .target(console::Target::Stdout)
                                     .encoder(Box::new(LogEncoder))
                                     .build())
        }
        "<stderr>" => {
            Box::new(ConsoleAppender::builder()
                                     .target(console::Target::Stderr)
                                     .encoder(Box::new(LogEncoder))
                                     .build())
        }
        "" => return Err(Error::from(
            "'target' can not be an empty string in log target"
        )),
        filename => {
            let append_mode = try!(config.get("append")
                                         .map(|x| x.as_bool())
                                         .unwrap_or(Some(false))
                                         .ok_or(Error::from(
                                             "'append' must be a boolean in log file target"
                                         )));
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
