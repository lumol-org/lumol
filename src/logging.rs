// Cymbalum, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

//! Logging configuration
//!
//! This module provides a `Logger` struct which can be used for logging
//! messages from Cymbalum. One and only one of the initilialization function
//! should be called, and it should only be called once. Any call to
//! initializing function after the first one will be ignored.
//!
//! Other [log](https://github.com/rust-lang/log) compatible crates can be used
//! for logging if you whant.

use std::io;
use std::io::prelude::*;
use std::fs::OpenOptions;
use std::path::Path;
use std::sync::{Arc, Mutex, RwLock};

use log::{Log, LogRecord, LogMetadata, set_logger};
pub use log::LogLevel;

/// Log an error, and then panic with the same message
macro_rules! fatal_error {
    ($($args:tt)*) => {{
        error!($($args)*);
        panic!($($args)*);
    }};
}

macro_rules! internal_error {
    ($arg: expr) => {{
        let message = format!("Internal error: {}", $arg);
        fatal_error!("{}", message);
    }};
}


macro_rules! warn_once {
    ($($args:tt)*) => {{
        use ::std::collections::BTreeSet;
        use ::std::sync::RwLock;
        lazy_static!{
            static ref WARNINGS: RwLock<BTreeSet<String>> = RwLock::new(BTreeSet::new());
        }
        let message = format!($($args)*);
        let already_logged = WARNINGS.read()
                                     .expect("RwLock was poisonned")
                                     .contains(&message);
        if !already_logged {
            WARNINGS.write()
                    .expect("RwLock was poisonned")
                    .insert(message.clone());
            warn!("{}", message);
        }
    }};
}

/// Target for the logging output
#[derive(Clone, Debug)]
pub enum Target {
    /// Log are written to the standard output stream
    StdOut,
    /// Log are written to the standard error stream
    StdErr,
    /// Log are written to a file, the path to the file is in the `String`
    File(String),
    /// Logging is not initialized yet
    None
}

/// Logger with capacity to write to the standard output stream, the standard
/// error stream or a file.
pub struct Logger {
    level: LogLevel,
    target: Target,
    writer: Arc<Mutex<Box<Write + Send + Sync>>>,
}

lazy_static!(
    static ref TARGET: RwLock<Target> = RwLock::new(Target::None);
);

impl Logger{
    fn new<T>(level: LogLevel, target: Target, handle: T) -> Logger where T: Write + Send + Sync + 'static {
        ::chemfiles::Logger::get().log_callback(|level, message| {
            match level {
                ::chemfiles::LogLevel::Error => error!("{}", message),
                ::chemfiles::LogLevel::Warning => warn!("{}", message),
                ::chemfiles::LogLevel::Info => info!("{}", message),
                ::chemfiles::LogLevel::Debug => debug!("{}", message),
            }
        }).expect("Could not set the chemfiles logging callback");
        Logger {
            level: level,
            target: target,
            writer: Arc::new(Mutex::new(Box::new(handle))),
        }
    }

    fn init(logger: Logger) {
        let res = set_logger(|max_log_level| {
            max_log_level.set(logger.level.to_log_level_filter());
            let mut target = TARGET.write().expect("Logging TARGET lock is poisonned");
            *target = logger.target.clone();
            Box::new(logger)
        });

        if let Err(_) = res {
            warn!("Logger was initialized more than once! This call is ignored.");
        }
    }

    /// Initialize the global logger to use the standard output stream.
    pub fn stdout() {
        Logger::stdout_at_level(levels::DEFAULT_LEVEL)
    }

    /// Initialize the global logger to use the standard output stream, with
    /// the maximum log level of `level`
    pub fn stdout_at_level(level: LogLevel) {
        let logger = Logger::new(level, Target::StdOut, io::stdout());
        Logger::init(logger);
    }

    /// Initialize the global logger to use the standard error stream.
    pub fn stderr() {
        Logger::stderr_at_level(levels::DEFAULT_LEVEL)
    }

    /// Initialize the global logger to use the standard error stream, with
    /// the maximum log level of `level`
    pub fn stderr_at_level(level: LogLevel) {
        let logger = Logger::new(level, Target::StdErr, io::stderr());
        Logger::init(logger);
    }

    /// Initialize the global logger to write to the file at `path`, creating
    /// the file if needed. If the file already exists, it is not removed.
    pub fn file<P: AsRef<Path>>(path: P) -> Result<(), io::Error> {
        Logger::file_at_level(path, levels::DEFAULT_LEVEL)
    }

    /// Initialize the global logger to write to the file at `path`, creating
    /// the file if needed. If the file already exists, it is not removed. The
    /// maximum log level is set to `level`.
    pub fn file_at_level<P: AsRef<Path>>(path: P, level: LogLevel) -> Result<(), io::Error> {
        let filename = format!("{}", path.as_ref().display());
        let file = try!(OpenOptions::new().write(true).create(true).append(true).open(path));
        let logger = Logger::new(level, Target::File(filename), file);
        Logger::init(logger);
        Ok(())
    }

    /// Get the target for log events
    pub fn target() -> Target {
        TARGET.read().expect("Logging TARGET lock is poisonned").clone()
    }

    /// Is the log events target the screen?
    pub fn is_screen() -> bool {
        match *TARGET.read().expect("Logging TARGET lock is poisonned") {
            Target::StdErr | Target::StdOut => true,
            Target::File(_) | Target::None  => false
        }
    }

    /// Is the log events target a file?
    pub fn is_file() -> bool {
        match *TARGET.read().expect("Logging TARGET lock is poisonned") {
            Target::File(_) => true,
            Target::StdErr | Target::StdOut | Target::None => false
        }
    }
 }


#[cfg(not(debug_assertions))]
 mod levels {
     use super::LogLevel;
     /// Default log level
     pub const DEFAULT_LEVEL: LogLevel = LogLevel::Warn;
 }

#[cfg(debug_assertions)]
 mod levels {
     use super::LogLevel;
     /// Default log level
     pub const DEFAULT_LEVEL: LogLevel = LogLevel::Debug;
 }

/// Implements Log trait for Logger
impl Log for Logger {
    fn enabled(&self, metadata: &LogMetadata) -> bool {
        metadata.level() <= self.level
    }

    fn log(&self, record: &LogRecord) {
        let mut out = self.writer.lock().expect("Could not lock the logger.");
        let write_res = match record.level() {
            LogLevel::Info => write!(&mut out, "{}\n", record.args()),
            LogLevel::Warn => write!(
                   &mut out, "Warning: {} from {}:{}\n",
                   record.args(),
                   record.location().file(), record.location().line()
            ),
            _ => write!(
                   &mut out, "{}: {} from {}:{}\n",
                   record.level(), record.args(),
                   record.location().file(), record.location().line()
            ),
        };


        if let Err(err) = write_res {
            let error_res = write!(& mut io::stderr(),
                "Error while writing log: {}\nMessage was: {} at {}:{} ({})\n",
                err, record.args(),
                record.location().file(), record.location().line(),
                record.level()
            );
            if let Err(err2) = error_res {
                panic!(format!(
                    "Error while logging logger error, aborting. Error chain: {} -> {} from {}:{}",
                    err, err2,
                    record.location().file(), record.location().line(),
                ))
            }
        }
    }
}
