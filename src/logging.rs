/* Cymbalum, Molecular Simulation in Rust - Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
 */

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
use std::sync::{Arc, Mutex};

extern crate log;
use log::{Log, LogRecord, LogMetadata, set_logger};
pub use log::LogLevel;

/// Logger with capacity to write to the standard output stream, the standard
/// error stream or a file.
 pub struct Logger {
     level: LogLevel,
     writer: Arc<Mutex<Box<Write + Send + Sync>>>,
 }

 impl Logger{
     fn new<T>(level: LogLevel, handle: T) -> Logger where T: Write + Send + Sync + 'static {
         Logger{
             level: level,
             writer: Arc::new(Mutex::new(Box::new(handle))),
         }
     }

     fn init(logger: Logger) {
         let res = set_logger(|max_log_level| {
             max_log_level.set(logger.level.to_log_level_filter());
             Box::new(logger)
         });

         if let Err(_) = res {
             warn!("Logger was initialized more than once! This call is ignored.");
         }
     }

     /// Initialize the global logger to use the standard output stream.
     pub fn stdout() {
         let logger = Logger::new(levels::DEFAULT_LEVEL, io::stdout());
         Logger::init(logger);
     }

     /// Initialize the global logger to use the standard output stream, with
     /// the maximum log level of `level`
     pub fn stdout_at_level(level: LogLevel) {
         let logger = Logger::new(level, io::stdout());
         Logger::init(logger);
     }

     /// Initialize the global logger to use the standard error stream.
     pub fn stderr() {
         let logger = Logger::new(levels::DEFAULT_LEVEL, io::stderr());
         Logger::init(logger);
     }

     /// Initialize the global logger to use the standard error stream, with
     /// the maximum log level of `level`
     pub fn stderr_at_level(level: LogLevel) {
         let logger = Logger::new(level, io::stderr());
         Logger::init(logger);
     }

     /// Initialize the global logger to write to the file at `path`, creating
     /// the file if needed. If the file already exists, it is not removed.
     pub fn file<P: AsRef<Path>>(path: P) -> Result<(), io::Error> {
         let file = try!(OpenOptions::new().write(true).append(true).open(path));
         let logger = Logger::new(levels::DEFAULT_LEVEL, file);
         Logger::init(logger);
         Ok(())
     }

     /// Initialize the global logger to write to the file at `path`, creating
     /// the file if needed. If the file already exists, it is not removed. The
     /// maximum log level is set to `level`.
     pub fn file_at_level<P: AsRef<Path>>(path: P, level: LogLevel) -> Result<(), io::Error> {
         let file = try!(OpenOptions::new().write(true).append(true).open(path));
         let logger = Logger::new(level, file);
         Logger::init(logger);
         Ok(())
     }
 }


#[cfg(not(debug_assertions))]
 mod levels {
     use super::log::LogLevel;
     /// Default log level
     pub const DEFAULT_LEVEL: LogLevel = LogLevel::Warn;
 }

#[cfg(debug_assertions)]
 mod levels {
     use super::log::LogLevel;
     /// Default log level
     pub const DEFAULT_LEVEL: LogLevel = LogLevel::Debug;
 }

/// Implements Log trait for Logger
impl Log for Logger {
    fn enabled(&self, metadata: &LogMetadata) -> bool {
        metadata.level() <= self.level
    }

    fn log(&self, record: &LogRecord) {
        let mut out = self.writer.lock().ok().expect("Could not lock the logger.");
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
                "Error while writing log: {:?}\nMessage was: {} at {}:{} ({})\n",
                err, record.args(),
                record.location().file(), record.location().line(),
                record.level()
            );
            if let Err(err2) = error_res {
                panic!(format!(
                    "Error while logging logger error, aborting. Error chain: {:?} -> {:?} from {}:{}",
                    err, err2,
                    record.location().file(), record.location().line(),
                ))
            }
        }
    }
}
