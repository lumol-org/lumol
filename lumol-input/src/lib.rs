// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

//! This crate provide a way to build a Lumol simulation using input files.
//!
//! Instead of building the `System` object, the `Simulation` object by hand
//! before being able to use them to run the simulation, this crate allow to
//! describe the simulation and the system in a TOML input file, using a simple
//! syntax.
//!
//! The main entry point of this crate is the `Input` struct, which allow to
//! read a whole simulation configuration. The easiest way to run a simulation
//! using Lumol is:
//!
//! ```no_run
//! extern crate lumol_input;
//! use lumol_input::Input;
//!
//! fn main() {
//!     let input = Input::new("simulation.toml").unwrap();
//!     let mut config = input.read().unwrap();
//!
//!     config.simulation.run(&mut config.system, config.nsteps);
//! }
//! ```
//!
//! This crate also provide an `InteractionsInput` for reading interactions
//! from a TOML file. It can be used to set the interactions in a system:
//!
//! ```no_run
//! extern crate lumol_core;
//! extern crate lumol_input;
//! use lumol_core::sys::System;
//! use lumol_input::InteractionsInput;
//!
//! fn main() {
//!     let mut system = System::new();
//!
//!     // ... Build the system by hand
//!
//!     // Read the interactions
//!     let input = InteractionsInput::new("potentials.toml").unwrap();
//!     input.read(&mut system).unwrap();
//! }
//! ```
//!
#![warn(missing_docs, trivial_casts, unused_import_braces, variant_size_differences)]
#![warn(unused_qualifications, unused_results)]
// Clippy configuration
#![allow(unknown_lints)]
#![warn(clippy, clippy_pedantic)]
// Not embed software, integer and float arithmeric are allowed
#![allow(float_arithmetic, integer_arithmetic, indexing_slicing)]
// Cast issues
#![allow(cast_possible_truncation, cast_precision_loss, cast_sign_loss, cast_possible_wrap)]
// Style issues
#![allow(shadow_reuse, shadow_same, shadow_unrelated)]
#![allow(use_self, redundant_field_names, or_fun_call, needless_return)]
#![allow(missing_docs_in_private_items, should_implement_trait)]

extern crate chemfiles;
extern crate lumol_core as lumol;
extern crate toml;

extern crate log4rs;
#[macro_use]
extern crate log;

use toml::value::Table;

macro_rules! try_io {
    ($expr: expr, $path: expr) => (
        match $expr {
            Ok(val) => val,
            Err(err) => {
                return Err(Error::from((err, $path)));
            }
        }
    );
}

mod extract;
mod error;
mod interactions;
mod simulations;

pub use self::error::{Error, Result};
pub use self::interactions::Input as InteractionsInput;
pub use self::simulations::{Config, Input};
pub use self::simulations::setup_default_logger;

/// Convert a TOML table to a Rust type.
pub trait FromToml: Sized {
    /// Do the conversion from `table` to Self.
    fn from_toml(table: &Table) -> Result<Self>;
}

/// Convert a TOML table and some additional owned data to a Rust type.
pub trait FromTomlWithData: Sized {
    /// The type of the additional data needed.
    type Data;
    /// Do the conversion from `table` and `data` to Self.
    fn from_toml(table: &Table, data: Self::Data) -> Result<Self>;
}

/// Convert a TOML table to a Rust type using information from an additional reference.
pub trait FromTomlWithRefData: Sized {
    /// The type of the additional data needed.
    type Data;
    /// Do the conversion from `table` and `data` to Self.
    fn from_toml(table: &Table, data: &Self::Data) -> Result<Self>;
}

fn validate(config: &Table) -> Result<()> {
    let input = config.get("input").ok_or(
        Error::from("Missing 'input' table")
    )?;

    let version = input.get("version").ok_or(
        Error::from("Missing 'version' key in 'input' table")
    )?;

    let version = version.as_integer().ok_or(
        Error::from("'input.version' must be an integer")
    )?;

    if version != 1 {
        return Err(Error::from(
            format!("Can only read version 1 of input, got version {}", version),
        ));
    }
    Ok(())
}
