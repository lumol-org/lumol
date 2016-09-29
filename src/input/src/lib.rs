// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

//! Input capacities for lumol
//!
//! This module provide input files reader for two types of files:
//!
//! * Configuration files uses the TOML format to define a `System` or a
//!   `Simulation` in an human-readable way, without writing any code. This
//!   type of input is further divided in potentials input files and whole
//!   simulation input files;
//! * Structures files defines the positions and the names of particles in a
//!   `System` or in a specific `Molecule`.
//!
//! # Reading configuration files
//!
//! The `read_interactions` function read interactions into a `System`, and the
//! `read_config` function reads a whole simulation (`Simulation` and `System`
//! objects).


#![warn(
    missing_docs, trivial_casts, unused_import_braces, variant_size_differences,
    unused_qualifications, unused_results
)]

#![warn(clippy, clippy_pedantic)]
#![allow(unknown_lints)]
// List of lints we allow in this code
#![allow(
    float_arithmetic, integer_arithmetic, indexing_slicing, needless_return,
    needless_range_loop, shadow_reuse, shadow_same, shadow_unrelated,
    cast_possible_truncation, cast_precision_loss, cast_sign_loss,
    cast_possible_wrap, float_cmp, or_fun_call, string_add, non_ascii_literal,
    doc_markdown, missing_docs_in_private_items, module_inception, stutter,
    unseparated_literal_suffix
)]

#[macro_use]
extern crate log;
extern crate toml;
extern crate chemfiles;
extern crate lumol;

mod error;
pub use self::error::{Error, Result};

#[macro_use]
mod macros;

mod interactions;
pub use self::interactions::read_interactions;
pub use self::interactions::FromTomlWithPairs;

mod simulations;
pub use self::simulations::read_config;

#[cfg(test)]
pub mod testing;

#[cfg(test)]
pub use self::interactions::read_interactions_string;

use toml::Table;

fn validate(config: &Table) -> Result<()> {
    let input = try!(config.get("input").ok_or(
        Error::from("Missing 'input' table")
    ));

    let version = try!(input.lookup("version").ok_or(
        Error::from("Missing 'version' key in 'input' table")
    ));

    let version = try!(version.as_integer().ok_or(
        Error::from("'input.version' must be an integer")
    ));

    if version != 1 {
        return Err(Error::from(
            format!("Only version 1 of input can be read, got {}", version)
        ))
    }
    Ok(())
}


/// Convert a TOML table to a Rust type.
pub trait FromToml: Sized {
    /// Do the conversion from `table` to Self.
    fn from_toml(table: &Table) -> Result<Self>;
}

/// Convert a TOML table and some additional data to a Rust type.
pub trait FromTomlWithData: Sized {
    /// The type of the additional data needed.
    type Data;
    /// Do the conversion from `table` and `data` to Self.
    fn from_toml(table: &Table, data: Self::Data) -> Result<Self>;
}
