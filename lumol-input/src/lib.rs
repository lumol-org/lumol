// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

//! Input system for lumol using TOML as a language

#![warn(missing_docs, trivial_casts, unused_import_braces, variant_size_differences)]
#![warn(unused_qualifications, unused_results, rust_2018_idioms)]
// Clippy configuration
#![warn(clippy::all, clippy::pedantic)]
// Not embed software, integer and float arithmetic are allowed
#![allow(clippy::float_arithmetic, clippy::indexing_slicing)]
// Cast issues
#![allow(clippy::cast_possible_truncation, clippy::cast_precision_loss)]
#![allow(clippy::cast_sign_loss, clippy::cast_possible_wrap)]
// Style issues
#![allow(clippy::shadow_reuse, clippy::shadow_same, clippy::shadow_unrelated)]
#![allow(clippy::use_self, clippy::redundant_field_names, clippy::or_fun_call)]
#![allow(clippy::needless_return, clippy::needless_range_loop, clippy::doc_markdown)]
#![allow(clippy::missing_docs_in_private_items, clippy::module_name_repetitions)]
#![allow(clippy::new_without_default, clippy::range_plus_one, clippy::missing_panics_doc)]
#![allow(clippy::if_not_else, clippy::redundant_closure_for_method_calls)]
#![allow(clippy::must_use_candidate, clippy::return_self_not_must_use)]

#![allow(clippy::missing_errors_doc)]

// deny(warnings) in doc tests
#![doc(test(attr(deny(warnings))))]
#![doc(test(attr(allow(unused_variables))))]

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
mod alternator;

pub use self::error::Error;
pub use self::interactions::InteractionsInput;
pub use self::simulations::{Config, Input};
pub use self::simulations::setup_default_logger;

/// Convert a TOML table to a Rust type.
pub trait FromToml: Sized {
    /// Do the conversion from `table` to Self.
    fn from_toml(table: &Table) -> Result<Self, Error>;
}

/// Convert a TOML table and some additional owned data to a Rust type.
pub trait FromTomlWithData: Sized {
    /// The type of the additional data needed.
    type Data;
    /// Do the conversion from `table` and `data` to Self.
    fn from_toml(table: &Table, data: Self::Data) -> Result<Self, Error>;
}

/// Convert a TOML table to a Rust type using information from an additional reference.
pub trait FromTomlWithRefData: Sized {
    /// The type of the additional data needed.
    type Data;
    /// Do the conversion from `table` and `data` to Self.
    fn from_toml(table: &Table, data: &Self::Data) -> Result<Self, Error>;
}

fn validate(config: &Table) -> Result<(), Error> {
    let input = config.get("input").ok_or(
        Error::from("missing 'input' table")
    )?;

    let version = input.get("version").ok_or(
        Error::from("missing 'version' key in 'input' table")
    )?;

    let version = version.as_integer().ok_or(
        Error::from("'input.version' must be an integer")
    )?;

    if version != 1 {
        return Err(Error::from(
            format!("can only read version 1 of input, got version {version}"),
        ));
    }
    Ok(())
}
