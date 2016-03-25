// Cymbalum, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

//! Cymbalum is a classical molecular simulation engine that provides a solid
//! base for developing new algorithms and methods.
//!
//! Using Cymbalum, you can customize the behavior of all the algorithms in a
//! simulation (from force fields to barostats and Monte-Carlo moves).
//!
//! Cymbalum goals are to be:
//!
//! - **Easy to extend**: the code is modular, object-oriented, well documented,
//!   well tested, open-source and readable;
//! - **Easy to use**: the user interface is nice, with human-oriented input
//!   files;
//! - **Stable**: it will never crash on a good input, and provides helpful
//!   error messages.

#![warn(
    missing_docs, trivial_casts, unused_import_braces, variant_size_differences,
    unused_qualifications, unused_results
)]

#![cfg_attr(feature="lint", plugin(herbie_lint))]
#![cfg_attr(feature="lint", feature(plugin))]
#![cfg_attr(feature="lint", plugin(clippy))]
#![cfg_attr(feature="lint", warn(clippy))]
// Additional lints from the Allow group in clippy
#![cfg_attr(feature="lint", warn(
    enum_glob_use, mut_mut, option_unwrap_used, print_stdout,
    result_unwrap_used, single_match_else, wrong_pub_self_convention
))]
// These are for readability
#![cfg_attr(feature="lint", allow(
    needless_return, needless_range_loop, or_fun_call, new_without_default,
    float_cmp,
))]

#[macro_use]
extern crate log;
#[macro_use]
extern crate bitflags;
#[macro_use]
extern crate lazy_static;

extern crate chemfiles;
extern crate ndarray;
extern crate num;
extern crate rand;
extern crate special;
extern crate yaml_rust as yaml;

#[macro_use]
mod tests;
#[macro_use]
mod utils;
#[macro_use]
mod logging;
pub use logging::{Logger, LogLevel};

pub mod units;
pub mod constants;

pub mod types;
pub mod potentials;
pub mod system;
pub mod simulation;

pub use types::*;
pub use potentials::*;
pub use system::*;
pub use simulation::*;

pub mod input;
