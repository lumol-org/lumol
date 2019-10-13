// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

//! Lumol is a classical molecular simulation engine that provides a solid
//! base for developing new algorithms and methods.
//!
//! Using Lumol, you can customize the behavior of all the algorithms in a
//! simulation (from force fields to barostats and Monte Carlo moves).
//!
//! Lumol goals are to be:
//!
//! - **Easy to extend**: the code is modular, object-oriented, well documented,
//!   well tested, open-source and readable;
//! - **Easy to use**: the user interface is nice, with human-oriented input
//!   files;
//! - **Stable**: it will never crash on a good input, and provides helpful
//!   error messages.

/// The full version of the crate, containing git state if available
pub static VERSION: &str = env!("LUMOL_FULL_GIT_VERSION");

pub mod sim;
pub mod input;

pub use lumol_core::*;
