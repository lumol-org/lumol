/* Cymbalum, Molecular Simulation in Rust - Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
 */

//! **cymbalum** is a molecular simulation library, which provides basic
//! building blocks to create your very own simulations. It is based on two main
//! ideas:
//!
//!  - any algorithm used in the simulation can be replaced by another one. This
//!    allow for modularity and easy developement of novel algorithms.
//!  - data and algorithms should be separated. This is accomplished through the
//!    use of the `Universe` (for data) and `Simulation` (for algorithms) types,
//!    interacting together to run the simulation on some data.
//!
//!  Rust provides a nice way to implement these two ideas with the concept of
//!  traits.

#![allow(non_snake_case)]
#![warn(
    missing_docs,
    trivial_casts,
    trivial_numeric_casts,
    unused_import_braces,
    variant_size_differences,
    unused_qualifications
)]

#[macro_use]
extern crate log;

#[macro_use]
mod tests;

mod logging;
pub use logging::Logger;

pub mod units;
pub mod constants;

pub mod types;
pub mod potentials;
pub mod universe;
pub mod simulation;

pub use types::*;
pub use potentials::*;
pub use universe::*;
pub use simulation::*;
