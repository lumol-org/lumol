// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

#![allow(clippy::needless_doctest_main)]

//! This module provide a way to build a Lumol simulation using input files.
//!
//! Instead of building the [`System`] and [`Simulation`] objects by hand before
//! being able to use them, this module provides a way to describe the
//! simulation and the system in a TOML input file, using a simpler syntax.
//!
//! The main entry point is the `Input` struct, which allow to read a whole
//! simulation configuration. The easiest way to run a simulation using Lumol
//! is then:
//!
//! ```no_run
//! use lumol::input::Input;
//!
//! fn main() {
//!     let input = Input::new("simulation.toml").unwrap();
//!     let mut config = input.read().unwrap();
//!
//!     config.simulation.run(&mut config.system, config.nsteps);
//! }
//! ```
//!
//! This crate also provide an [`InteractionsInput`] for reading interactions
//! from a TOML file. It can be used to set the interactions in a system:
//!
//! ```no_run
//! use lumol::System;
//! use lumol::input::InteractionsInput;
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
//! [`System`]: ../sys/struct.System.html
//! [`Simulation`]: ../sim/struct.Simulation.html
//! [`InteractionsInput`]: struct.InteractionsInput.html

pub use lumol_input::*;
