// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

//! Simulation algorithms for lumol

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
#![allow(use_self, redundant_field_names, or_fun_call, needless_return, needless_range_loop)]
#![allow(doc_markdown, stutter, missing_docs_in_private_items, non_ascii_literal)]
#![allow(new_without_default, new_without_default_derive, should_implement_trait)]
#![allow(needless_pass_by_value, unreadable_literal, redundant_field_names, range_plus_one)]

extern crate lumol_core as core;

#[macro_use]
extern crate log;
#[macro_use]
extern crate log_once;
#[macro_use]
extern crate soa_derive;

#[cfg(test)]
#[macro_use]
extern crate approx;

extern crate caldyn;
extern crate rand;
extern crate num_traits as num;

mod propagator;
pub use self::propagator::Propagator;
pub use self::propagator::TemperatureStrategy;

pub mod output;
pub mod md;
pub mod mc;
pub mod min;

mod simulations;
pub use self::mc::MonteCarlo;
pub use self::md::MolecularDynamics;
pub use self::min::Minimization;
pub use self::simulations::Simulation;

mod velocities;
pub use self::velocities::{InitVelocities, BoltzmannVelocities, UniformVelocities};
