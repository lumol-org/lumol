// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

//! Core types and definitions for lumol.

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
#![allow(clippy::type_complexity)]

#![allow(clippy::missing_errors_doc)]

// Tests lints
#![cfg_attr(test, allow(clippy::float_cmp))]

// deny(warnings) in doc tests
#![doc(test(attr(deny(warnings))))]
#![doc(test(attr(allow(unused_variables))))]

// Helper modules
#[macro_use]
mod utils;
mod math;

// Main modules
pub mod units;
pub mod consts;
pub mod types;
pub mod energy;
pub mod sys;

pub use self::types::*;
pub use self::energy::*;
pub use self::sys::*;
