// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

//! Core types and definitions for lumol.

#![warn(missing_docs, trivial_casts, unused_import_braces, variant_size_differences)]
#![warn(unused_qualifications, unused_results)]
// Clippy configuration
#![warn(clippy::all, clippy::pedantic)]
// Not embed software, integer and float arithmeric are allowed
#![allow(clippy::float_arithmetic, clippy::integer_arithmetic, clippy::indexing_slicing)]
// Cast issues
#![allow(clippy::cast_possible_truncation, clippy::cast_precision_loss)]
#![allow(clippy::cast_sign_loss, clippy::cast_possible_wrap)]
// Style issues
#![allow(clippy::shadow_reuse, clippy::shadow_same, clippy::shadow_unrelated)]
#![allow(clippy::use_self, clippy::redundant_field_names, clippy::or_fun_call)]
#![allow(clippy::needless_return, clippy::needless_range_loop, clippy::doc_markdown)]
#![allow(clippy::missing_docs_in_private_items, clippy::non_ascii_literal)]
#![allow(clippy::new_without_default, clippy::new_without_default_derive)]
#![allow(clippy::should_implement_trait, clippy::needless_pass_by_value)]
#![allow(clippy::unreadable_literal, clippy::redundant_field_names)]
#![allow(clippy::range_plus_one, clippy::filter_map, clippy::stutter)]

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
