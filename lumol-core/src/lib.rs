// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

//! Core types and definitions for lumol.

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

// deny(warnings) in doc tests
#![doc(test(attr(deny(warnings))))]
#![doc(test(attr(allow(unused_variables))))]

#[macro_use]
extern crate bitflags;
#[macro_use]
extern crate lazy_static;
#[macro_use]
extern crate log;
#[macro_use]
extern crate log_once;
#[macro_use]
extern crate soa_derive;

#[cfg(test)]
#[macro_use]
extern crate approx;

extern crate chemfiles;
extern crate ndarray;
extern crate num_traits as num;
extern crate rayon;
extern crate special;
extern crate thread_local;

macro_rules! zip {
    (@map $pattern:pat => $tuple:expr) => {
        |$pattern| $tuple
    };
    (@map $pattern:pat => ( $($tuple:tt)* ) , $_iter:expr $(, $tail:expr )*) => {
        zip!(@map ($pattern, b) => ( $($tuple)*, b ) $( , $tail )*)
    };
    ($first:expr $( , $rest:expr )* $(,)*) => {
        ::std::iter::IntoIterator::into_iter($first)
            $(.zip($rest))*
            .map(
                zip!(@map a => (a) $( , $rest )*)
            )
    };
}

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
