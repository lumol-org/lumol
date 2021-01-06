// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

//! Access usual math function directly, without having to use a `f64::` prefix,
//! or to resort to method style call.
#![allow(clippy::inline_always)]

use special::Error;

#[inline(always)]
pub fn erf(value: f64) -> f64 {
    f64::error(value)
}

#[inline(always)]
pub fn erfc(value: f64) -> f64 {
    f64::compl_error(value)
}
