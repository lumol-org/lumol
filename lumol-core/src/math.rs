// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

//! Access usual math function directly, without having to use a `f64::` prefix,
//! or to resort to method style call.
#![allow(clippy::inline_always)]

use special::Error;

macro_rules! make_math_fn {
    ($name: ident) => (
        #[inline(always)]
        pub fn $name(value: f64) -> f64 {
            f64::$name(value)
        }
    );
}

make_math_fn!(sqrt);
make_math_fn!(exp);
make_math_fn!(abs);
make_math_fn!(cos);
make_math_fn!(sin);
make_math_fn!(acos);
make_math_fn!(floor);
make_math_fn!(round);

#[inline(always)]
pub fn erf(value: f64) -> f64 {
    f64::error(value)
}

#[inline(always)]
pub fn erfc(value: f64) -> f64 {
    f64::compl_error(value)
}
