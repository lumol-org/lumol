// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 Lumol's contributors — BSD license

//! Various internal utilities, which do not have there own module
#[macro_use]
mod macros;

mod parallel_shortcuts;
pub use self::parallel_shortcuts::ParallelShortcuts;

#[cfg(test)]
mod xyz;
#[cfg(test)]
pub use self::xyz::system_from_xyz;


/// Internal version of `units::from`, where the unit is assumed to be correct
pub fn unit_from(value: f64, unit: &str) -> f64 {
    ::units::from(value, unit).expect("Internal unit error. This is a bug.")
}

/// Internal version of `units::to`, where the unit is assumed to be correct
pub fn unit_to(value: f64, unit: &str) -> f64 {
    ::units::to(value, unit).expect("Internal unit error. This is a bug.")
}
