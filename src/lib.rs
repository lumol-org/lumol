extern crate lumol_core;

/// The full version of the crate, containing git state if available
pub static VERSION: &'static str = env!("LUMOL_FULL_GIT_VERSION");

pub use lumol_core::*;
