pub extern crate lumol_core as core;
pub extern crate lumol_input as input;
pub extern crate lumol_sim as sim;

/// The full version of the crate, containing git state if available
pub static VERSION: &'static str = env!("LUMOL_FULL_GIT_VERSION");

pub use self::core::*;
