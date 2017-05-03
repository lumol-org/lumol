// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 Lumol's contributors — BSD license

//! The `system` module provide a way to store data about a simulated system.

mod config;
pub use self::config::*;

mod system;
pub use self::system::System;
