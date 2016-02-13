// Cymbalum, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

//! Reading input files for cymbalum.
//!
//! Configuration files uses the YAML format to define the `System` or the
//! `Simulation` in an human-readable way, and without writing any code.
//!
//!  * The `read_interactions` function read interaction into an `System`.
//!  * The `molecule_from_file` function read the first molecule of the first
//!    frame of a file. This is be useful to define Monte-Carlo moves operating
//!    on a specific molecule.
mod interactions;
pub use self::interactions::{read_interactions, read_interactions_string};
pub use self::interactions::Error;

mod molecules;
pub use self::molecules::molecule_from_file;
