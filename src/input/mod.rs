/* Cymbalum, Molecular Simulation in Rust - Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
 */
//! Reading input files for cymbalum.
//!
//! Configuration files uses the YAML format to define the `Universe` or the
//! `Simulation` in an human-readable way, and without writing any code.
//!
//!  * The `read_interactions` function read interaction into an `Universe`.
//!  * The `molecule_from_file` function read the first molecule of the first
//!    frame of a file. This is be useful to define Monte-Carlo moves operating
//!    on a specific molecule.

mod interactions;
pub use self::interactions::read_interactions;
pub use self::interactions::Error;

mod molecules;
pub use self::molecules::molecule_from_file;
