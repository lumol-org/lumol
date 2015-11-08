/* Cymbalum, Molecular Simulation in Rust - Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
 */

//! Some basic types used in all the other modules

mod vectors;
pub use self::vectors::Vector3D;

mod matrix;
pub use self::matrix::Matrix3;

mod arrays;
pub use self::arrays::{Array2, Array3};
