// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

//! This module provides complexe numbers; 3D vectors and matrix; and
//! multidimensional arrays for use in all other modules.
mod vectors;
pub use self::vectors::Vector3D;

mod matrix;
pub use self::matrix::Matrix3;

mod complex;
pub use self::complex::Complex;

mod arrays;
pub use self::arrays::{Array2, Array3};
