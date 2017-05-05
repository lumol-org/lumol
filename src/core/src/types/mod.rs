// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 Lumol's contributors — BSD license

//! Linear algebra types for Lumol.
//!
//! This module provides vectors, matrix, 2D and 3D tensors for use in all other
//! modules.
pub use num::{One, Zero};

mod vectors;
pub use self::vectors::Vector3D;

mod matrix;
pub use self::matrix::Matrix3;

mod complex;
pub use self::complex::Complex;

mod arrays;
pub use self::arrays::{Array2, Array3};
