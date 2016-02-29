// Cymbalum, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

//! Some basic types used in all the other modules

extern crate ndarray;
use self::ndarray::{OwnedArray, Ix};

/// Two dimmensional array, based on ndarray
pub type Array2<T> = OwnedArray<T, (Ix, Ix)>;
/// Three dimmensional array, based on ndarray
pub type Array3<T> = OwnedArray<T, (Ix, Ix, Ix)>;

mod vectors;
pub use self::vectors::Vector3D;

mod matrix;
pub use self::matrix::Matrix3;

mod complex;
pub use self::complex::Complex;
