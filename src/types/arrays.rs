// Cymbalum, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

//! Multi-dimmensional arrays based on ndarray
extern crate ndarray;
use self::ndarray::{OwnedArray, Ix};

use std::ops::{Index, IndexMut};
use std::cmp::PartialEq;
use types::Zero;

/// Two dimmensional array, based on ndarray.
///
/// Most of the methods are simply forwarded to ndarray, so also look the doc
/// for this crate.
#[derive(Debug, Clone)]
pub struct Array2<T>(OwnedArray<T, (Ix, Ix)>);

impl<T: Zero + Clone> Array2<T> {
    /// Create a new `Array2` of the specified `size` filled with the
    /// `Zero::zero` return value.
    pub fn zeros(size: (Ix, Ix)) -> Array2<T> {
        Array2(OwnedArray::<T, (Ix, Ix)>::zeros(size))
    }

    /// Resize the array if the current size is not `size`, and fill the
    /// new array with zeros.
    pub fn resize_if_different(&mut self, size: (Ix, Ix)) {
        if self.0.shape() != &[size.0, size.1] {
            *self = Array2::zeros(size);
        }
    }
}

impl<T> Array2<T> {
    /// Get the shape of the array
    pub fn shape(&self) -> (Ix, Ix) {
        let shape = self.0.shape();
        (shape[0], shape[1])
    }
}

impl<T: Clone> Array2<T> {
    /// Assign the given scalar to all entry in this array
    pub fn assign(&mut self, value: T) {
        self.0.assign_scalar(&value);
    }
}

impl<T: Default> Array2<T> {
    /// Create a new `Array2` of the specified `size` filled with the
    /// `Default::default` return value.
    pub fn default(size: (Ix, Ix)) -> Array2<T>{
        Array2(OwnedArray::<T, (Ix, Ix)>::default(size))
    }
}

impl<T: PartialEq> PartialEq for Array2<T> {
    fn eq(&self, other: &Array2<T>) -> bool {
        self == other
    }
}

impl<T> Index<(Ix, Ix)> for Array2<T> {
    type Output = T;
    fn index(&self, index: (Ix, Ix)) -> &T {
        unsafe {
            // ndarray does the check for us in debug builds
            self.0.uget(index)
        }
    }
}

impl<T> IndexMut<(Ix, Ix)> for Array2<T> {
    fn index_mut(&mut self, index: (Ix, Ix)) -> &mut T {
        unsafe {
            // ndarray does the check for us in debug builds
            self.0.uget_mut(index)
        }
    }
}

/******************************************************************************/

/// Three dimmensional array, based on ndarray
///
/// Most of the methods are simply forwarded to ndarray, so also look the doc
/// for this crate.
#[derive(Debug, Clone)]
pub struct Array3<T>(OwnedArray<T, (Ix, Ix, Ix)>);

impl<T: Zero + Clone> Array3<T> {
    /// Create a new `Array3` of the specified `size` filled with the
    /// `Zero::zero` return value.
    pub fn zeros(size: (Ix, Ix, Ix)) -> Array3<T> {
        Array3(OwnedArray::<T, (Ix, Ix, Ix)>::zeros(size))
    }

    /// Resize the array if the current size is not `size`, and fill the
    /// new array with zeros.
    pub fn resize_if_different(&mut self, size: (Ix, Ix, Ix)) {
        if self.0.shape() != &[size.0, size.1, size.2] {
            *self = Array3::zeros(size);
        }
    }
}

impl<T> Array3<T> {
    /// Get the shape of the array
    pub fn shape(&self) -> (Ix, Ix, Ix) {
        let shape = self.0.shape();
        (shape[0], shape[1], shape[2])
    }
}

impl<T: Clone> Array3<T> {
    /// Assign the given scalar to all entry in this array
    pub fn assign(&mut self, value: T) {
        self.0.assign_scalar(&value);
    }
}

impl<T: Default> Array3<T> {
    /// Create a new `Array3` of the specified `size` filled with the
    /// `Default::default` return value.
    pub fn default(size: (Ix, Ix, Ix)) -> Array3<T>{
        Array3(OwnedArray::<T, (Ix, Ix, Ix)>::default(size))
    }
}

impl<T> Index<(Ix, Ix, Ix)> for Array3<T> {
    type Output = T;
    fn index(&self, index: (Ix, Ix, Ix)) -> &T {
        unsafe {
            // ndarray does the check for us in debug builds
            self.0.uget(index)
        }
    }
}

impl<T> IndexMut<(Ix, Ix, Ix)> for Array3<T> {
    fn index_mut(&mut self, index: (Ix, Ix, Ix)) -> &mut T {
        unsafe {
            // ndarray does the check for us in debug builds
            self.0.uget_mut(index)
        }
    }
}

#[cfg(test)]
mod tests {
    mod array2 {
        use super::super::Array2;

        #[test]
        fn zeros() {
            let a: Array2<f64> = Array2::zeros((3, 4));
            for i in 0..3 {
                for j in 0..4 {
                    assert_eq!(a[(i, j)], 0.0);
                }
            }
        }

        #[test]
        fn default() {
            #[derive(Clone)] struct F64(f64);
            impl Default for F64 {fn default() -> F64 { F64(42.0) }}

            let a: Array2<F64> = Array2::default((3, 4));
            for i in 0..3 {
                for j in 0..4 {
                    assert_eq!(a[(i, j)].0, 42.0);
                }
            }
        }

        #[test]
        fn resize() {
            let mut a: Array2<f64> = Array2::zeros((3, 4));
            assert_eq!(a.shape(), (3, 4));
            a[(1, 1)] = 42.0;

            a.resize_if_different((7, 90));
            assert_eq!(a.shape(), (7, 90));
            assert_eq!(a[(1, 1)], 0.0);

            a[(1, 1)] = 42.0;
            a.resize_if_different((7, 90));
            assert_eq!(a[(1, 1)], 42.0);
        }

        #[test]
        fn index() {
            let mut a: Array2<f64> = Array2::zeros((3, 4));

            assert_eq!(a[(1, 3)], 0.0);

            a[(1, 3)] = 45.0;
            assert_eq!(a[(1, 3)], 45.0);
        }

        #[test]
        #[should_panic]
        #[cfg(debug_assertions)]
        fn out_of_bound_1() {
            let a: Array2<f64> = Array2::zeros((3, 4));
            let _ = a[(5, 1)];
        }

        #[test]
        #[should_panic]
        #[cfg(debug_assertions)]
        fn out_of_bound_2() {
            let a: Array2<f64> = Array2::zeros((3, 4));
            let _ = a[(2, 7)];
        }
    }

    mod array3 {
        use super::super::Array3;

        #[test]
        fn zeros() {
            let a: Array3<f64> = Array3::zeros((3, 4, 8));
            for i in 0..3 {
                for j in 0..4 {
                    for k in 0..8 {
                        assert_eq!(a[(i, j, k)], 0.0);
                    }
                }
            }
        }

        #[test]
        fn default() {
            #[derive(Clone)] struct F64(f64);
            impl Default for F64 {fn default() -> F64 { F64(42.0) }}

            let a: Array3<F64> = Array3::default((3, 4, 8));
            for i in 0..3 {
                for j in 0..4 {
                    for k in 0..8 {
                        assert_eq!(a[(i, j, k)].0, 42.0);
                    }
                }
            }
        }

        #[test]
        fn resize() {
            let mut a: Array3<f64> = Array3::zeros((3, 4, 5));
            assert_eq!(a.shape(), (3, 4, 5));
            a[(1, 1, 1)] = 42.0;

            a.resize_if_different((7, 90, 8));
            assert_eq!(a.shape(), (7, 90, 8));
            assert_eq!(a[(1, 1, 1)], 0.0);

            a[(1, 1, 1)] = 42.0;
            a.resize_if_different((7, 90, 8));
            assert_eq!(a[(1, 1, 1)], 42.0);
        }

        #[test]
        fn index() {
            let mut a: Array3<f64> = Array3::zeros((3, 4, 5));
            assert_eq!(a[(1, 3, 2)], 0.0);

            a[(1, 3, 2)] = 45.0;
            assert_eq!(a[(1, 3, 2)], 45.0);
        }

        #[test]
        #[should_panic]
        #[cfg(debug_assertions)]
        fn out_of_bound_1() {
            let a: Array3<f64> = Array3::zeros((3, 4, 89));
            let _ = a[(5, 1, 6)];
        }

        #[test]
        #[should_panic]
        #[cfg(debug_assertions)]
        fn out_of_bound_2() {
            let a: Array3<f64> = Array3::zeros((3, 4, 89));
            let _ = a[(2, 67, 6)];
        }

        #[test]
        #[should_panic]
        #[cfg(debug_assertions)]
        fn out_of_bound_3() {
            let a: Array3<f64> = Array3::zeros((3, 4, 89));
            let _ = a[(2, 1, 600)];
        }
    }
}
