// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

//! Multi-dimensional arrays based on ndarray
use ndarray;

use std::ops::{Deref, DerefMut, Index, IndexMut};
use types::Zero;

/// Two dimensional tensors, based on ndarray.
///
/// Most of the methods are simply forwarded to ndarray, so also look the doc
/// for this crate. This array type mainly supports indexing using tuples as
/// indices and is though as storage backend for multi-dimensional data.
///
/// ```
/// # use lumol_core::types::Array2;
/// let mut a = Array2::zeros((3, 5));
///
/// assert_eq!(a[(0, 4)], 0.0);
///
/// a[(0, 4)] = 7.0;
/// assert_eq!(a[(0, 4)], 7.0);
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct Array2<T>(ndarray::Array2<T>);

impl<T: Zero + Clone> Array2<T> {
    /// Create a new `Array2` of the specified `size` filled with the
    /// `Zero::zero` return value.
    ///
    /// # Examples
    ///
    /// ```
    /// # use lumol_core::types::Array2;
    /// let a: Array2<f64> = Array2::zeros((8, 5));
    /// assert_eq!(a[(6, 2)], 0.0);
    /// ```
    pub fn zeros(size: (usize, usize)) -> Array2<T> {
        Array2(ndarray::Array2::zeros(size))
    }

    /// Resize the array if the current size is not `size`, and fill the
    /// new array with zeros.
    ///
    /// # Examples
    ///
    /// ```
    /// # use lumol_core::types::Array2;
    /// let mut a = Array2::zeros((8, 5));
    ///
    /// a[(3, 3)] = 42.0;
    ///
    /// // This does nothing
    /// a.resize_if_different((8, 5));
    /// assert_eq!(a[(3, 3)], 42.0);
    ///
    /// // This allocates a new array
    /// a.resize_if_different((8, 9));
    /// assert_eq!(a[(3, 3)], 0.0);
    /// ```
    pub fn resize_if_different(&mut self, size: (usize, usize)) {
        if self.dim() != size {
            *self = Array2::zeros(size);
        }
    }
}

impl<T: Default> Array2<T> {
    /// Create a new `Array2` of the specified `size` filled with the
    /// `Default::default` return value.
    ///
    /// # Examples
    ///
    /// ```
    /// # use lumol_core::types::Array2;
    /// let a: Array2<f64> = Array2::zeros((8, 5));
    /// let b: Array2<f64> = Array2::default((8, 5));
    ///
    /// assert_eq!(a, b);
    /// ```
    pub fn default(size: (usize, usize)) -> Array2<T> {
        Array2(ndarray::Array2::default(size))
    }
}

impl<T> Index<(usize, usize)> for Array2<T> {
    type Output = T;
    fn index(&self, index: (usize, usize)) -> &T {
        unsafe {
            // ndarray does the check for us in debug builds
            self.0.uget(index)
        }
    }
}

impl<T> IndexMut<(usize, usize)> for Array2<T> {
    fn index_mut(&mut self, index: (usize, usize)) -> &mut T {
        unsafe {
            // ndarray does the check for us in debug builds
            self.0.uget_mut(index)
        }
    }
}

impl<T> Deref for Array2<T> {
    type Target = ndarray::Array2<T>;

    fn deref(&self) -> &ndarray::Array2<T> {
        &self.0
    }
}

impl<T> DerefMut for Array2<T> {
    fn deref_mut(&mut self) -> &mut ndarray::Array2<T> {
        &mut self.0
    }
}

/// Three dimensional tensors, based on ndarray
///
/// Most of the methods are simply forwarded to ndarray, so also look the doc
/// for this crate. This array type mainly supports indexing using tuples as
/// indices and is though as storage backend for multi-dimensional data.
///
/// ```
/// # use lumol_core::types::Array3;
/// let mut a = Array3::zeros((3, 5, 2));
///
/// assert_eq!(a[(0, 4, 1)], 0.0);
///
/// a[(0, 4, 1)] = 7.0;
/// assert_eq!(a[(0, 4, 1)], 7.0);
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct Array3<T>(ndarray::Array3<T>);

impl<T> Array3<T> {
    /// Create a new `Array3` of the specified `size` filled with the
    /// `Zero::zero` return value.
    ///
    /// # Examples
    ///
    /// ```
    /// # use lumol_core::types::Array3;
    /// let a: Array3<f64> = Array3::zeros((8, 5, 2));
    /// assert_eq!(a[(6, 2, 0)], 0.0);
    /// ```
    pub fn zeros(size: (usize, usize, usize)) -> Array3<T>
    where
        T: Zero + Clone,
    {
        Array3(ndarray::Array3::zeros(size))
    }

    /// Resize the array if the current size is not `size`, and fill the
    /// new array with zeros.
    ///
    /// # Examples
    ///
    /// ```
    /// # use lumol_core::types::Array3;
    /// let mut a = Array3::zeros((8, 5, 7));
    ///
    /// a[(3, 3, 3)] = 42.0;
    ///
    /// // This does nothing
    /// a.resize_if_different((8, 5, 7));
    /// assert_eq!(a[(3, 3, 3)], 42.0);
    ///
    /// // This allocates a new array
    /// a.resize_if_different((8, 5, 6));
    /// assert_eq!(a[(3, 3, 3)], 0.0);
    /// ```
    pub fn resize_if_different(&mut self, size: (usize, usize, usize))
    where
        T: Zero + Clone,
    {
        if self.0.shape() != [size.0, size.1, size.2] {
            *self = Array3::zeros(size);
        }
    }

    /// Create a new `Array3` of the specified `size` filled with the
    /// `Default::default` return value.
    /// `Default::default` return value.
    ///
    /// # Examples
    ///
    /// ```
    /// # use lumol_core::types::Array3;
    /// let a: Array3<f64> = Array3::zeros((8, 5, 2));
    /// let b: Array3<f64> = Array3::default((8, 5, 2));
    ///
    /// assert_eq!(a, b);
    /// ```
    pub fn default(size: (usize, usize, usize)) -> Array3<T>
    where
        T: Default,
    {
        Array3(ndarray::Array3::default(size))
    }
}

impl<T> Index<(usize, usize, usize)> for Array3<T> {
    type Output = T;
    fn index(&self, index: (usize, usize, usize)) -> &T {
        unsafe {
            // ndarray does the check for us in debug builds
            self.0.uget(index)
        }
    }
}

impl<T> IndexMut<(usize, usize, usize)> for Array3<T> {
    fn index_mut(&mut self, index: (usize, usize, usize)) -> &mut T {
        unsafe {
            // ndarray does the check for us in debug builds
            self.0.uget_mut(index)
        }
    }
}

impl<T> Deref for Array3<T> {
    type Target = ndarray::Array3<T>;

    fn deref(&self) -> &ndarray::Array3<T> {
        &self.0
    }
}

impl<T> DerefMut for Array3<T> {
    fn deref_mut(&mut self) -> &mut ndarray::Array3<T> {
        &mut self.0
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
            #[derive(Clone)]
            struct F64(f64);
            impl Default for F64 {
                fn default() -> F64 {
                    F64(42.0)
                }
            }

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
            assert_eq!(a.dim(), (3, 4));
            a[(1, 1)] = 42.0;

            a.resize_if_different((7, 90));
            assert_eq!(a.dim(), (7, 90));
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
            #[derive(Clone)]
            struct F64(f64);
            impl Default for F64 {
                fn default() -> F64 {
                    F64(42.0)
                }
            }

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
            assert_eq!(a.dim(), (3, 4, 5));
            a[(1, 1, 1)] = 42.0;

            a.resize_if_different((7, 90, 8));
            assert_eq!(a.dim(), (7, 90, 8));
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
