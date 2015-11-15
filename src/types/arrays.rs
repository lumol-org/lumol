/* Cymbalum, Molecular Simulation in Rust - Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
 */
//! Array types, for storage of data in multidimentional arrays. The types in
//! this module are only intented for data storage, and do not provide any
//! mathematical operations.
use std::ops::{Index, IndexMut};
use std::cmp;

/// Trait bound for data in the array types
pub trait Data: Clone + PartialEq + Default {}
impl<T: Clone + PartialEq + Default> Data for T {}

/******************************************************************************/
/// Bidimensional array, with cache locality and runtime size. This type is
/// intended for storage of data only, and do not provide any mathematical
/// operations. Data is stored in the row-major mode.
#[derive(Debug, Clone, PartialEq)]
pub struct Array2<T: Data> {
    /// Data storage
    data: Vec<T>,
    /// Size of the array
    size: (usize, usize)
}

impl<T: Data> Array2<T> {
    /// Get the size of the array
    pub fn size(&self) -> (usize, usize) {
        self.size
    }

    /// Create a new empty `Array2`.
    pub fn new() -> Array2<T> {
        Array2{
            data: Vec::new(),
            size: (0, 0)
        }
    }

    /// Create a new `Array2` of size `size` filed with `T::default()` values.
    pub fn with_size(size: (usize, usize)) -> Array2<T> {
        Array2::<T>::check_size(size);
        let data_size = size.0 * size.1;
        let vec = vec![T::default(); data_size];
        Array2{
            data: vec,
            size: size
        }
    }

    /// Get the index in the `data` vector corresponding to the point `i, j`
    #[inline] fn get_index(&self, i: usize, j: usize) -> usize {
        Array2::<T>::get_index_for_size(i, j, self.size)
    }

    /// Get the index in the `data` vector corresponding to the point `i, j`,
    /// and the size `size`.
    #[inline] fn get_index_for_size(i: usize, j: usize, size: (usize, usize)) -> usize {
        i * size.1 + j
    }

    fn check_size(size: (usize, usize)) {
        assert!(size.0 != 0, "size can not be null in dimmension 0");
        assert!(size.1 != 0, "size can not be null in dimmension 1");
    }

    /// Resize the array to be of size `new_size * new_size`. If the new size
    /// is bigger than the old size, the new values inserted are `T::default()`.
    /// Old values are preseved.
    pub fn resize(&mut self, size: (usize, usize)) {
        Array2::<T>::check_size(size);
        if size == self.size {return;}
        let mut new_data = vec![T::default(); size.0 * size.1];
        let n0 = cmp::min(self.size.0, size.0);
        let n1 = cmp::min(self.size.1, size.1);
        for i0 in 0..n0 {
            for i1 in 0..n1 {
                let new_idx = Array2::<T>::get_index_for_size(i0, i1, size);
                new_data[new_idx] = self[(i0, i1)].clone();
            }
        }
        self.data = new_data;
        self.size = size;
    }
}

impl<T: Data> Index<(usize, usize)> for Array2<T> {
    type Output = T;
    #[inline]
    fn index(&self, index: (usize, usize)) -> &T {
        debug_assert!(index.0 < self.size.0, format!("index out of bounds in dimmension 0: the len is {} but the index is {}", self.size.0, index.0));
        debug_assert!(index.1 < self.size.1, format!("index out of bounds in dimmension 1: the len is {} but the index is {}", self.size.1, index.1));
        let idx = self.get_index(index.0, index.1);
        unsafe {
            return self.data.get_unchecked(idx);
        }
    }
}

impl<T: Data> IndexMut<(usize, usize)> for Array2<T> {
    #[inline]
    fn index_mut(&mut self, index: (usize, usize)) -> &mut T {
        debug_assert!(index.0 < self.size.0, format!("index out of bounds in dimmension 0: the len is {} but the index is {}", self.size.0, index.0));
        debug_assert!(index.1 < self.size.1, format!("index out of bounds in dimmension 1: the len is {} but the index is {}", self.size.1, index.1));
        let idx = self.get_index(index.0, index.1);
        unsafe {
            return self.data.get_unchecked_mut(idx);
        }
    }
}

/******************************************************************************/
/// 3 dimmensional array type with cache locality, and runtime size. This type
/// is intended for storage of data only, and do not provide any mathematical
/// operations. Data is stored in the row-major mode.
#[derive(Debug, Clone, PartialEq)]
pub struct Array3<T: Data> {
    data: Vec<T>,
    size: (usize, usize, usize),
}

impl<T: Data> Array3<T> {
    /// Create a new empty tensor
    pub fn new() -> Array3<T> {
        Array3 {
            data: Vec::new(),
            size: (0, 0, 0),
        }
    }

    /// Get the current size of the array
    pub fn size(&self) -> (usize, usize, usize) {
        self.size
    }

    /// Create a new `Array3` of size `size` filed with `T::default()` values.
    pub fn with_size(size: (usize, usize, usize)) -> Array3<T> {
        Array3::<T>::check_size(size);
        let data_size = size.0 * size.1 * size.2;
        let vec = vec![T::default(); data_size];
        Array3 {
            data: vec,
            size: size,
        }
    }

    /// Get the index in the `data` vector corresponding to the point `i, j, k`
    #[inline] fn get_index(&self, i: usize, j: usize, k: usize) -> usize {
        Array3::<T>::get_index_for_size(i, j, k, self.size)
    }

    /// Get the index in the `data` vector corresponding to the point `i, j`,
    /// and the size `size`.
    #[inline] fn get_index_for_size(i: usize, j: usize, k:usize, size: (usize, usize, usize)) -> usize {
        i * size.1 * size.2 + j * size.2 + k
    }

    fn check_size(size: (usize, usize, usize)) {
        assert!(size.0 != 0, "size can not be null in dimmension 0");
        assert!(size.1 != 0, "size can not be null in dimmension 1");
        assert!(size.2 != 0, "size can not be null in dimmension 2");
    }

    /// Resize the array to be of size `new_size * new_size`. If the new size
    /// is bigger than the old size, the new values inserted are `T::default()`.
    /// Old values are preseved.
    pub fn resize(&mut self, size: (usize, usize, usize)) {
        Array3::<T>::check_size(size);
        if size == self.size {return;}
        let mut new_data = vec![T::default(); size.0 * size.1 * size.2];
        let n0 = cmp::min(self.size.0, size.0);
        let n1 = cmp::min(self.size.1, size.1);
        let n2 = cmp::min(self.size.2, size.2);
        for i0 in 0..n0 {
            for i1 in 0..n1 {
                for i2 in 0..n2 {
                    let new_idx = Array3::<T>::get_index_for_size(i0, i1, i2, size);
                    new_data[new_idx] = self[(i0, i1, i2)].clone();
                }
            }
        }
        self.data = new_data;
        self.size = size;
    }
}

impl<T: Data> Index<(usize, usize, usize)> for Array3<T> {
    type Output = T;
    #[inline]
    fn index(&self, index: (usize, usize, usize)) -> &T {
        debug_assert!(index.0 < self.size.0, format!("index out of bounds in dimmension 0: the len is {} but the index is {}", self.size.0, index.0));
        debug_assert!(index.1 < self.size.1, format!("index out of bounds in dimmension 1: the len is {} but the index is {}", self.size.1, index.1));
        debug_assert!(index.2 < self.size.2, format!("index out of bounds in dimmension 2: the len is {} but the index is {}", self.size.2, index.2));
        let idx = self.get_index(index.0, index.1, index.2);
        unsafe {
            return self.data.get_unchecked(idx);
        }
    }
}

impl<T: Data> IndexMut<(usize, usize, usize)> for Array3<T> {
    #[inline]
    fn index_mut(&mut self, index: (usize, usize, usize)) -> &mut T {
        debug_assert!(index.0 < self.size.0, format!("index out of bounds in dimmension 0: the len is {} but the index is {}", self.size.0, index.0));
        debug_assert!(index.1 < self.size.1, format!("index out of bounds in dimmension 1: the len is {} but the index is {}", self.size.1, index.1));
        debug_assert!(index.2 < self.size.2, format!("index out of bounds in dimmension 2: the len is {} but the index is {}", self.size.2, index.2));
        let idx = self.get_index(index.0, index.1, index.2);
        unsafe {
            return self.data.get_unchecked_mut(idx);
        }
    }
}

#[cfg(test)]
mod tests {
    pub use super::*;

    mod array_2 {
        use super::*;

        #[test]
        #[should_panic]
        fn out_of_bounds_dim_1() {
            let a = Array2::<f64>::with_size((4, 5));
            let _ = a[(6, 1)];
        }

        #[test]
        #[should_panic]
        fn out_of_bounds_dim_2() {
            let a = Array2::<f64>::with_size((4, 5));
            let _ = a[(1, 6)];
        }

        #[test]
        fn size() {
            let a = Array2::<f64>::with_size((6, 4));
            assert_eq!(a.size(), (6, 4));
        }

        #[test]
        fn index() {
            let mut a = Array2::<f64>::with_size((4, 8));
            for i in 0..4 {
                for j in 0..8 {
                    assert_eq!(a[(i, j)], 0.0);
                }
            }

            a[(2, 3)] = 7.0;
            assert_eq!(a[(2, 3)], 7.0);
        }

        #[test]
        fn resize() {
            let mut a = Array2::<f64>::with_size((6, 5));
            assert_eq!(a.size(), (6, 5));
            a[(2, 3)] = 42.0;

            a.resize((9, 7));
            assert_eq!(a.size(), (9, 7));
            assert_eq!(a[(2, 3)], 42.0);

            a.resize((3, 4));
            assert_eq!(a.size(), (3, 4));
            assert_eq!(a[(2, 3)], 42.0);
        }

        #[test]
        #[should_panic]
        fn resize_to_zero() {
            let mut a = Array2::<f64>::with_size((6, 5));
            a.resize((0, 3));
        }
    }

    mod array_3 {
        use super::*;

        #[test]
        #[should_panic]
        fn out_of_bounds_dim_1() {
            let a = Array3::<f64>::with_size((4, 5, 6));
            let _ = a[(6, 1, 3)];
        }

        #[test]
        #[should_panic]
        fn out_of_bounds_dim_2() {
            let a = Array3::<f64>::with_size((4, 5, 6));
            let _ = a[(1, 7, 3)];
        }

        #[test]
        #[should_panic]
        fn out_of_bounds_dim_3() {
            let a = Array3::<f64>::with_size((4, 5, 6));
            let _ = a[(0, 1, 8)];
        }

        #[test]
        fn size() {
            let a = Array3::<f64>::with_size((6, 4, 9));
            assert_eq!(a.size(), (6, 4, 9));
        }

        #[test]
        fn index() {
            let mut a = Array3::<f64>::with_size((4, 5, 3));
            for i in 0..4 {
                for j in 0..5 {
                    for k in 0..3 {
                        assert_eq!(a[(i, j, k)], 0.0);
                    }
                }
            }

            a[(2, 3, 0)] = 7.0;
            assert_eq!(a[(2, 3, 0)], 7.0);
        }

        #[test]
        fn resize() {
            let mut a = Array3::<f64>::with_size((6, 6, 8));
            assert_eq!(a.size(), (6, 6, 8));
            a[(2, 3, 0)] = 42.0;

            a.resize((9, 10, 2));
            assert_eq!(a.size(), (9, 10, 2));
            assert_eq!(a[(2, 3, 0)], 42.0);

            a.resize((3, 4, 1));
            assert_eq!(a.size(), (3, 4, 1));
            assert_eq!(a[(2, 3, 0)], 42.0);
        }

        #[test]
        #[should_panic]
        fn resize_to_zero() {
            let mut a = Array3::<f64>::with_size((6, 5, 7));
            a.resize((0, 3, 2));
        }
    }
}
