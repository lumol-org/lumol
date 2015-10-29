/* Cymbalum, Molecular Simulation in Rust - Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
 */

//! Matrix types: generic square matrix, and 3x3 matrix.
use std::slice;
use std::ops::{Add, Sub, Mul, Index, IndexMut};
use super::vectors::Vector3D;

/// 3x3 dimensional matrix type, implementing all usual operations
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Matrix3 {
    data : [[f64; 3]; 3]
}

impl Matrix3 {
    /// Create a new `Matrix3` specifying all its components
    pub fn new(m00: f64, m01: f64, m02: f64,
               m10: f64, m11: f64, m12: f64,
               m20: f64, m21: f64, m22: f64) -> Matrix3 {
        Matrix3{data: [[m00, m01, m02],
                       [m10, m11, m12],
                       [m20, m21, m22]]}
    }

    /// Create an identity matrix
    pub fn one() -> Matrix3 {
        Matrix3::new(1.0, 0.0, 0.0,
                     0.0, 1.0, 0.0,
                     0.0, 0.0, 1.0)
    }

    /// Create a new null matrix
    pub fn zero() -> Matrix3 {
        Matrix3::new(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    }

    /// Computes the inverse of a matrix, which is assumed to exist
    pub fn inverse(&self) -> Matrix3 {
        let mut determinant = 0.0;
        determinant += self[(0, 0)] * (self[(1, 1)] * self[(2, 2)] - self[(2, 1)] * self[(1, 2)]);
        determinant -= self[(0, 1)] * (self[(1, 0)] * self[(2, 2)] - self[(1, 2)] * self[(2, 0)]);
        determinant += self[(0, 2)] * (self[(1, 0)] * self[(2, 1)] - self[(1, 1)] * self[(2, 0)]);;

        assert!(determinant != 0.0, "The matrix is not inversible!");
        let invdet = 1.0 / determinant;
        let mut res = Matrix3::zero();
        res[(0, 0)] = (self[(1, 1)] * self[(2, 2)] - self[(2, 1)] * self[(1, 2)]) * invdet;
        res[(0, 1)] = (self[(0, 2)] * self[(2, 1)] - self[(0, 1)] * self[(2, 2)]) * invdet;
        res[(0, 2)] = (self[(0, 1)] * self[(1, 2)] - self[(0, 2)] * self[(1, 1)]) * invdet;
        res[(1, 0)] = (self[(1, 2)] * self[(2, 0)] - self[(1, 0)] * self[(2, 2)]) * invdet;
        res[(1, 1)] = (self[(0, 0)] * self[(2, 2)] - self[(0, 2)] * self[(2, 0)]) * invdet;
        res[(1, 2)] = (self[(1, 0)] * self[(0, 2)] - self[(0, 0)] * self[(1, 2)]) * invdet;
        res[(2, 0)] = (self[(1, 0)] * self[(2, 1)] - self[(2, 0)] * self[(1, 1)]) * invdet;
        res[(2, 1)] = (self[(2, 0)] * self[(0, 1)] - self[(0, 0)] * self[(2, 1)]) * invdet;
        res[(2, 2)] = (self[(0, 0)] * self[(1, 1)] - self[(1, 0)] * self[(0, 1)]) * invdet;
        return res;
    }
}

impl Index<usize> for Matrix3 {
    type Output = [f64; 3];
    #[inline]
    fn index(&self, index: usize) -> &[f64; 3] {
        &self.data[index]
    }
}

impl Index<(usize, usize)> for Matrix3 {
    type Output = f64;
    #[inline]
    fn index(&self, index: (usize, usize)) -> &f64 {
        &self.data[index.0][index.1]
    }
}

impl IndexMut<usize> for Matrix3 {
    #[inline]
    fn index_mut(&mut self, index: usize) -> &mut [f64; 3] {
        &mut self.data[index]
    }
}

impl IndexMut<(usize, usize)> for Matrix3 {
    #[inline]
    fn index_mut(&mut self, index: (usize, usize)) -> &mut f64 {
        &mut self.data[index.0][index.1]
    }
}

/// Add two matrix
impl Add for Matrix3 {
    type Output = Matrix3;
    fn add(self, other: Matrix3) -> Matrix3 {
        Matrix3::new(self[0][0] + other[0][0], self[0][1] + other[0][1], self[0][2] + other[0][2],
                     self[1][0] + other[1][0], self[1][1] + other[1][1], self[1][2] + other[1][2],
                     self[2][0] + other[2][0], self[2][1] + other[2][1], self[2][2] + other[2][2])
    }
}

/// Substract two matrix
impl Sub for Matrix3 {
    type Output = Matrix3;
    fn sub(self, other: Matrix3) -> Matrix3 {
        Matrix3::new(self[0][0] - other[0][0], self[0][1] - other[0][1], self[0][2] - other[0][2],
                     self[1][0] - other[1][0], self[1][1] - other[1][1], self[1][2] - other[1][2],
                     self[2][0] - other[2][0], self[2][1] - other[2][1], self[2][2] - other[2][2])
    }
}

/// Multiply by a scalar on the right hand side
impl Mul<f64> for Matrix3 {
    type Output = Matrix3;
    fn mul(self, other: f64) -> Matrix3 {
        Matrix3::new(self[0][0] * other, self[0][1] * other, self[0][2] * other,
                     self[1][0] * other, self[1][1] * other, self[1][2] * other,
                     self[2][0] * other, self[2][1] * other, self[2][2] * other)
    }
}

/// Multiply by a scalar on the left hand side
impl Mul<Matrix3> for f64 {
    type Output = Matrix3;
    fn mul(self, other: Matrix3) -> Matrix3 {
        Matrix3::new(self * other[0][0], self * other[0][1], self * other[0][2],
                     self * other[1][0], self * other[1][1], self * other[1][2],
                     self * other[2][0], self * other[2][1], self * other[2][2])
    }
}

/// Product of the two matrix
impl Mul<Matrix3> for Matrix3 {
    type Output = Matrix3;
    fn mul(self, other: Matrix3) -> Matrix3 {
        let m00 = self[(0, 0)] * other[(0, 0)] + self[(0, 1)] * other[(1, 0)] + self[(0, 2)] * other[(2, 0)];
        let m01 = self[(0, 0)] * other[(0, 1)] + self[(0, 1)] * other[(1, 1)] + self[(0, 2)] * other[(2, 1)];
        let m02 = self[(0, 0)] * other[(0, 2)] + self[(0, 1)] * other[(1, 2)] + self[(0, 2)] * other[(2, 2)];

        let m10 = self[(1, 0)] * other[(0, 0)] + self[(1, 1)] * other[(1, 0)] + self[(1, 2)] * other[(2, 0)];
        let m11 = self[(1, 0)] * other[(0, 1)] + self[(1, 1)] * other[(1, 1)] + self[(1, 2)] * other[(2, 1)];
        let m12 = self[(1, 0)] * other[(0, 2)] + self[(1, 1)] * other[(1, 2)] + self[(1, 2)] * other[(2, 2)];

        let m20 = self[(2, 0)] * other[(0, 0)] + self[(2, 1)] * other[(1, 0)] + self[(2, 2)] * other[(2, 0)];
        let m21 = self[(2, 0)] * other[(0, 1)] + self[(2, 1)] * other[(1, 1)] + self[(2, 2)] * other[(2, 1)];
        let m22 = self[(2, 0)] * other[(0, 2)] + self[(2, 1)] * other[(1, 2)] + self[(2, 2)] * other[(2, 2)];

        Matrix3::new(m00, m01, m02, m10, m11, m12, m20, m21, m22)
    }
}

/// Multiply by a Vector3D
impl Mul<Vector3D> for Matrix3 {
    type Output = Vector3D;
    fn mul(self, vec: Vector3D) -> Vector3D {
        let x = self[0][0] * vec.x + self[0][1] * vec.y + self[0][2] * vec.z;
        let y = self[1][0] * vec.x + self[1][1] * vec.y + self[1][2] * vec.z;
        let z = self[2][0] * vec.x + self[2][1] * vec.y + self[2][2] * vec.z;
        Vector3D::new(x, y, z)
    }
}

/******************************************************************************/
/// Bidimensional square matrix, with cache locality and runtime size. This type
/// only provide an indexable bidimensional storage.
#[derive(Debug)]
pub struct Matrix<T> {
    /// Data storage
    data: Vec<T>,
    /// Size of the matrix (it is an NxN matrix)
    n: usize
}

impl<T> Matrix<T> {
    /// Get the size of the matrix
    pub fn size(&self) -> usize {
        self.n
    }
}

impl<T: Clone> Matrix<T> {
    /// Create a new `Matrix` of size `n x n` filed with `value`.
    pub fn new(value: T, n: usize) -> Matrix<T> {
        let mut vec = Vec::with_capacity(n*n);
        for _ in 0..n*n {
            vec.push(value.clone())
        }
        Matrix{data: vec, n: n}
    }
}

impl<T: Default + Clone> Matrix<T> {
    /// Create a new `Matrix` of size `n x n` filed with the default value for
    /// `T`.
    pub fn with_size(n: usize) -> Matrix<T> {
        let vec = vec![T::default(); n*n];
        Matrix{data: vec, n: n}
    }

    /// Resize the matrix to be of size `new_size * new_size`. If the new size
    /// is bigger than the old size, the new values inserted are T::default().
    /// Old values are preseved.
    pub fn resize(&mut self, new_size: usize) {
        if new_size == self.n { return;}
        if new_size > self.n {
            let mut new_data = vec![T::default(); new_size * new_size];
            for i in 1..self.n {
                for j in 1..self.n {
                    new_data[i*new_size + j] = self[i][j].clone();
                }
            }
            self.data = new_data;
        } else {
            for i in 1..new_size {
                for j in 1..new_size {
                    self.data[i*new_size + j] = self.data[i*self.n + j].clone();
                }
            }
            self.data.truncate(new_size * new_size);
        }
        self.n = new_size;
    }
}

impl<T: PartialEq> PartialEq for Matrix<T> {
    fn eq(&self, rhs: &Matrix<T>) -> bool {
        if self.n != rhs.n {
            return false;
        }
        return self.data == rhs.data;
    }
}

impl<T: Clone> Clone for Matrix<T> {
    fn clone(&self) -> Matrix<T> {
        Matrix{
            data: self.data.clone(),
            n: self.n
        }
    }
}

impl<T> Index<usize> for Matrix<T> {
    type Output = [T];
    fn index<'a>(&'a self, i: usize) -> &'a [T] {
        assert!(i < self.n, format!("index out of bounds: the len is {} but the index is {}", self.n, i));
        unsafe {
            let ptr = self.data.as_ptr().offset((i * self.n) as isize);
            slice::from_raw_parts(ptr,self.n)
        }
    }
}

impl<T> IndexMut<usize> for Matrix<T> {
    fn index_mut<'a>(&'a mut self, i: usize) -> &'a mut [T] {
        assert!(i < self.n, format!("index out of bounds: the len is {} but the index is {}", self.n, i));
        unsafe {
            let ptr = self.data.as_mut_ptr().offset((i * self.n) as isize);
            slice::from_raw_parts_mut(ptr, self.n)
        }
    }
}

/******************************************************************************/

#[cfg(test)]
mod tests {
    pub use super::*;

    mod matrix {
        use super::*;

        #[test]
        #[should_panic]
        fn out_of_bounds() {
            let a = Matrix::<f64>::with_size(4);
            let _ = a[6][1];
        }

        #[test]
        fn size() {
            let mut a = Matrix::<f64>::with_size(4);
            assert_eq!(a.size(), 4);

            for i in 0..4 {
                for j in 0..4 {
                    assert_eq!(a[i][j], 0.0);
                }
            }

            a[2][3] = 7.0;
            assert_eq!(a[2][3], 7.0);

            a = Matrix::<f64>::new(42.0, 8);
            assert_eq!(a.size(), 8);
            for i in 0..8 {
                for j in 0..8 {
                    assert_eq!(a[i][j], 42.0);
                }
            }
        }

        #[test]
        fn resize() {
            let mut a = Matrix::<f64>::with_size(6);
            assert_eq!(a.size(), 6);
            a[2][3] = 42.0;

            a.resize(9);
            assert_eq!(a.size(), 9);
            assert_eq!(a[2][3], 42.0);

            a.resize(4);
            assert_eq!(a[2][3], 42.0);
        }
    }

    mod matrix3 {
        use super::*;
        use types::Vector3D;

        #[test]
        #[should_panic]
        fn out_of_bounds() {
            let a = Matrix3::zero();
            let _ = a[(3, 1)];
        }

        #[test]
        fn specials_matrix() {
            let a = Matrix3::zero();
            let b = Matrix3::one();

            for i in (0..2) {
                for j in (0..2) {
                    assert_eq!(a[(i, j)], 0.0);
                    if i == j {
                        assert_eq!(b[(i, j)], 1.0);
                    } else {
                        assert_eq!(b[(i, j)], 0.0);
                    }
                }
            }
        }


        #[test]
        fn add() {
            let a = Matrix3::new(1.0, 2.0, 3.0,
                                 4.0, 5.0, 6.0,
                                 8.0, 9.0, 10.0);
            let b = Matrix3::one();
            let add = a + b;

            let res = Matrix3::new(2.0, 2.0, 3.0,
                                   4.0, 6.0, 6.0,
                                   8.0, 9.0, 11.0);

            for i in (0..3) {
                for j in (0..3) {
                    assert_eq!(add[(i, j)], res[(i, j)]);
                }
            }
        }

        #[test]
        fn sub() {
            let a = Matrix3::new(1.0, 2.0, 3.0,
                                 4.0, 5.0, 6.0,
                                 8.0, 9.0, 10.0);
            let b = Matrix3::one();
            let sub = a - b;

            let res = Matrix3::new(0.0, 2.0, 3.0,
                                   4.0, 4.0, 6.0,
                                   8.0, 9.0, 9.0);

            for i in (0..3) {
                for j in (0..3) {
                    assert_eq!(sub[(i, j)], res[(i, j)]);
                }
            }
        }

        #[test]
        fn mul_scalar() {
            let a = Matrix3::new(1.0, 2.0, 3.0,
                                 4.0, 5.0, 6.0,
                                 8.0, 9.0, 10.0);
            let b = 2.0;
            let mul_r = a * b;
            let mul_l = b * a;

            let res = Matrix3::new(2.0, 4.0, 6.0,
                                   8.0, 10.0, 12.0,
                                   16.0, 18.0, 20.0);

            for i in (0..3) {
                for j in (0..3) {
                    assert_eq!(mul_r[(i, j)], res[(i, j)]);
                    assert_eq!(mul_l[(i, j)], res[(i, j)]);
                }
            }
        }

        #[test]
        fn mul_matrix() {
            let unit = Matrix3::one();
            let unit_sq = unit * unit;

            let a = Matrix3::new(2.0, 4.0, 6.0,
                                 8.0, 10.0, 12.0,
                                 16.0, 18.0, 20.0);
            let mul_r = unit * a;
            let mul_l = a * unit;

            for i in 0..3 {
                for j in 0..3 {
                    // 1-Matrix is the same as (1-Matrix)*(1-Matrix)
                    assert_eq!(unit[(i, j)], unit_sq[(i, j)]);
                    // (1-Matrix)*A is the same as A and A*(1-Matrix)
                    assert_eq!(mul_l[(i, j)], a[(i, j)]);
                    assert_eq!(mul_r[(i, j)], a[(i, j)]);
                }
            }
        }

        #[test]
        fn mul_vector() {
            let A = Matrix3::new(1.0, 2.0, 3.0,
                                 4.0, 5.0, 6.0,
                                 8.0, 9.0, 10.0);

            let vec = Vector3D::new(1.0, 1.0, 1.0);
            let mul = A * vec;
            let res = Vector3D::new(6.0, 15.0, 27.0);

            assert_eq!(mul.x, res.x);
            assert_eq!(mul.y, res.y);
            assert_eq!(mul.z, res.z);

            let unit = Matrix3::one();
            let vec = Vector3D::new(567.45, 356.8, 215673.12);
            let mul = unit * vec;
            let res = vec;

            assert_eq!(mul.x, res.x);
            assert_eq!(mul.y, res.y);
            assert_eq!(mul.z, res.z);
        }

        #[test]
        fn inverse() {
            let I = Matrix3::one();
            let res = I.inverse();

            let A = Matrix3::new(1.0, 2.0, 3.0,
                                 2.0, 5.0, 3.0,
                                 1.0, 3.0, 8.0,);
            let B = A.inverse();
            let C = A*B;

            for i in 0..3 {
                for j in 0..3 {
                    assert_eq!(I[(i, j)], res[(i, j)]);
                    assert_eq!(I[(i, j)], C[(i, j)]);
                }
            }
        }
    }
}
