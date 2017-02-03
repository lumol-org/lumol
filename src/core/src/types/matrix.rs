// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 Lumol's contributors — BSD license

//! 3x3 matrix type.
use std::ops::{Add, Sub, Mul, Div, Index, IndexMut};
use std::ops::{AddAssign, SubAssign, MulAssign, DivAssign};
use types::{Vector3D, Zero, One};

/// A 3x3 square matrix type.
///
/// `Matrix3` implements all the usual arithmetic operations:
///
/// ```
/// use lumol::types::{Matrix3, Vector3D, One};
///
/// let one = Matrix3::one();
/// let a = Matrix3::new(
///     1.0, 0.0, 3.0,
///     0.0, 2.0, 5.6,
///     0.0, 0.0, 8.0
/// );
///
/// let v = Vector3D::new(3.0, 2.0, 1.0);
///
/// // Indexing
/// assert_eq!(a[0][0], 1.0);
/// assert_eq!(a[1][2], 5.6);
///
/// // Addition
/// let c = a + one;
/// assert_eq!(c, Matrix3::new(
///     2.0, 0.0, 3.0,
///     0.0, 3.0, 5.6,
///     0.0, 0.0, 9.0
/// ));
///
/// // Subtraction
/// let c = a - one;
/// assert_eq!(c, Matrix3::new(
///     0.0, 0.0, 3.0,
///     0.0, 1.0, 5.6,
///     0.0, 0.0, 7.0
/// ));
///
/// // Multiplication
/// let c = a * one;  // matrix - matrix
/// assert_eq!(c, a);
///
/// let c = a * v;  // matrix - vector
/// assert_eq!(c, Vector3D::new(6.0, 9.6, 8.0));
///
/// let c = 42.0 * one;  // matrix - scalar
/// assert_eq!(c, Matrix3::new(
///     42., 0.0, 0.0,
///     0.0, 42., 0.0,
///     0.0, 0.0, 42.
/// ));
///
/// // Division
/// let c = a / 2.0;
/// assert_eq!(c, Matrix3::new(
///     0.5, 0.0, 1.5,
///     0.0, 1.0, 2.8,
///     0.0, 0.0, 4.0
/// ));
/// ```
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Matrix3 {
    data : [[f64; 3]; 3]
}

impl Matrix3 {
    /// Create a new `Matrix3` specifying all its components
    /// # Examples
    /// ```
    /// # use lumol::types::Matrix3;
    /// let matrix = Matrix3::new(
    ///     0.0, 0.0, 3.0,
    ///     0.0, 1.0, 5.6,
    ///     0.0, 0.0, 7.0
    /// );
    /// assert_eq!(matrix[0][2], 3.0);
    /// ```
    #[allow(too_many_arguments)]
    pub fn new(m00: f64, m01: f64, m02: f64,
               m10: f64, m11: f64, m12: f64,
               m20: f64, m21: f64, m22: f64) -> Matrix3 {
        Matrix3{data: [[m00, m01, m02],
                       [m10, m11, m12],
                       [m20, m21, m22]]}
    }

    /// Compute the trace of the matrix
    /// # Examples
    /// ```
    /// # use lumol::types::Matrix3;
    /// let matrix = Matrix3::new(
    ///     0.0, 0.0, 3.0,
    ///     0.0, 1.0, 5.6,
    ///     0.0, 0.0, 7.0
    /// );
    /// assert_eq!(matrix.trace(), 8.0);
    /// ```
    pub fn trace(&self) -> f64 {
        return self[(0, 0)] + self[(1, 1)] + self[(2, 2)];
    }

    /// Computes the inverse of a matrix
    ///
    /// # Examples
    ///
    /// ```
    /// # use lumol::types::{Matrix3, One};
    /// // A diagonal matrix is trivially invertible
    /// let matrix = Matrix3::new(
    ///     4.0, 0.0, 0.0,
    ///     0.0, 1.0, 0.0,
    ///     0.0, 0.0, 7.0
    /// );
    ///
    /// let inverted = Matrix3::new(
    ///     1.0 / 4.0, 0.0,    0.0,
    ///        0.0,    1.0,    0.0,
    ///        0.0,    0.0, 1.0 / 7.0
    /// );
    ///
    /// assert_eq!(matrix.inverse(), inverted);
    /// assert_eq!(matrix * matrix.inverse(), Matrix3::one());
    /// ```
    ///
    /// # Panics
    ///
    /// If the matrix is not invertible, *i.e.* if the matrix determinant
    /// equals zero.
    pub fn inverse(&self) -> Matrix3 {
        let determinant = self.determinant();
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

    /// Computes the [determinant][Wiki] of a matrix
    ///
    /// [Wiki]: https://en.wikipedia.org/wiki/Determinant
    ///
    /// # Examples
    ///
    /// ```
    /// # use lumol::types::Matrix3;
    /// let matrix = Matrix3::new(
    ///     4.0, 0.0, 0.0,
    ///     0.0, 1.5, 0.0,
    ///     0.0, 0.0, 7.0
    /// );
    ///
    /// assert_eq!(matrix.determinant(), 4.0 * 1.5 * 7.0);
    /// ```
    pub fn determinant(&self) -> f64 {
        ( self[(0, 0)] * (self[(1, 1)] * self[(2, 2)] - self[(2, 1)] * self[(1, 2)])
        - self[(0, 1)] * (self[(1, 0)] * self[(2, 2)] - self[(1, 2)] * self[(2, 0)])
        + self[(0, 2)] * (self[(1, 0)] * self[(2, 1)] - self[(1, 1)] * self[(2, 0)]))
    }

    /// Transpose this matrix into a new matrix
    /// # Examples
    ///
    /// ```
    /// # use lumol::types::Matrix3;
    /// let matrix = Matrix3::new(
    ///     1.0, 2.0, 4.0,
    ///     0.0, 1.0, 3.0,
    ///     0.0, 0.0, 1.0
    /// );
    ///
    /// let transposed = Matrix3::new(
    ///     1.0, 0.0, 0.0,
    ///     2.0, 1.0, 0.0,
    ///     4.0, 3.0, 1.0
    /// );
    ///
    /// assert_eq!(matrix.transposed(), transposed);
    /// ```
    pub fn transposed(&self) -> Matrix3 {
        Matrix3::new(
            self[(0, 0)], self[(1, 0)], self[(2, 0)],
            self[(0, 1)], self[(1, 1)], self[(2, 1)],
            self[(0, 2)], self[(1, 2)], self[(2, 2)]
        )
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

impl_arithmetic!(
    Matrix3, Matrix3, Add, add, Matrix3,
    self, other,
    Matrix3::new(self[0][0] + other[0][0], self[0][1] + other[0][1], self[0][2] + other[0][2],
                 self[1][0] + other[1][0], self[1][1] + other[1][1], self[1][2] + other[1][2],
                 self[2][0] + other[2][0], self[2][1] + other[2][1], self[2][2] + other[2][2])
);

impl_inplace_arithmetic!(
    Matrix3, Matrix3, AddAssign, add_assign,
    self, other,
    {
        self[0][0] += other[0][0]; self[0][1] += other[0][1]; self[0][2] += other[0][2];
        self[1][0] += other[1][0]; self[1][1] += other[1][1]; self[1][2] += other[1][2];
        self[2][0] += other[2][0]; self[2][1] += other[2][1]; self[2][2] += other[2][2];
    }
);


impl_arithmetic!(
    Matrix3, Matrix3, Sub, sub, Matrix3,
    self, other,
    Matrix3::new(self[0][0] - other[0][0], self[0][1] - other[0][1], self[0][2] - other[0][2],
                 self[1][0] - other[1][0], self[1][1] - other[1][1], self[1][2] - other[1][2],
                 self[2][0] - other[2][0], self[2][1] - other[2][1], self[2][2] - other[2][2])
);

impl_inplace_arithmetic!(
    Matrix3, Matrix3, SubAssign, sub_assign,
    self, other,
    {
        self[0][0] -= other[0][0]; self[0][1] -= other[0][1]; self[0][2] -= other[0][2];
        self[1][0] -= other[1][0]; self[1][1] -= other[1][1]; self[1][2] -= other[1][2];
        self[2][0] -= other[2][0]; self[2][1] -= other[2][1]; self[2][2] -= other[2][2];
    }
);

impl_arithmetic!(
    Matrix3, Matrix3, Mul, mul, Matrix3,
    self, other,
    {
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
);

impl_inplace_arithmetic!(
    Matrix3, Matrix3, MulAssign, mul_assign,
    self, other,
    {
        let m00 = self[(0, 0)] * other[(0, 0)] + self[(0, 1)] * other[(1, 0)] + self[(0, 2)] * other[(2, 0)];
        let m01 = self[(0, 0)] * other[(0, 1)] + self[(0, 1)] * other[(1, 1)] + self[(0, 2)] * other[(2, 1)];
        let m02 = self[(0, 0)] * other[(0, 2)] + self[(0, 1)] * other[(1, 2)] + self[(0, 2)] * other[(2, 2)];

        let m10 = self[(1, 0)] * other[(0, 0)] + self[(1, 1)] * other[(1, 0)] + self[(1, 2)] * other[(2, 0)];
        let m11 = self[(1, 0)] * other[(0, 1)] + self[(1, 1)] * other[(1, 1)] + self[(1, 2)] * other[(2, 1)];
        let m12 = self[(1, 0)] * other[(0, 2)] + self[(1, 1)] * other[(1, 2)] + self[(1, 2)] * other[(2, 2)];

        let m20 = self[(2, 0)] * other[(0, 0)] + self[(2, 1)] * other[(1, 0)] + self[(2, 2)] * other[(2, 0)];
        let m21 = self[(2, 0)] * other[(0, 1)] + self[(2, 1)] * other[(1, 1)] + self[(2, 2)] * other[(2, 1)];
        let m22 = self[(2, 0)] * other[(0, 2)] + self[(2, 1)] * other[(1, 2)] + self[(2, 2)] * other[(2, 2)];

        *self = Matrix3::new(m00, m01, m02, m10, m11, m12, m20, m21, m22)
    }
);

impl_arithmetic!(
    Matrix3, Vector3D, Mul, mul, Vector3D,
    self, other,
    {
        let x = self[0][0] * other[0] + self[0][1] * other[1] + self[0][2] * other[2];
        let y = self[1][0] * other[0] + self[1][1] * other[1] + self[1][2] * other[2];
        let z = self[2][0] * other[0] + self[2][1] * other[1] + self[2][2] * other[2];
        Vector3D::new(x, y, z)
    }
);

/******************************************************************************/

lsh_scal_arithmetic!(
    Matrix3, Mul, mul, Matrix3,
    self, other,
    Matrix3::new(self[0][0] * other, self[0][1] * other, self[0][2] * other,
                 self[1][0] * other, self[1][1] * other, self[1][2] * other,
                 self[2][0] * other, self[2][1] * other, self[2][2] * other)
);

rhs_scal_arithmetic!(
    Matrix3, Mul, mul, Matrix3,
    self, other,
    Matrix3::new(self * other[0][0], self * other[0][1], self * other[0][2],
                 self * other[1][0], self * other[1][1], self * other[1][2],
                 self * other[2][0], self * other[2][1], self * other[2][2])
);

impl_inplace_arithmetic!(
    Matrix3, f64, MulAssign, mul_assign,
    self, other,
    {
        let other = other.clone();
        self[0][0] *= other; self[0][1] *= other; self[0][2] *= other;
        self[1][0] *= other; self[1][1] *= other; self[1][2] *= other;
        self[2][0] *= other; self[2][1] *= other; self[2][2] *= other;
    }
);

lsh_scal_arithmetic!(
    Matrix3, Div, div, Matrix3,
    self, other,
    Matrix3::new(self[0][0] / other, self[0][1] / other, self[0][2] / other,
                 self[1][0] / other, self[1][1] / other, self[1][2] / other,
                 self[2][0] / other, self[2][1] / other, self[2][2] / other)
);

impl_inplace_arithmetic!(
    Matrix3, f64, DivAssign, div_assign,
    self, other,
    {
        let other = other.clone();
        self[0][0] /= other; self[0][1] /= other; self[0][2] /= other;
        self[1][0] /= other; self[1][1] /= other; self[1][2] /= other;
        self[2][0] /= other; self[2][1] /= other; self[2][2] /= other;
    }
);

impl Zero for Matrix3 {
    fn zero() -> Matrix3 {
        Matrix3::new(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    }

    fn is_zero(&self) -> bool {
        *self == Matrix3::zero()
    }
}

impl One for Matrix3 {
    /// Create an identity matrix
    fn one() -> Matrix3 {
        Matrix3::new(1.0, 0.0, 0.0,
                     0.0, 1.0, 0.0,
                     0.0, 0.0, 1.0)
    }
}

#[cfg(test)]
mod tests {
    // TODO: remove most of these when we can generate coverage from doc tests
    // see https://github.com/rust-lang/rust/issues/36956
    use super::*;
    use types::{Vector3D, Zero, One};

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

        for i in 0..3 {
            for j in 0..3 {
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
        let mut a = Matrix3::new(1.0, 2.0, 3.0,
                                 4.0, 5.0, 6.0,
                                 8.0, 9.0, 10.0);
        let res = Matrix3::new(2.0, 2.0, 3.0,
                               4.0, 6.0, 6.0,
                               8.0, 9.0, 11.0);

        assert_eq!(a + Matrix3::one(), res);

        a += Matrix3::one();
        assert_eq!(a, res);
    }

    #[test]
    fn sub() {
        let mut a = Matrix3::new(1.0, 2.0, 3.0,
                                 4.0, 5.0, 6.0,
                                 8.0, 9.0, 10.0);
        let res = Matrix3::new(0.0, 2.0, 3.0,
                               4.0, 4.0, 6.0,
                               8.0, 9.0, 9.0);

        assert_eq!(a - Matrix3::one(), res);

        a -= Matrix3::one();
        assert_eq!(a, res);
    }

    #[test]
    fn mul_scalar() {
        let mut a = Matrix3::new(1.0, 2.0, 3.0,
                                 4.0, 5.0, 6.0,
                                 8.0, 9.0, 10.0);
        let res = Matrix3::new(2.0, 4.0, 6.0,
                               8.0, 10.0, 12.0,
                               16.0, 18.0, 20.0);
        assert_eq!(a * 2.0, res);
        assert_eq!(2.0 * a, res);

        a *= 2.0;
        assert_eq!(a, res);
    }

    #[test]
    fn div_scalar() {
        let mut a = Matrix3::new(1.0, 2.0, 3.0,
                                 4.0, 5.0, 6.0,
                                 8.0, 9.0, 10.0);
        let res = Matrix3::new(0.5, 1.0, 1.5,
                               2.0, 2.5, 3.0,
                               4.0, 4.5, 5.0);
        assert_eq!(a / 2.0, res);

        a /= 2.0;
        assert_eq!(a, res);
    }

    #[test]
    fn mul_matrix() {
        let unit = Matrix3::one();
        let mut a = Matrix3::new(2.0, 4.0, 6.0,
                                 8.0, 10.0, 12.0,
                                 16.0, 18.0, 20.0);

        assert_eq!(unit, unit * unit);
        assert_eq!(unit * a, a);
        assert_eq!(a * unit, a);

        let previous = a;
        a *= unit;
        assert_eq!(previous, a);
    }

    #[test]
    fn mul_vector() {
        let a = Matrix3::new(1.0, 2.0, 3.0,
                             4.0, 5.0, 6.0,
                             8.0, 9.0, 10.0);

        let vec = Vector3D::new(1.0, 1.0, 1.0);
        assert_eq!(a * vec, Vector3D::new(6.0, 15.0, 27.0));

        let unit = Matrix3::one();
        let vec = Vector3D::new(567.45, 356.8, 215673.12);
        assert_eq!(unit * vec, vec);
    }

    #[test]
    fn inverse() {
        let one = Matrix3::one();
        assert_eq!(one, one.inverse());

        let a = Matrix3::new(1.0, 2.0, 3.0,
                             2.0, 5.0, 3.0,
                             1.0, 3.0, 8.0,);
        assert!(a.determinant() != 0.0);
        assert_eq!(one, a * a.inverse());
    }

    #[test]
    fn determinant() {
        let one = Matrix3::one();
        assert_eq!(one.determinant(), 1.0);

        let a = Matrix3::new(1.0, 2.0, 3.0,
                             2.0, 5.0, 3.0,
                             1.0, 3.0, 8.0,);
        assert_eq!(a.determinant(), 8.0);
    }

    #[test]
    fn trace() {
        let one = Matrix3::one();
        assert_eq!(one.trace(), 3.0);

        let a = Matrix3::new(1.0, 2.0, 3.0,
                             2.0, 5.0, 3.0,
                             1.0, 3.0, 8.0,);
        assert_eq!(a.trace(), 14.0);
    }

    #[test]
    fn transposed() {
        let matrix = Matrix3::new(
            1.0, 2.0, 4.0,
            0.0, 1.0, 3.0,
            0.0, 0.0, 1.0
        );

        let transposed = Matrix3::new(
            1.0, 0.0, 0.0,
            2.0, 1.0, 0.0,
            4.0, 3.0, 1.0
        );

        assert_eq!(matrix.transposed(), transposed);
    }
}
