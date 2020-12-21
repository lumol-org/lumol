// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

//! 3x3 matrix type.
use std::iter;
use std::ops::{Add, Div, Mul, Sub};
use std::ops::{AddAssign, DivAssign, MulAssign, SubAssign};
use std::ops::{Deref, DerefMut};

use num_traits::{One, Zero};

use crate::Vector3D;

/// A 3x3 square matrix type.
///
/// `Matrix3` implements all the usual arithmetic operations:
///
/// ```
/// use lumol_core::types::{Matrix3, Vector3D};
///
/// let one = Matrix3::one();
/// let a = Matrix3::new([
///     [1.0, 0.0, 3.0],
///     [0.0, 2.0, 5.6],
///     [0.0, 0.0, 8.0]
/// ]);
///
/// let v = Vector3D::new(3.0, 2.0, 1.0);
///
/// // Indexing
/// assert_eq!(a[0][0], 1.0);
/// assert_eq!(a[1][2], 5.6);
///
/// // Addition
/// let c = a + one;
/// assert_eq!(c, Matrix3::new([
///     [2.0, 0.0, 3.0],
///     [0.0, 3.0, 5.6],
///     [0.0, 0.0, 9.0]
/// ]));
///
/// // Subtraction
/// let c = a - one;
/// assert_eq!(c, Matrix3::new([
///     [0.0, 0.0, 3.0],
///     [0.0, 1.0, 5.6],
///     [0.0, 0.0, 7.0]
/// ]));
///
/// // Multiplication
/// let c = a * one;  // matrix - matrix
/// assert_eq!(c, a);
///
/// let c = a * v;  // matrix - vector
/// assert_eq!(c, Vector3D::new(6.0, 9.6, 8.0));
///
/// let c = 42.0 * one;  // matrix - scalar
/// assert_eq!(c, Matrix3::new([
///     [42., 0.0, 0.0],
///     [0.0, 42., 0.0],
///     [0.0, 0.0, 42.]
/// ]));
///
/// // Division
/// let c = a / 2.0;
/// assert_eq!(c, Matrix3::new([
///     [0.5, 0.0, 1.5],
///     [0.0, 1.0, 2.8],
///     [0.0, 0.0, 4.0]
/// ]));
/// ```
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Matrix3([[f64; 3]; 3]);

impl Matrix3 {
    /// Create a new `Matrix3` specifying all its components
    /// # Examples
    /// ```
    /// # use lumol_core::types::Matrix3;
    /// let matrix = Matrix3::new([
    ///     [0.0, 0.0, 3.0],
    ///     [0.0, 1.0, 5.6],
    ///     [0.0, 0.0, 7.0],
    /// ]);
    /// assert_eq!(matrix[0][2], 3.0);
    /// ```
    pub fn new(data: [[f64; 3]; 3]) -> Matrix3 {
        Matrix3(data)
    }

    /// Create a new `Matrix3` with components set to `0`
    ///
    /// # Examples
    ///
    /// ```
    /// # use lumol_core::types::Matrix3;
    /// let matrix = Matrix3::zero();
    ///
    /// for i in 0..3 {
    ///     for j in 0..3 {
    ///         assert_eq!(matrix[i][j], 0.0);
    ///     }
    /// }
    /// ```
    pub fn zero() -> Matrix3 {
        <Matrix3 as Zero>::zero()
    }

    /// Create a new `Vector3D` with components all components set to `0`,
    /// except the diagonal components set to 1.
    ///
    /// # Examples
    ///
    /// ```
    /// # use lumol_core::types::Matrix3;
    /// let matrix = Matrix3::one();
    ///
    /// for i in 0..3 {
    ///     for j in 0..3 {
    ///         if i == j {
    ///             assert_eq!(matrix[i][j], 1.0);
    ///         } else {
    ///             assert_eq!(matrix[i][j], 0.0);
    ///         }
    ///     }
    /// }
    /// ```
    pub fn one() -> Matrix3 {
        <Matrix3 as One>::one()
    }

    /// Returns rotation matrix given a rotation angle and an axis.
    ///
    /// # Examples
    ///
    /// ```
    /// # use lumol_core::types::{Vector3D, Matrix3};
    /// let e1 = Vector3D::new(1.0f64, 0.0, 0.0);
    /// let e3 = Vector3D::new(0.0f64, 0.0, 1.0);
    /// let angle = 90f64.to_radians();
    ///
    /// let rotation = Matrix3::rotation(&(e1 * 2.0), angle);
    /// let mut rot_vec = rotation * e3;
    /// for i in 0..3 {
    ///     rot_vec[i] = (rot_vec[i] * 1.0e8).round() / 1.0e8
    /// };
    /// assert_eq!(rot_vec, Vector3D::new(0.0, 1.0, 0.0));
    /// ```
    pub fn rotation(axis: &Vector3D, angle: f64) -> Matrix3 {
        let n = axis.normalized();
        let sin = f64::sin(angle);
        let cos = f64::cos(angle);

        let x_sin = n[0] * sin;
        let y_sin = n[1] * sin;
        let z_sin = n[2] * sin;
        let one_minus_cos = 1.0 - cos;
        let xy = n[0] * n[1] * one_minus_cos;
        let xz = n[0] * n[2] * one_minus_cos;
        let yz = n[1] * n[2] * one_minus_cos;

        // Build the rotation matrix
        Matrix3::new([
            [(n[0] * n[0]) * one_minus_cos + cos, xy + z_sin, xz - y_sin],
            [xy - z_sin, (n[1] * n[1]) * one_minus_cos + cos, yz + x_sin],
            [xz + y_sin, yz - x_sin, (n[2] * n[2]) * one_minus_cos + cos],
        ])
    }

    /// Compute the trace of the matrix
    /// # Examples
    /// ```
    /// # use lumol_core::types::Matrix3;
    /// let matrix = Matrix3::new([
    ///     [0.0, 0.0, 3.0],
    ///     [0.0, 1.0, 5.6],
    ///     [0.0, 0.0, 7.0]
    /// ]);
    /// assert_eq!(matrix.trace(), 8.0);
    /// ```
    pub fn trace(&self) -> f64 {
        return self[0][0] + self[1][1] + self[2][2];
    }

    /// Computes the inverse of a matrix
    ///
    /// # Examples
    ///
    /// ```
    /// # use lumol_core::types::Matrix3;
    /// // A diagonal matrix is trivially invertible
    /// let matrix = Matrix3::new([
    ///     [4.0, 0.0, 0.0],
    ///     [0.0, 1.0, 0.0],
    ///     [0.0, 0.0, 7.0]
    /// ]);
    ///
    /// let inverted = Matrix3::new([
    ///     [1.0/4.0, 0.0, 0.0],
    ///     [0.0, 1.0, 0.0],
    ///     [0.0, 0.0, 1.0/7.0]
    /// ]);
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
        assert!(determinant.abs() > 1e-30, "The matrix is not invertible!");
        let inverse_determinant = 1.0 / determinant;
        let mut res = Matrix3::zero();
        res[0][0] = (self[1][1] * self[2][2] - self[2][1] * self[1][2]) * inverse_determinant;
        res[0][1] = (self[0][2] * self[2][1] - self[0][1] * self[2][2]) * inverse_determinant;
        res[0][2] = (self[0][1] * self[1][2] - self[0][2] * self[1][1]) * inverse_determinant;
        res[1][0] = (self[1][2] * self[2][0] - self[1][0] * self[2][2]) * inverse_determinant;
        res[1][1] = (self[0][0] * self[2][2] - self[0][2] * self[2][0]) * inverse_determinant;
        res[1][2] = (self[1][0] * self[0][2] - self[0][0] * self[1][2]) * inverse_determinant;
        res[2][0] = (self[1][0] * self[2][1] - self[2][0] * self[1][1]) * inverse_determinant;
        res[2][1] = (self[2][0] * self[0][1] - self[0][0] * self[2][1]) * inverse_determinant;
        res[2][2] = (self[0][0] * self[1][1] - self[1][0] * self[0][1]) * inverse_determinant;
        return res;
    }

    /// Computes the [determinant][Wiki] of a matrix
    ///
    /// [Wiki]: https://en.wikipedia.org/wiki/Determinant
    ///
    /// # Examples
    ///
    /// ```
    /// # use lumol_core::types::Matrix3;
    /// let matrix = Matrix3::new([
    ///     [4.0, 0.0, 0.0],
    ///     [0.0, 1.5, 0.0],
    ///     [0.0, 0.0, 7.0]
    /// ]);
    ///
    /// assert_eq!(matrix.determinant(), 4.0 * 1.5 * 7.0);
    /// ```
    pub fn determinant(&self) -> f64 {
          self[0][0] * (self[1][1] * self[2][2] - self[2][1] * self[1][2])
        - self[0][1] * (self[1][0] * self[2][2] - self[1][2] * self[2][0])
        + self[0][2] * (self[1][0] * self[2][1] - self[1][1] * self[2][0])
    }

    /// Transpose this matrix into a new matrix
    ///
    /// # Examples
    ///
    /// ```
    /// # use lumol_core::types::Matrix3;
    /// let matrix = Matrix3::new([
    ///     [1.0, 2.0, 4.0],
    ///     [0.0, 1.0, 3.0],
    ///     [0.0, 0.0, 1.0]
    /// ]);
    ///
    /// let transposed = Matrix3::new([
    ///     [1.0, 0.0, 0.0],
    ///     [2.0, 1.0, 0.0],
    ///     [4.0, 3.0, 1.0]
    /// ]);
    ///
    /// assert_eq!(matrix.transposed(), transposed);
    /// ```
    pub fn transposed(&self) -> Matrix3 {
        Matrix3::new([
            [self[0][0], self[1][0], self[2][0]],
            [self[0][1], self[1][1], self[2][1]],
            [self[0][2], self[1][2], self[2][2]],
        ])
    }

    /// Compute the (Frobenius) norm of the matrix
    ///
    /// # Examples
    ///
    /// ```
    /// # use lumol_core::types::Matrix3;
    /// let matrix = Matrix3::new([
    ///     [1.0, 2.0, 4.0],
    ///     [0.0, 1.0, 3.0],
    ///     [0.0, 0.0, 1.0]
    /// ]);
    ///
    /// assert_eq!(matrix.norm(), 5.656854249492381);
    /// ```
    pub fn norm(&self) -> f64 {
        f64::sqrt(
            self[0][0] * self[0][0] + self[1][0] * self[1][0] + self[2][0] * self[2][0] +
            self[0][1] * self[0][1] + self[1][1] * self[1][1] + self[2][1] * self[2][1] +
            self[0][2] * self[0][2] + self[1][2] * self[1][2] + self[2][2] * self[2][2]
        )
    }
}

impl_arithmetic!(
    Matrix3, Matrix3, Add, add, Matrix3,
    self, other,
    Matrix3::new([
        [self[0][0] + other[0][0], self[0][1] + other[0][1], self[0][2] + other[0][2]],
        [self[1][0] + other[1][0], self[1][1] + other[1][1], self[1][2] + other[1][2]],
        [self[2][0] + other[2][0], self[2][1] + other[2][1], self[2][2] + other[2][2]],
    ])
);

impl_in_place_arithmetic!(
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
    Matrix3::new([
        [self[0][0] - other[0][0], self[0][1] - other[0][1], self[0][2] - other[0][2]],
        [self[1][0] - other[1][0], self[1][1] - other[1][1], self[1][2] - other[1][2]],
        [self[2][0] - other[2][0], self[2][1] - other[2][1], self[2][2] - other[2][2]],
    ])
);

impl_in_place_arithmetic!(
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
        let m00 = self[0][0] * other[0][0] + self[0][1] * other[1][0] + self[0][2] * other[2][0];
        let m01 = self[0][0] * other[0][1] + self[0][1] * other[1][1] + self[0][2] * other[2][1];
        let m02 = self[0][0] * other[0][2] + self[0][1] * other[1][2] + self[0][2] * other[2][2];

        let m10 = self[1][0] * other[0][0] + self[1][1] * other[1][0] + self[1][2] * other[2][0];
        let m11 = self[1][0] * other[0][1] + self[1][1] * other[1][1] + self[1][2] * other[2][1];
        let m12 = self[1][0] * other[0][2] + self[1][1] * other[1][2] + self[1][2] * other[2][2];

        let m20 = self[2][0] * other[0][0] + self[2][1] * other[1][0] + self[2][2] * other[2][0];
        let m21 = self[2][0] * other[0][1] + self[2][1] * other[1][1] + self[2][2] * other[2][1];
        let m22 = self[2][0] * other[0][2] + self[2][1] * other[1][2] + self[2][2] * other[2][2];

        Matrix3::new([
            [m00, m01, m02],
            [m10, m11, m12],
            [m20, m21, m22],
        ])
    }
);

impl_in_place_arithmetic!(
    Matrix3, Matrix3, MulAssign, mul_assign,
    self, other,
    {
        let m00 = self[0][0] * other[0][0] + self[0][1] * other[1][0] + self[0][2] * other[2][0];
        let m01 = self[0][0] * other[0][1] + self[0][1] * other[1][1] + self[0][2] * other[2][1];
        let m02 = self[0][0] * other[0][2] + self[0][1] * other[1][2] + self[0][2] * other[2][2];

        let m10 = self[1][0] * other[0][0] + self[1][1] * other[1][0] + self[1][2] * other[2][0];
        let m11 = self[1][0] * other[0][1] + self[1][1] * other[1][1] + self[1][2] * other[2][1];
        let m12 = self[1][0] * other[0][2] + self[1][1] * other[1][2] + self[1][2] * other[2][2];

        let m20 = self[2][0] * other[0][0] + self[2][1] * other[1][0] + self[2][2] * other[2][0];
        let m21 = self[2][0] * other[0][1] + self[2][1] * other[1][1] + self[2][2] * other[2][1];
        let m22 = self[2][0] * other[0][2] + self[2][1] * other[1][2] + self[2][2] * other[2][2];

        *self = Matrix3::new([
            [m00, m01, m02],
            [m10, m11, m12],
            [m20, m21, m22],
        ]);
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

lsh_scalar_arithmetic!(
    Matrix3, Mul, mul, Matrix3,
    self, other,
    Matrix3::new([
        [self[0][0] * other, self[0][1] * other, self[0][2] * other],
        [self[1][0] * other, self[1][1] * other, self[1][2] * other],
        [self[2][0] * other, self[2][1] * other, self[2][2] * other],
    ])
);

rhs_scalar_arithmetic!(
    Matrix3, Mul, mul, Matrix3,
    self, other,
    Matrix3::new([
        [self * other[0][0], self * other[0][1], self * other[0][2]],
        [self * other[1][0], self * other[1][1], self * other[1][2]],
        [self * other[2][0], self * other[2][1], self * other[2][2]],
    ])
);

impl_in_place_arithmetic!(
    Matrix3, f64, MulAssign, mul_assign,
    self, other,
    {
        #[allow(clippy::clone_on_copy)]
        let other = other.clone();
        self[0][0] *= other; self[0][1] *= other; self[0][2] *= other;
        self[1][0] *= other; self[1][1] *= other; self[1][2] *= other;
        self[2][0] *= other; self[2][1] *= other; self[2][2] *= other;
    }
);

lsh_scalar_arithmetic!(
    Matrix3, Div, div, Matrix3,
    self, other,
    Matrix3::new([
        [self[0][0] / other, self[0][1] / other, self[0][2] / other],
        [self[1][0] / other, self[1][1] / other, self[1][2] / other],
        [self[2][0] / other, self[2][1] / other, self[2][2] / other],
    ])
);

impl_in_place_arithmetic!(
    Matrix3, f64, DivAssign, div_assign,
    self, other,
    {
        #[allow(clippy::clone_on_copy)]
        let other = other.clone();
        self[0][0] /= other; self[0][1] /= other; self[0][2] /= other;
        self[1][0] /= other; self[1][1] /= other; self[1][2] /= other;
        self[2][0] /= other; self[2][1] /= other; self[2][2] /= other;
    }
);

impl Zero for Matrix3 {
    fn zero() -> Matrix3 {
        Matrix3::new([[0.0; 3]; 3])
    }

    fn is_zero(&self) -> bool {
        *self == Matrix3::zero()
    }
}

impl One for Matrix3 {
    /// Create an identity matrix
    fn one() -> Matrix3 {
        Matrix3::new([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    }
}

impl iter::Sum for Matrix3 {
    fn sum<I>(iter: I) -> Self
    where
        I: Iterator<Item = Matrix3>,
    {
        iter.fold(Matrix3::zero(), |a, b| a + b)
    }
}

impl Deref for Matrix3 {
    type Target = [[f64; 3]; 3];
    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for Matrix3 {
    #[inline]
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl From<[[f64; 3]; 3]> for Matrix3 {
    fn from(data: [[f64; 3]; 3]) -> Matrix3 {
        Matrix3(data)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Vector3D;

    use approx::{AbsDiffEq, RelativeEq, UlpsEq, assert_ulps_eq};

    impl AbsDiffEq for Matrix3 {
        type Epsilon = <f64 as AbsDiffEq>::Epsilon;

        fn default_epsilon() -> Self::Epsilon {
            f64::default_epsilon()
        }

        fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
            f64::abs_diff_eq(&self[0][0], &other[0][0], epsilon) &&
            f64::abs_diff_eq(&self[0][1], &other[0][1], epsilon) &&
            f64::abs_diff_eq(&self[0][2], &other[0][2], epsilon) &&

            f64::abs_diff_eq(&self[1][0], &other[1][0], epsilon) &&
            f64::abs_diff_eq(&self[1][1], &other[1][1], epsilon) &&
            f64::abs_diff_eq(&self[1][2], &other[1][2], epsilon) &&

            f64::abs_diff_eq(&self[2][0], &other[2][0], epsilon) &&
            f64::abs_diff_eq(&self[2][1], &other[2][1], epsilon) &&
            f64::abs_diff_eq(&self[2][2], &other[2][2], epsilon)
        }
    }

    impl RelativeEq for Matrix3 {
        fn default_max_relative() -> Self::Epsilon {
            f64::default_max_relative()
        }

        fn relative_eq(&self, other: &Self, epsilon: Self::Epsilon, max_relative: Self::Epsilon) -> bool {
            f64::relative_eq(&self[0][0], &other[0][0], epsilon, max_relative) &&
            f64::relative_eq(&self[0][1], &other[0][1], epsilon, max_relative) &&
            f64::relative_eq(&self[0][2], &other[0][2], epsilon, max_relative) &&

            f64::relative_eq(&self[1][0], &other[1][0], epsilon, max_relative) &&
            f64::relative_eq(&self[1][1], &other[1][1], epsilon, max_relative) &&
            f64::relative_eq(&self[1][2], &other[1][2], epsilon, max_relative) &&

            f64::relative_eq(&self[2][0], &other[2][0], epsilon, max_relative) &&
            f64::relative_eq(&self[2][1], &other[2][1], epsilon, max_relative) &&
            f64::relative_eq(&self[2][2], &other[2][2], epsilon, max_relative)
        }
    }

    impl UlpsEq for Matrix3 {
        fn default_max_ulps() -> u32 {
            f64::default_max_ulps()
        }

        fn ulps_eq(&self, other: &Self, epsilon: Self::Epsilon, max_ulps: u32) -> bool {
            f64::ulps_eq(&self[0][0], &other[0][0], epsilon, max_ulps) &&
            f64::ulps_eq(&self[0][1], &other[0][1], epsilon, max_ulps) &&
            f64::ulps_eq(&self[0][2], &other[0][2], epsilon, max_ulps) &&

            f64::ulps_eq(&self[1][0], &other[1][0], epsilon, max_ulps) &&
            f64::ulps_eq(&self[1][1], &other[1][1], epsilon, max_ulps) &&
            f64::ulps_eq(&self[1][2], &other[1][2], epsilon, max_ulps) &&

            f64::ulps_eq(&self[2][0], &other[2][0], epsilon, max_ulps) &&
            f64::ulps_eq(&self[2][1], &other[2][1], epsilon, max_ulps) &&
            f64::ulps_eq(&self[2][2], &other[2][2], epsilon, max_ulps)
        }
    }

    #[test]
    fn specials_matrix() {
        let a = Matrix3::zero();
        let b = Matrix3::one();

        for i in 0..3 {
            for j in 0..3 {
                assert_eq!(a[i][j], 0.0);
                if i == j {
                    assert_eq!(b[i][j], 1.0);
                } else {
                    assert_eq!(b[i][j], 0.0);
                }
            }
        }
    }

    #[test]
    fn add() {
        let mut a = Matrix3::new([
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [8.0, 9.0, 10.0]
        ]);
        let res = Matrix3::new([
            [2.0, 2.0, 3.0],
            [4.0, 6.0, 6.0],
            [8.0, 9.0, 11.0]
        ]);

        assert_eq!(a + Matrix3::one(), res);

        a += Matrix3::one();
        assert_eq!(a, res);
    }

    #[test]
    fn sub() {
        let mut a = Matrix3::new([
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [8.0, 9.0, 10.0]
        ]);
        let res = Matrix3::new([
            [0.0, 2.0, 3.0],
            [4.0, 4.0, 6.0],
            [8.0, 9.0, 9.0]
        ]);

        assert_eq!(a - Matrix3::one(), res);

        a -= Matrix3::one();
        assert_eq!(a, res);
    }

    #[test]
    fn mul_scalar() {
        let mut a = Matrix3::new([
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [8.0, 9.0, 10.0],
        ]);
        let res = Matrix3::new([
            [2.0, 4.0, 6.0],
            [8.0, 10.0, 12.0],
            [16.0, 18.0, 20.0],
        ]);
        assert_eq!(a * 2.0, res);
        assert_eq!(2.0 * a, res);

        a *= 2.0;
        assert_eq!(a, res);
    }

    #[test]
    fn div_scalar() {
        let mut a = Matrix3::new([
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [8.0, 9.0, 10.0],
        ]);
        let res = Matrix3::new([
            [0.5, 1.0, 1.5],
            [2.0, 2.5, 3.0],
            [4.0, 4.5, 5.0]
        ]);
        assert_eq!(a / 2.0, res);

        a /= 2.0;
        assert_eq!(a, res);
    }

    #[test]
    fn mul_matrix() {
        let unit = Matrix3::one();
        let mut a = Matrix3::new([
            [2.0, 4.0, 6.0],
            [8.0, 10.0, 12.0],
            [16.0, 18.0, 20.0],
        ]);

        assert_eq!(unit, unit * unit);
        assert_eq!(unit * a, a);
        assert_eq!(a * unit, a);

        let previous = a;
        a *= unit;
        assert_eq!(previous, a);
    }

    #[test]
    #[allow(clippy::unreadable_literal)]
    fn mul_vector() {
        let a = Matrix3::new([
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [8.0, 9.0, 10.0],
        ]);

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

        let a = Matrix3::new([
            [1.0, 2.0, 3.0],
            [2.0, 5.0, 3.0],
            [1.0, 3.0, 8.0],
        ]);
        assert!(a.determinant() != 0.0);
        assert_eq!(one, a * a.inverse());
    }

    #[test]
    fn determinant() {
        let one = Matrix3::one();
        assert_eq!(one.determinant(), 1.0);

        let a = Matrix3::new([
            [1.0, 2.0, 3.0],
            [2.0, 5.0, 3.0],
            [1.0, 3.0, 8.0],
        ]);
        assert_eq!(a.determinant(), 8.0);
    }

    #[test]
    fn trace() {
        let one = Matrix3::one();
        assert_eq!(one.trace(), 3.0);

        let a = Matrix3::new([
            [1.0, 2.0, 3.0],
            [2.0, 5.0, 3.0],
            [1.0, 3.0, 8.0],
        ]);
        assert_eq!(a.trace(), 14.0);
    }

    #[test]
    fn transposed() {
        let matrix = Matrix3::new([
            [1.0, 2.0, 4.0],
            [0.0, 1.0, 3.0],
            [0.0, 0.0, 1.0]
        ]);

        let transposed = Matrix3::new([
            [1.0, 0.0, 0.0],
            [2.0, 1.0, 0.0],
            [4.0, 3.0, 1.0]
        ]);

        assert_eq!(matrix.transposed(), transposed);
    }

    #[test]
    fn iter_sum() {
        let a = Matrix3::new([
            [1.0, 2.0, 4.0],
            [0.0, 1.0, 3.0],
            [0.0, 0.0, 1.0],
        ]);
        let b = Matrix3::new([
            [7.0, -2.0, 5.0],
            [1.0, 1.0, 3.0],
            [0.0, -8.0, 3.0],
        ]);

        let vec = vec![a, b, a];
        let sum : Matrix3 = vec.into_iter().sum();
        assert_eq!(sum, a + a + b);
    }

    #[test]
    #[allow(clippy::unreadable_literal)]
    fn norm() {
        let matrix = Matrix3::new([
            [1.0, 2.0, 4.0],
            [0.0, 1.0, 3.0],
            [-2.0, 2.2, 1.0]
        ]);

        assert_ulps_eq!(matrix.norm(), 6.390618123468183);

        assert_eq!(Matrix3::zero().norm(), 0.0);
        assert_eq!(Matrix3::one().norm(), f64::sqrt(3.0));
    }
}
