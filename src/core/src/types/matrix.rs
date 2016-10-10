// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

//! 3x3 matrix type.
use std::ops::{Add, Sub, Mul, Div, Index, IndexMut};
use std::ops::{AddAssign, SubAssign, MulAssign, DivAssign};
use types::{Vector3D, Zero, One};

/// 3x3 dimensional matrix type, implementing all usual operations
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Matrix3 {
    data : [[f64; 3]; 3]
}

impl Matrix3 {
    /// Create a new `Matrix3` specifying all its components
    #[allow(too_many_arguments)]
    pub fn new(m00: f64, m01: f64, m02: f64,
               m10: f64, m11: f64, m12: f64,
               m20: f64, m21: f64, m22: f64) -> Matrix3 {
        Matrix3{data: [[m00, m01, m02],
                       [m10, m11, m12],
                       [m20, m21, m22]]}
    }

    /// Compute the trace of the matrix
    pub fn trace(&self) -> f64 {
        return self[(0, 0)] + self[(1, 1)] + self[(2, 2)];
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

/******************************************************************************/

#[cfg(test)]
mod tests {
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

        for i in 0..2 {
            for j in 0..2 {
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
        let b = Matrix3::one();
        let add = a + b;

        let res = Matrix3::new(2.0, 2.0, 3.0,
                               4.0, 6.0, 6.0,
                               8.0, 9.0, 11.0);

        for i in 0..3 {
            for j in 0..3 {
                assert_eq!(add[(i, j)], res[(i, j)]);
            }
        }

        a += b;
        for i in 0..3 {
            for j in 0..3 {
                assert_eq!(a[(i, j)], res[(i, j)]);
            }
        }
    }

    #[test]
    fn sub() {
        let mut a = Matrix3::new(1.0, 2.0, 3.0,
                                 4.0, 5.0, 6.0,
                                 8.0, 9.0, 10.0);
        let b = Matrix3::one();
        let sub = a - b;

        let res = Matrix3::new(0.0, 2.0, 3.0,
                               4.0, 4.0, 6.0,
                               8.0, 9.0, 9.0);

        for i in 0..3 {
            for j in 0..3 {
                assert_eq!(sub[(i, j)], res[(i, j)]);
            }
        }

        a -= b;
        for i in 0..3 {
            for j in 0..3 {
                assert_eq!(a[(i, j)], res[(i, j)]);
            }
        }
    }

    #[test]
    fn mul_scalar() {
        let mut a = Matrix3::new(1.0, 2.0, 3.0,
                                 4.0, 5.0, 6.0,
                                 8.0, 9.0, 10.0);
        let b = 2.0;
        let mul_r = a * b;
        let mul_l = b * a;

        let res = Matrix3::new(2.0, 4.0, 6.0,
                               8.0, 10.0, 12.0,
                               16.0, 18.0, 20.0);

        for i in 0..3 {
            for j in 0..3 {
                assert_eq!(mul_r[(i, j)], res[(i, j)]);
                assert_eq!(mul_l[(i, j)], res[(i, j)]);
            }
        }

        a *= b;
        for i in 0..3 {
            for j in 0..3 {
                assert_eq!(a[(i, j)], res[(i, j)]);
            }
        }
    }

    #[test]
    fn div_scalar() {
        let mut a = Matrix3::new(1.0, 2.0, 3.0,
                                 4.0, 5.0, 6.0,
                                 8.0, 9.0, 10.0);
        let b = 2.0;
        let div = a / b;

        let res = Matrix3::new(0.5, 1.0, 1.5,
                               2.0, 2.5, 3.0,
                               4.0, 4.5, 5.0);

        for i in 0..3 {
            for j in 0..3 {
                assert_eq!(div[(i, j)], res[(i, j)]);
            }
        }

        a /= b;
        for i in 0..3 {
            for j in 0..3 {
                assert_eq!(a[(i, j)], res[(i, j)]);
            }
        }
    }

    #[test]
    fn mul_matrix() {
        let unit = Matrix3::one();
        let unit_sq = unit * unit;

        let mut a = Matrix3::new(2.0, 4.0, 6.0,
                                 8.0, 10.0, 12.0,
                                 16.0, 18.0, 20.0);
        let mul_r = unit * a;
        let mul_l = a * unit;

        for i in 0..3 {
            for j in 0..3 {
                // One is the same as (One)*(One)
                assert_eq!(unit[(i, j)], unit_sq[(i, j)]);
                // (One)*A is the same as A and A*(One)
                assert_eq!(mul_l[(i, j)], a[(i, j)]);
                assert_eq!(mul_r[(i, j)], a[(i, j)]);
            }
        }

        let res = a.clone();
        a *= unit;
        for i in 0..3 {
            for j in 0..3 {
                assert_eq!(a[(i, j)], res[(i, j)]);
            }
        }
    }

    #[test]
    fn mul_vector() {
        let a = Matrix3::new(1.0, 2.0, 3.0,
                             4.0, 5.0, 6.0,
                             8.0, 9.0, 10.0);

        let vec = Vector3D::new(1.0, 1.0, 1.0);
        let mul = a * vec;
        let res = Vector3D::new(6.0, 15.0, 27.0);

        assert_eq!(mul[0], res[0]);
        assert_eq!(mul[1], res[1]);
        assert_eq!(mul[2], res[2]);

        let unit = Matrix3::one();
        let vec = Vector3D::new(567.45, 356.8, 215673.12);
        let mul = unit * vec;
        let res = vec;

        assert_eq!(mul[0], res[0]);
        assert_eq!(mul[1], res[1]);
        assert_eq!(mul[2], res[2]);
    }

    #[test]
    fn inverse() {
        let one = Matrix3::one();
        let res = one.inverse();

        for i in 0..3 {
            for j in 0..3 {
                assert_eq!(one[(i, j)], res[(i, j)]);
            }
        }

        let a = Matrix3::new(1.0, 2.0, 3.0,
                             2.0, 5.0, 3.0,
                             1.0, 3.0, 8.0,);
        let b = a.inverse();

        let c = a*b;
        for i in 0..3 {
            for j in 0..3 {
                assert_eq!(one[(i, j)], c[(i, j)]);
            }
        }
    }
}