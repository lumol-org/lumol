// Cymbalum, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

//! 3-dimmensional vector type
use std::ops::{Add, Sub, Neg, Mul, Div, BitXor, Index, IndexMut};
use std::cmp::PartialEq;
use super::matrix::Matrix3;

/// 3 dimensional vector type, implementing all usual operations
#[derive(Copy, Clone, Debug)]
pub struct Vector3D([f64; 3]);

impl Vector3D {
    /// Create a new Vector3D with components `x`, `y`, `z`
    pub fn new(x: f64, y: f64, z: f64) -> Vector3D {
        Vector3D([x, y, z])
    }
    /// Return the squared euclidean norm of a Vector3D
    #[inline] pub fn norm2(&self) -> f64 {
        (*self) * (*self)
    }
    /// Return the euclidean norm of a Vector3D
    #[inline] pub fn norm(&self) -> f64 {
        f64::sqrt(self.norm2())
    }
    /// Normalize a Vector3D
    #[inline] pub fn normalized(&self) -> Vector3D {
        *self / self.norm()
    }
    /// Tensorial product between vectors
    pub fn tensorial(&self, other: &Vector3D) -> Matrix3 {
        Matrix3::new(self[0] * other[0], self[0] * other[1], self[0] * other[2],
                     self[1] * other[0], self[1] * other[1], self[1] * other[2],
                     self[2] * other[0], self[2] * other[1], self[2] * other[2])
    }

    /// Convert the vector to a slice
    #[inline] pub fn as_slice(&self) -> &[f64] {
        &self.0
    }
}

/// Add two vectors
impl Add for Vector3D {
    type Output = Vector3D;
    #[inline] fn add(self, other: Vector3D) -> Vector3D {
        Vector3D::new(self[0] + other[0], self[1] + other[1], self[2] + other[2])
    }
}

/// Substract two vectors
impl Sub for Vector3D {
    type Output = Vector3D;
    #[inline] fn sub(self, other: Vector3D) -> Vector3D {
        Vector3D::new(self[0] - other[0], self[1] - other[1], self[2] - other[2])
    }
}

/// Unary - operator
impl Neg for Vector3D {
    type Output = Vector3D;
    #[inline] fn neg(self) -> Vector3D {
        Vector3D::new(-self[0], -self[1], -self[2])
    }
}

/// Multiply by a scalar on the right hand side
impl Mul<f64> for Vector3D {
    type Output = Vector3D;
    #[inline] fn mul(self, other: f64) -> Vector3D {
        Vector3D::new(self[0] * other, self[1] * other, self[2] * other)
    }
}

/// Multiply by a scalar on the left hand side
impl Mul<Vector3D> for f64 {
    type Output = Vector3D;
    #[inline] fn mul(self, other: Vector3D) -> Vector3D {
        Vector3D::new(self * other[0], self * other[1], self * other[2])
    }
}

/// Scalar product between vectors
impl Mul<Vector3D> for Vector3D {
    type Output = f64;
    #[inline] fn mul(self, other: Vector3D) -> f64 {
        self[0] * other[0] + self[1] * other[1] + self[2] * other[2]
    }
}

/// Vectorial product will use the a^b notation.
impl BitXor<Vector3D> for Vector3D {
    type Output = Vector3D;
    fn bitxor(self, other: Vector3D) -> Vector3D {
        let x = self[1] * other[2] - self[2] * other[1];
        let y = self[2] * other[0] - self[0] * other[2];
        let z = self[0] * other[1] - self[1] * other[0];
        Vector3D::new(x, y, z)
    }
}

/// Dividing a vector by a scalar
impl Div<f64> for Vector3D {
    type Output = Vector3D;
    #[inline] fn div(self, other: f64) -> Vector3D {
        Vector3D::new(self[0] / other, self[1] / other, self[2] / other)
    }
}

/// Comparing two vectors
impl PartialEq for Vector3D {
    #[inline] fn eq(&self, other: &Vector3D) -> bool {
        self[0] == other[0] && self[1] == other[1] && self[2] == other[2]
    }
}

/// This is provided for convenience only, and is slower than direct field access
impl Index<usize> for Vector3D {
    type Output = f64;
    #[inline]
    fn index(&self, index: usize) -> &f64 {
        &self.0[index]
    }
}

/// This is provided for convenience only, and is slower than direct field access
impl IndexMut<usize> for Vector3D {
    #[inline]
    fn index_mut(&mut self, index: usize) -> &mut f64 {
        &mut self.0[index]
    }
}

/******************************************************************************/

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn add() {
        let a = Vector3D::new(2.0, 3.5, 4.8);
        let b = Vector3D::new(6.1, -8.5, 7.3);

        let c = a + b;
        assert_eq!(c, Vector3D::new(8.1, -5.0, 12.1));
    }

    #[test]
    fn sub() {
        let a = Vector3D::new(2.0, 3.5, 4.8);
        let b = Vector3D::new(6.1, -8.5, 7.3);

        let c = a - b;
        assert_eq!(c, Vector3D::new(-4.1, 12.0, -2.5));

        let d = -c;
        assert_eq!(d, Vector3D::new(4.1, -12.0, 2.5));
    }

    #[test]
    fn mul() {
        let a = Vector3D::new(2.0, 3.5, 4.8);
        let b = 2.0;

        let c = b * a;
        assert_eq!(c, Vector3D::new(4.0, 7.0, 9.6));

        let b = 1.5;
        let c = a * b;
        assert_eq!(c, Vector3D::new(3.0, 5.25, 7.199999999999999));
    }

    #[test]
    fn dot_product() {
        let a = Vector3D::new(2.1, 3.5, 4.8);
        let b = Vector3D::new(6.1, -8.5, 7.3);

        let c = a * b;
        assert_eq!(c, 18.1);
    }

    #[test]
    fn cross_product() {
        let a = Vector3D::new(2.1, 3.5, 4.8);
        let b = Vector3D::new(6.1, -8.5, 7.3);

        let c = a ^ b;
        assert_eq!(c*a, 0.0);

        let a = Vector3D::new(1.0, 0.0, 0.0);
        let b = Vector3D::new(0.0, 1.0, 0.0);

        let c = a ^ b;
        assert_eq!(c, Vector3D::new(0.0, 0.0, 1.0));
    }

    #[test]
    fn index() {
        let mut a = Vector3D::new(2.1, 3.5, 4.8);

        assert_eq!(a[0], a[0]);
        assert_eq!(a[1], a[1]);
        assert_eq!(a[2], a[2]);

        a[0] = 1.0;
        a[1] = 1.0;
        a[2] = 1.0;

        assert_eq!(a[0], 1.0);
        assert_eq!(a[1], 1.0);
        assert_eq!(a[2], 1.0);
    }

    #[test]
    #[should_panic]
    fn index_out_of_bounds() {
        let mut a = Vector3D::new(2.1, 3.5, 4.8);
        a[3] += 4.0;
    }
}
