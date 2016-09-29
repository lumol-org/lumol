// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

//! 3-dimmensional vector type
use std::ops::{Add, Sub, Neg, Mul, Div, BitXor, Index, IndexMut};
use std::ops::{AddAssign, SubAssign, MulAssign, DivAssign};
use std::cmp::PartialEq;
use types::{Matrix3, Zero};

/// A simple 3-dimensional vector type, storing three `f64`.
///
/// A `Vector3D` implement all the arithmetics operations:
/// ```rust
/// // `u`, `v` and `w` are vectors. `a` and `b` are scalar
/// let u = Vector3D::new(1.0, 2.0, 3.0);
/// let v = Vector3D::new(1.0, 2.0, 3.0);
/// let a = 42.0;
///
/// // Addition
/// let w = u + v;
/// // Substraction
/// let w = u - v;
/// // Negation
/// let w = -u;
/// // Cross product
/// let w = u ^ v;
///
/// // Multiplication
/// let w = a * u;
/// let w = v * a;
/// // Division
/// let w = u / a;
///
/// // Dot product
/// let b = u * v;
/// ```
#[derive(Copy, Clone, Debug)]
pub struct Vector3D([f64; 3]);

impl Vector3D {
    /// Create a new `Vector3D` with components `x`, `y`, `z`
    pub fn new(x: f64, y: f64, z: f64) -> Vector3D {
        Vector3D([x, y, z])
    }
    /// Return the squared euclidean norm of a Vector3D
    #[inline] pub fn norm2(&self) -> f64 {
        self * self
    }
    /// Return the euclidean norm of a Vector3D
    #[inline] pub fn norm(&self) -> f64 {
        f64::sqrt(self.norm2())
    }
    /// Normalize a Vector3D
    #[inline] pub fn normalized(&self) -> Vector3D {
        self / self.norm()
    }
    /// Tensorial product between vectors
    pub fn tensorial(&self, other: &Vector3D) -> Matrix3 {
        Matrix3::new(self[0] * other[0], self[0] * other[1], self[0] * other[2],
                     self[1] * other[0], self[1] * other[1], self[1] * other[2],
                     self[2] * other[0], self[2] * other[1], self[2] * other[2])
    }
}

impl_arithmetic!(
    Vector3D, Vector3D, Add, add, Vector3D,
    self, other,
    Vector3D::new(self[0] + other[0], self[1] + other[1], self[2] + other[2])
);

impl_inplace_arithmetic!(
    Vector3D, Vector3D, AddAssign, add_assign,
    self, other,
    {self[0] += other[0]; self[1] += other[1]; self[2] += other[2]}
);

impl_arithmetic!(
    Vector3D, Vector3D, Sub, sub, Vector3D,
    self, other,
    Vector3D::new(self[0] - other[0], self[1] - other[1], self[2] - other[2])
);

impl_inplace_arithmetic!(
    Vector3D, Vector3D, SubAssign, sub_assign,
    self, other,
    {self[0] -= other[0]; self[1] -= other[1]; self[2] -= other[2]}
);

// Dot product
impl_arithmetic!(
    Vector3D, Vector3D, Mul, mul, f64,
    self, other,
    self[0] * other[0] + self[1] * other[1] + self[2] * other[2]
);

// Cross product
impl_arithmetic!(
    Vector3D, Vector3D, BitXor, bitxor, Vector3D,
    self, other,
    {let x = self[1] * other[2] - self[2] * other[1];
    let y = self[2] * other[0] - self[0] * other[2];
    let z = self[0] * other[1] - self[1] * other[0];
    Vector3D::new(x, y, z)}
);

/******************************************************************************/

lsh_scal_arithmetic!(
    Vector3D, Mul, mul, Vector3D,
    self, other,
    Vector3D::new(self[0] * other, self[1] * other, self[2] * other)
);

rhs_scal_arithmetic!(
    Vector3D, Mul, mul, Vector3D,
    self, other,
    Vector3D::new(self * other[0], self * other[1], self * other[2])
);

impl_inplace_arithmetic!(
    Vector3D, f64, MulAssign, mul_assign,
    self, other,
    {let other = other.clone(); self[0] *= other; self[1] *= other; self[2] *= other}
);

lsh_scal_arithmetic!(
    Vector3D, Div, div, Vector3D,
    self, other,
    Vector3D::new(self[0] / other, self[1] / other, self[2] / other)
);

impl_inplace_arithmetic!(
    Vector3D, f64, DivAssign, div_assign,
    self, other,
    {let other = other.clone(); self[0] /= other; self[1] /= other; self[2] /= other}
);

/******************************************************************************/

impl Neg for Vector3D {
    type Output = Vector3D;
    #[inline] fn neg(self) -> Vector3D {
        Vector3D::new(-self[0], -self[1], -self[2])
    }
}

impl<'a> Neg for &'a Vector3D {
    type Output = Vector3D;
    #[inline] fn neg(self) -> Vector3D {
        Vector3D::new(-self[0], -self[1], -self[2])
    }
}

impl<'a> Neg for &'a mut Vector3D {
    type Output = Vector3D;
    #[inline] fn neg(self) -> Vector3D {
        Vector3D::new(-self[0], -self[1], -self[2])
    }
}

/******************************************************************************/

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

impl Zero for Vector3D {
    fn zero() -> Vector3D {
        Vector3D::new(0.0, 0.0, 0.0)
    }

    fn is_zero(&self) -> bool {
        self.norm2() == 0.0
    }
}

/******************************************************************************/

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn add() {
        let mut a = Vector3D::new(2.0, 3.5, 4.8);
        let mut b = Vector3D::new(6.1, -8.5, 7.3);

        let c = a + b;
        assert_eq!(c, Vector3D::new(8.1, -5.0, 12.1));

        a += b;
        assert_eq!(a, Vector3D::new(8.1, -5.0, 12.1));

        // Just checking that everything compile
        let _ = &a + b;
        let _ = a + &b;
        let _ = &a + &b;
        let _ = a + &mut b;
        let _ = &mut a + b;
        let _ = &mut a + &mut b;
        let _ = &mut a + &b;
        let _ = &a + &mut b;
        a += &b;
        a += &mut b;
        let _ = a;
    }

    #[test]
    fn sub() {
        let mut a = Vector3D::new(2.0, 3.5, 4.8);
        let mut b = Vector3D::new(6.1, -8.5, 7.3);

        let c = a - b;
        assert_eq!(c, Vector3D::new(-4.1, 12.0, -2.5));

        a -= b;
        assert_eq!(a, Vector3D::new(-4.1, 12.0, -2.5));

        // Just checking that everything compile
        let _ = &a - b;
        let _ = a - &b;
        let _ = &a - &b;
        let _ = a - &mut b;
        let _ = &mut a - b;
        let _ = &mut a - &mut b;
        let _ = &mut a - &b;
        let _ = &a - &mut b;
        a -= &b;
        a -= &mut b;
        let _ = a;
    }

    #[test]
    fn neg() {
        let mut a = Vector3D::new(6.1, -8.5, 7.3);

        let b = -a;
        assert_eq!(b, Vector3D::new(-6.1, 8.5, -7.3));

        let _ = -&a;
        let _ = -&mut a;
    }

    #[test]
    fn mul() {
        let mut a = Vector3D::new(2.0, 3.5, 4.8);
        let b = 2.0;

        let c = b * a;
        assert_eq!(c, Vector3D::new(4.0, 7.0, 9.6));

        let _ = b * &a;
        let _ = b * &mut a;

        let mut b = 1.5;
        let c = a * b;
        assert_eq!(c, Vector3D::new(3.0, 5.25, 7.199999999999999));

        a *= b;
        assert_eq!(a, Vector3D::new(3.0, 5.25, 7.199999999999999));

        // Just checking that everything compile
        let _ = &a * b;
        let _ = &mut a * b;
        a *= &b;
        a *= &mut b;
        let _ = a;
    }

    #[test]
    fn div() {
        let mut a = Vector3D::new(2.0, 3.5, 4.8);
        let mut b = 2.0;
        let c = a / b;
        assert_eq!(c, Vector3D::new(1.0, 1.75, 2.4));

        a /= b;
        assert_eq!(a, Vector3D::new(1.0, 1.75, 2.4));

        // Just checking that everything compile
        let _ = &a / b;
        let _ = &mut a / b;
        a /= &b;
        a /= &mut b;
        let _ = a;
    }

    #[test]
    fn dot_product() {
        let mut a = Vector3D::new(2.1, 3.5, 4.8);
        let mut b = Vector3D::new(6.1, -8.5, 7.3);

        let c = a * b;
        assert_eq!(c, 18.1);

        // Just checking that everything compile
        let _ = &a * b;
        let _ = a * &b;
        let _ = &a * &b;
        let _ = a * &mut b;
        let _ = &mut a * b;
        let _ = &mut a * &mut b;
        let _ = &mut a * &b;
        let _ = &a * &mut b;
    }

    #[test]
    fn cross_product() {
        let a = Vector3D::new(2.1, 3.5, 4.8);
        let b = Vector3D::new(6.1, -8.5, 7.3);

        let c = a ^ b;
        assert_eq!(c*a, 0.0);

        let mut a = Vector3D::new(1.0, 0.0, 0.0);
        let mut b = Vector3D::new(0.0, 1.0, 0.0);

        let c = a ^ b;
        assert_eq!(c, Vector3D::new(0.0, 0.0, 1.0));

        // Just checking that everything compile
        let _ = &a ^ b;
        let _ = a ^ &b;
        let _ = &a ^ &b;
        let _ = a ^ &mut b;
        let _ = &mut a ^ b;
        let _ = &mut a ^ &mut b;
        let _ = &mut a ^ &b;
        let _ = &a ^ &mut b;
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
