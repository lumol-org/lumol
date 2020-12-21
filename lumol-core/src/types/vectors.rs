// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

//! 3-dimensional vector type
use std::ops::{Add, BitXor, Div, Mul, Neg, Sub};
use std::ops::{AddAssign, DivAssign, MulAssign, SubAssign};
use std::ops::{Deref, DerefMut};

use num_traits::Zero;

use crate::Matrix3;

/// A 3-dimensional vector type
///
/// A `Vector3D` implement all the arithmetic operations:
///
/// ```
/// # use lumol_core::types::Vector3D;
/// let u = Vector3D::new(1.0, 2.0, 3.0);
/// let v = Vector3D::new(4.0, -2.0, 1.0);
///
/// // Indexing
/// assert_eq!(u[0], 1.0);
/// assert_eq!(u[1], 2.0);
/// assert_eq!(u[2], 3.0);
///
/// // Addition
/// let w = u + v;
/// assert_eq!(w, Vector3D::new(5.0, 0.0, 4.0));
///
/// // Subtraction
/// let w = u - v;
/// assert_eq!(w, Vector3D::new(-3.0, 4.0, 2.0));
///
/// // Negation
/// let w = -u;
/// assert_eq!(w, Vector3D::new(-1.0, -2.0, -3.0));
///
/// // Cross product
/// let w = u ^ v;
/// assert_eq!(w, Vector3D::new(8.0, 11.0, -10.0));
///
/// // Multiplication
/// let w = 2.0 * u;
/// assert_eq!(w, Vector3D::new(2.0, 4.0, 6.0));
///
/// let w = u * 3.0;
/// assert_eq!(w, Vector3D::new(3.0, 6.0, 9.0));
///
/// // Division
/// let w = u / 2.0;
/// assert_eq!(w, Vector3D::new(0.5, 1.0, 1.5));
///
/// // Dot product
/// let a = u * v;
/// assert_eq!(a, 3.0);
/// ```
#[derive(Copy, Clone, Debug, PartialEq, PartialOrd)]
pub struct Vector3D([f64; 3]);

impl Vector3D {
    /// Create a new `Vector3D` with components `x`, `y`, `z`
    ///
    /// # Examples
    ///
    /// ```
    /// # use lumol_core::types::Vector3D;
    /// let vector = Vector3D::new(1.0, 0.0, -42.0);
    ///
    /// assert_eq!(vector[0], 1.0);
    /// assert_eq!(vector[1], 0.0);
    /// assert_eq!(vector[2], -42.0);
    /// ```
    pub fn new(x: f64, y: f64, z: f64) -> Vector3D {
        Vector3D([x, y, z])
    }

    /// Create a new `Vector3D` with components `0`, `0`, `0`
    ///
    /// # Examples
    ///
    /// ```
    /// # use lumol_core::types::Vector3D;
    /// let vector = Vector3D::zero();
    /// assert_eq!(vector[0], 0.0);
    /// assert_eq!(vector[1], 0.0);
    /// assert_eq!(vector[2], 0.0);
    /// ```
    pub fn zero() -> Vector3D {
        <Vector3D as Zero>::zero()
    }

    /// Return the squared euclidean norm of a `Vector3D`
    ///
    /// # Examples
    /// ```
    /// # use lumol_core::types::Vector3D;
    /// let vec = Vector3D::new(1.0, 0.0, -4.0);
    /// assert_eq!(vec.norm2(), 17.0);
    /// ```
    #[inline]
    pub fn norm2(&self) -> f64 {
        self * self
    }

    /// Return the euclidean norm of a `Vector3D`
    /// # Examples
    /// ```
    /// # use lumol_core::types::Vector3D;
    /// # use std::f64;
    /// let vec = Vector3D::new(1.0, 0.0, -4.0);
    /// assert_eq!(vec.norm(), f64::sqrt(17.0));
    /// ```
    #[inline]
    pub fn norm(&self) -> f64 {
        f64::sqrt(self.norm2())
    }

    /// Normalize a `Vector3D`.
    /// # Examples
    /// ```
    /// # use lumol_core::types::Vector3D;
    /// let vec = Vector3D::new(1.0, 0.0, -4.0);
    /// let n = vec.normalized();
    /// assert_eq!(n.norm(), 1.0);
    /// ```
    #[inline]
    pub fn normalized(&self) -> Vector3D {
        self / self.norm()
    }

    /// Tensorial product between vectors. The tensorial product between the
    /// vectors `a` and `b` creates a `Matrix3` with component (i, j) equals to
    /// `a[i] * b[j]`.
    ///
    /// # Examples
    ///
    /// ```
    /// # use lumol_core::types::Vector3D;
    /// # use lumol_core::types::Matrix3;
    /// let a = Vector3D::new(1.0, 0.0, -4.0);
    /// let b = Vector3D::new(1.0, 2.0, 3.0);
    /// let matrix = Matrix3::new([
    ///     [1.0, 2.0, 3.0],
    ///     [0.0, 0.0, 0.0],
    ///     [-4.0, -8.0, -12.0]
    /// ]);
    /// assert_eq!(a.tensorial(&b), matrix);
    /// ```
    pub fn tensorial(&self, other: &Vector3D) -> Matrix3 {
        Matrix3::new([
            [self[0] * other[0], self[0] * other[1], self[0] * other[2]],
            [self[1] * other[0], self[1] * other[1], self[1] * other[2]],
            [self[2] * other[0], self[2] * other[1], self[2] * other[2]],
        ])
    }

    /// Get the minimal value in this vector, using `std::f64::min` for
    /// comparison.
    ///
    /// # Examples
    ///
    /// ```
    /// # use lumol_core::types::Vector3D;
    /// let vector = Vector3D::new(1.0, 0.0, -4.0);
    ///
    /// assert_eq!(vector.min(), -4.0);
    /// ```
    pub fn min(&self) -> f64 {
        f64::min(f64::min(self[0], self[1]), self[2])
    }

    /// Get the maximal value in this vector, using `std::f64::max` for
    /// comparison.
    ///
    /// # Examples
    ///
    /// ```
    /// # use lumol_core::types::Vector3D;
    /// let vector = Vector3D::new(1.0, 0.0, -4.0);
    ///
    /// assert_eq!(vector.max(), 1.0);
    /// ```
    pub fn max(&self) -> f64 {
        f64::max(f64::max(self[0], self[1]), self[2])
    }
}

impl_arithmetic!(
    Vector3D, Vector3D, Add, add, Vector3D,
    self, other,
    Vector3D::new(self[0] + other[0], self[1] + other[1], self[2] + other[2])
);

impl_in_place_arithmetic!(
    Vector3D, Vector3D, AddAssign, add_assign,
    self, other,
    {self[0] += other[0]; self[1] += other[1]; self[2] += other[2]}
);

impl_arithmetic!(
    Vector3D, Vector3D, Sub, sub, Vector3D,
    self, other,
    Vector3D::new(self[0] - other[0], self[1] - other[1], self[2] - other[2])
);

impl_in_place_arithmetic!(
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

lsh_scalar_arithmetic!(
    Vector3D, Mul, mul, Vector3D,
    self, other,
    Vector3D::new(self[0] * other, self[1] * other, self[2] * other)
);

rhs_scalar_arithmetic!(
    Vector3D, Mul, mul, Vector3D,
    self, other,
    Vector3D::new(self * other[0], self * other[1], self * other[2])
);

impl_in_place_arithmetic!(
    Vector3D, f64, MulAssign, mul_assign,
    self, other,
    {
        #[allow(clippy::clone_on_copy)]
        let other = other.clone();
        self[0] *= other;
        self[1] *= other;
        self[2] *= other;
    }
);

lsh_scalar_arithmetic!(
    Vector3D, Div, div, Vector3D,
    self, other,
    Vector3D::new(self[0] / other, self[1] / other, self[2] / other)
);

impl_in_place_arithmetic!(
    Vector3D, f64, DivAssign, div_assign,
    self, other,
    {
        #[allow(clippy::clone_on_copy)]
        let other = other.clone();
        self[0] /= other;
        self[1] /= other;
        self[2] /= other;
    }
);

impl Neg for Vector3D {
    type Output = Vector3D;
    #[inline]
    fn neg(self) -> Vector3D {
        Vector3D::new(-self[0], -self[1], -self[2])
    }
}

impl<'a> Neg for &'a Vector3D {
    type Output = Vector3D;
    #[inline]
    fn neg(self) -> Vector3D {
        Vector3D::new(-self[0], -self[1], -self[2])
    }
}

impl<'a> Neg for &'a mut Vector3D {
    type Output = Vector3D;
    #[inline]
    fn neg(self) -> Vector3D {
        Vector3D::new(-self[0], -self[1], -self[2])
    }
}

impl Deref for Vector3D {
    type Target = [f64; 3];

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for Vector3D {
    #[inline]
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
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

impl From<[f64; 3]> for Vector3D {
    fn from(data: [f64; 3]) -> Self {
        Vector3D(data)
    }
}

impl Default for Vector3D {
    fn default() -> Self {
        Vector3D::zero()
    }
}

#[cfg(test)]
#[allow(clippy::op_ref)]
mod tests {
    use crate::{Matrix3, Vector3D};
    use std::f64;

    use approx::{AbsDiffEq, RelativeEq, UlpsEq};

    impl AbsDiffEq for Vector3D {
        type Epsilon = <f64 as AbsDiffEq>::Epsilon;

        fn default_epsilon() -> Self::Epsilon {
            f64::default_epsilon()
        }

        fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
            f64::abs_diff_eq(&self[0], &other[0], epsilon) &&
            f64::abs_diff_eq(&self[1], &other[1], epsilon) &&
            f64::abs_diff_eq(&self[2], &other[2], epsilon)
        }
    }

    impl RelativeEq for Vector3D {
        fn default_max_relative() -> Self::Epsilon {
            f64::default_max_relative()
        }

        fn relative_eq(&self, other: &Self, epsilon: Self::Epsilon, max_relative: Self::Epsilon) -> bool {
            f64::relative_eq(&self[0], &other[0], epsilon, max_relative) &&
            f64::relative_eq(&self[1], &other[1], epsilon, max_relative) &&
            f64::relative_eq(&self[2], &other[2], epsilon, max_relative)
        }
    }

    impl UlpsEq for Vector3D {
        fn default_max_ulps() -> u32 {
            f64::default_max_ulps()
        }

        fn ulps_eq(&self, other: &Self, epsilon: Self::Epsilon, max_ulps: u32) -> bool {
            f64::ulps_eq(&self[0], &other[0], epsilon, max_ulps) &&
            f64::ulps_eq(&self[1], &other[1], epsilon, max_ulps) &&
            f64::ulps_eq(&self[2], &other[2], epsilon, max_ulps)
        }
    }

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
    #[allow(clippy::unreadable_literal)]
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
        assert_eq!(c * a, 0.0);

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
    fn tensorial() {
        let a = Vector3D::new(1.0, 0.0, -4.0);
        let b = Vector3D::new(1.0, 2.0, 3.0);
        let matrix = Matrix3::new([[1.0, 2.0, 3.0], [0.0, 0.0, 0.0], [-4.0, -8.0, -12.0]]);
        assert_eq!(a.tensorial(&b), matrix);
        assert_eq!(b.tensorial(&a), matrix.transposed());
    }

    #[test]
    fn min() {
        assert_eq!(Vector3D::new(1.0, 0.0, -4.0).min(), -4.0);
        assert_eq!(Vector3D::new(1.0, f64::MIN, 4.0).min(), f64::MIN);

        assert_eq!(Vector3D::new(f64::NAN, 8.0, f64::NAN).min(), 8.0);
        assert_eq!(Vector3D::new(-4.0e89, 0.0, -f64::INFINITY).min(), -f64::INFINITY);
        assert_eq!(Vector3D::new(f64::MIN, 0.0, -f64::INFINITY).min(), -f64::INFINITY);
    }

    #[test]
    fn max() {
        assert_eq!(Vector3D::new(1.0, 0.0, -4.0).max(), 1.0);
        assert_eq!(Vector3D::new(1.0, f64::MAX, -4.0).max(), f64::MAX);

        assert_eq!(Vector3D::new(f64::NAN, 8.0, f64::NAN).max(), 8.0);
        assert_eq!(Vector3D::new(4.0e89, 0.0, f64::INFINITY).max(), f64::INFINITY);
        assert_eq!(Vector3D::new(f64::MAX, 0.0, f64::INFINITY).max(), f64::INFINITY);
    }
}
