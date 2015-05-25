/*
 * Cymbalum, Molecular Simulation in Rust
 * Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
*/

use std::ops::{Add, Sub, Mul, Div};
use super::matrix::Matrix3;

/// 3 dimensional vector type, implementing all usual operations
#[derive(Copy, Clone)]
pub struct Vector3D {
    /// First component of the vector
    pub x: f64,
    /// Second component of the vector
    pub y: f64,
    /// Third component of the vector
    pub z: f64,
}

impl Vector3D {
    /// Create a new Vector3D with components `x`, `y`, `z`
    pub fn new(x: f64, y: f64, z: f64) -> Vector3D {
        Vector3D{x: x, y: y, z: z}
    }
    /// Return the squared euclidean norm of a Vector3D
    #[inline]
    pub fn norm2(&self) -> f64 {
        (*self) * (*self)
    }
    /// Return the euclidean norm of a Vector3D
    #[inline]
    pub fn norm(&self) -> f64 {
        f64::sqrt(self.norm2())
    }
    /// Normalize a Vector3D
    #[inline]
    pub fn normalize(&self) -> Vector3D {
        *self / self.norm()
    }
    /// Tensorial product between vectors
    pub fn tensorial(&self, other: &Vector3D) -> Matrix3 {
        Matrix3::new(self.x * other.x, self.x * other.y, self.x * other.z,
                     self.y * other.x, self.y * other.y, self.y * other.z,
                     self.z * other.x, self.z * other.y, self.z * other.z)
    }
}

/// Add two vectors
impl Add for Vector3D {
    type Output = Vector3D;
    fn add(self, other: Vector3D) -> Vector3D {
        Vector3D::new(self.x + other.x, self.y + other.y, self.z + other.z)
    }
}

/// Substract two vectors
impl Sub for Vector3D {
    type Output = Vector3D;
    fn sub(self, other: Vector3D) -> Vector3D {
        Vector3D::new(self.x - other.x, self.y - other.y, self.z - other.z)
    }
}

/// Multiply by a scalar on the right hand side
impl Mul<f64> for Vector3D {
    type Output = Vector3D;
    fn mul(self, other: f64) -> Vector3D {
        Vector3D::new(self.x * other, self.y * other, self.z * other)
    }
}

/// Multiply by a scalar on the left hand side
impl Mul<Vector3D> for f64 {
    type Output = Vector3D;
    fn mul(self, other: Vector3D) -> Vector3D {
        Vector3D::new(self * other.x, self * other.y, self * other.z)
    }
}

/// Scalar product between vectors
impl Mul<Vector3D> for Vector3D {
    type Output = f64;
    fn mul(self, other: Vector3D) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }
}

impl Div<f64> for Vector3D {
    type Output = Vector3D;
    fn div(self, other: f64) -> Vector3D {
        Vector3D::new(self.x / other, self.y / other, self.z / other)
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
        assert_eq!(c.x, 8.1);
        assert_eq!(c.y, -5.0);
        assert_eq!(c.z, 12.1);
    }

    #[test]
    fn sub() {
        let a = Vector3D::new(2.0, 3.5, 4.8);
        let b = Vector3D::new(6.1, -8.5, 7.3);

        let c = a - b;
        assert_eq!(c.x, -4.1);
        assert_eq!(c.y, 12.0);
        assert_eq!(c.z, -2.5);
    }

    #[test]
    fn mul() {
        let a = Vector3D::new(2.0, 3.5, 4.8);
        let b = 2.0;

        let c = b * a;
        assert_eq!(c.x, 4.0);
        assert_eq!(c.y, 7.0);
        assert_eq!(c.z, 9.6);

        let b = 1.5;
        let c = a * b;
        assert_eq!(c.x, 3.0);
        assert_eq!(c.y, 5.25);
        assert_eq!(c.z, 7.199999999999999);
    }

    #[test]
    fn dot_product() {
        let a = Vector3D::new(2.1, 3.5, 4.8);
        let b = Vector3D::new(6.1, -8.5, 7.3);

        let c = a * b;
        assert_eq!(c, 18.1);
    }
}
