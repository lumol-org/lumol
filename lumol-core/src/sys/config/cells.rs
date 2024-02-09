// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

//! Simulations in computational chemistry are often made using periodic
//! boundaries conditions. The `UnitCell` type represents the enclosing box of
//! a simulated system, with some type of periodic condition.
use std::f64;
use std::f64::consts::PI;

use crate::{Matrix3, Vector3D};

/// The shape of a cell determine how we will be able to compute the periodic
/// boundaries condition.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CellShape {
    /// Infinite unit cell, with no boundaries
    Infinite,
    /// Orthorhombic unit cell, with cuboid shape
    Orthorhombic,
    /// Triclinic unit cell, with arbitrary parallelepipedic shape
    Triclinic,
}

/// An `UnitCell` defines the system physical boundaries.
///
/// The shape of the cell can be any of the [`CellShape`][CellShape], and will
/// influence how periodic boundary conditions are applied.
///
/// [CellShape]: enum.CellShape.html
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct UnitCell {
    /// Unit cell matrix
    cell: Matrix3,
    /// Inverse of the unit cell matrix. This is cached for performance reason,
    /// and MUST be updated as needed.
    inv: Matrix3,
    /// Unit cell shape
    shape: CellShape,
}

impl UnitCell {
    /// Create an infinite unit cell
    pub fn infinite() -> UnitCell {
        UnitCell {
            cell: Matrix3::zero(),
            inv: Matrix3::zero(),
            shape: CellShape::Infinite,
        }
    }
    /// Create an orthorhombic unit cell, with side lengths `a, b, c`.
    pub fn ortho(a: f64, b: f64, c: f64) -> UnitCell {
        assert!(a > 0.0 && b > 0.0 && c > 0.0, "Cell lengths must be positive");
        let cell = Matrix3::new([[a, 0.0, 0.0], [0.0, b, 0.0], [0.0, 0.0, c]]);
        UnitCell {
            cell: cell,
            inv: cell.inverse(),
            shape: CellShape::Orthorhombic,
        }
    }
    /// Create a cubic unit cell, with side lengths `length, length, length`.
    pub fn cubic(length: f64) -> UnitCell {
        assert!(length > 0.0, "Cell lengths must be positive");
        let cell = Matrix3::new([[length, 0.0, 0.0], [0.0, length, 0.0], [0.0, 0.0, length]]);
        UnitCell {
            cell: cell,
            inv: cell.inverse(),
            shape: CellShape::Orthorhombic,
        }
    }
    /// Create a triclinic unit cell, with side lengths `a, b, c` and angles
    /// `alpha, beta, gamma`.
    pub fn triclinic(a: f64, b: f64, c: f64, alpha: f64, beta: f64, gamma: f64) -> UnitCell {
        assert!(a > 0.0 && b > 0.0 && c > 0.0, "Cell lengths must be positive");
        let cos_alpha = alpha.to_radians().cos();
        let cos_beta = beta.to_radians().cos();
        let (sin_gamma, cos_gamma) = gamma.to_radians().sin_cos();

        let b_x = b * cos_gamma;
        let b_y = b * sin_gamma;

        let c_x = c * cos_beta;
        let c_y = c * (cos_alpha - cos_beta * cos_gamma) / sin_gamma;
        let c_z = f64::sqrt(c * c - c_y * c_y - c_x * c_x);

        let cell = Matrix3::new([[a, b_x, c_x], [0.0, b_y, c_y], [0.0, 0.0, c_z]]);

        UnitCell {
            cell: cell,
            inv: cell.inverse(),
            shape: CellShape::Triclinic,
        }
    }

    /// Get the cell shape
    #[inline]
    pub fn shape(&self) -> CellShape {
        self.shape
    }

    /// Check if this unit cell is infinite, *i.e.* if it does not have
    /// periodic boundary conditions.
    pub fn is_infinite(&self) -> bool {
        self.shape() == CellShape::Infinite
    }

    /// Get the first length of the cell (i.e. the norm of the first vector of
    /// the cell)
    pub fn a(&self) -> f64 {
        match self.shape {
            CellShape::Triclinic => self.vect_a().norm(),
            CellShape::Orthorhombic | CellShape::Infinite => self.cell[0][0],
        }
    }

    /// Get the second length of the cell (i.e. the norm of the second vector of
    /// the cell)
    pub fn b(&self) -> f64 {
        match self.shape {
            CellShape::Triclinic => self.vect_b().norm(),
            CellShape::Orthorhombic | CellShape::Infinite => self.cell[1][1],
        }
    }

    /// Get the third length of the cell (i.e. the norm of the third vector of
    /// the cell)
    pub fn c(&self) -> f64 {
        match self.shape {
            CellShape::Triclinic => self.vect_c().norm(),
            CellShape::Orthorhombic | CellShape::Infinite => self.cell[2][2],
        }
    }

    /// Get the distances between faces of the unit cell
    pub fn lengths(&self) -> Vector3D {
        if self.shape == CellShape::Infinite {
            return Vector3D::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
        }

        let (a, b, c) = (self.vect_a(), self.vect_b(), self.vect_c());
        // Plans normal vectors
        let na = (b ^ c).normalized();
        let nb = (c ^ a).normalized();
        let nc = (a ^ b).normalized();

        Vector3D::new(f64::abs(na * a), f64::abs(nb * b), f64::abs(nc * c))
    }

    /// Get the first angle of the cell
    pub fn alpha(&self) -> f64 {
        match self.shape {
            CellShape::Triclinic => {
                let b = self.vect_b();
                let c = self.vect_c();
                angle(b, c).to_degrees()
            }
            CellShape::Orthorhombic | CellShape::Infinite => 90.0,
        }
    }

    /// Get the second angle of the cell
    pub fn beta(&self) -> f64 {
        match self.shape {
            CellShape::Triclinic => {
                let a = self.vect_a();
                let c = self.vect_c();
                angle(a, c).to_degrees()
            }
            CellShape::Orthorhombic | CellShape::Infinite => 90.0,
        }
    }

    /// Get the third angle of the cell
    pub fn gamma(&self) -> f64 {
        match self.shape {
            CellShape::Triclinic => {
                let a = self.vect_a();
                let b = self.vect_b();
                angle(a, b).to_degrees()
            }
            CellShape::Orthorhombic | CellShape::Infinite => 90.0,
        }
    }

    /// Get the volume of the cell
    pub fn volume(&self) -> f64 {
        let volume = match self.shape {
            CellShape::Infinite => 0.0,
            CellShape::Orthorhombic => self.a() * self.b() * self.c(),
            CellShape::Triclinic => {
                // The volume is the mixed product of the three cell vectors
                let a = self.vect_a();
                let b = self.vect_b();
                let c = self.vect_c();
                a * (b ^ c)
            }
        };
        assert!(volume >= 0.0, "Volume is not positive!");
        return volume;
    }

    /// Scale this unit cell in-place by multiplying the cell matrix by `factor`.
    #[inline]
    pub fn scale_mut(&mut self, factor: Matrix3) {
        assert!(self.shape() != CellShape::Infinite, "can not scale infinite cells");
        self.cell *= factor;
        self.inv = self.cell.inverse();
    }

    /// Scale this unit cell by multiplying the cell matrix by `s`, and return a
    /// new scaled unit cell
    #[inline]
    pub fn scale(&self, s: Matrix3) -> UnitCell {
        assert!(self.shape() != CellShape::Infinite, "can not scale infinite cells");
        let cell = s * self.cell;
        UnitCell {
            cell: cell,
            inv: cell.inverse(),
            shape: self.shape,
        }
    }

    /// Get the reciprocal vector with the given `index`. This vector is null
    /// for infinite cells.
    pub fn k_vector(&self, index: [f64; 3]) -> Vector3D {
        return 2.0 * PI * self.inv * Vector3D::from(index);
    }

    /// Get the matricial representation of the unit cell
    pub fn matrix(&self) -> Matrix3 {
        self.cell
    }

    /// Get the first vector of the cell
    fn vect_a(&self) -> Vector3D {
        let x = self.cell[0][0];
        let y = self.cell[1][0];
        let z = self.cell[2][0];
        Vector3D::new(x, y, z)
    }

    /// Get the second vector of the cell
    fn vect_b(&self) -> Vector3D {
        let x = self.cell[0][1];
        let y = self.cell[1][1];
        let z = self.cell[2][1];
        Vector3D::new(x, y, z)
    }

    /// Get the third vector of the cell
    fn vect_c(&self) -> Vector3D {
        let x = self.cell[0][2];
        let y = self.cell[1][2];
        let z = self.cell[2][2];
        Vector3D::new(x, y, z)
    }
}

/// Geometric operations using periodic boundary conditions
impl UnitCell {
    /// Wrap a vector in the unit cell, obeying the periodic boundary conditions.
    /// For a cubic cell of side length `L`, this produce a vector with all
    /// components in `[0, L)`.
    pub fn wrap_vector(&self, vect: &mut Vector3D) {
        match self.shape {
            CellShape::Infinite => (),
            CellShape::Orthorhombic => {
                vect[0] -= f64::floor(vect[0] / self.a()) * self.a();
                vect[1] -= f64::floor(vect[1] / self.b()) * self.b();
                vect[2] -= f64::floor(vect[2] / self.c()) * self.c();
            }
            CellShape::Triclinic => {
                let mut fractional = self.fractional(vect);
                fractional[0] -= f64::floor(fractional[0]);
                fractional[1] -= f64::floor(fractional[1]);
                fractional[2] -= f64::floor(fractional[2]);
                *vect = self.cartesian(&fractional);
            }
        }
    }

    /// Find the image of a vector in the unit cell, obeying the periodic
    /// boundary conditions. For a cubic cell of side length `L`, this produce a
    /// vector with all components in `[-L/2, L/2)`.
    pub fn vector_image(&self, vect: &mut Vector3D) {
        match self.shape {
            CellShape::Infinite => (),
            CellShape::Orthorhombic => {
                vect[0] -= f64::round(vect[0] / self.a()) * self.a();
                vect[1] -= f64::round(vect[1] / self.b()) * self.b();
                vect[2] -= f64::round(vect[2] / self.c()) * self.c();
            }
            CellShape::Triclinic => {
                let mut fractional = self.fractional(vect);
                fractional[0] -= f64::round(fractional[0]);
                fractional[1] -= f64::round(fractional[1]);
                fractional[2] -= f64::round(fractional[2]);
                *vect = self.cartesian(&fractional);
            }
        }
    }

    /// Get the fractional representation of the `vector` in this cell
    #[inline]
    pub fn fractional(&self, vector: &Vector3D) -> Vector3D {
        return self.inv * vector;
    }

    /// Get the Cartesian representation of the `fractional` vector in this
    /// cell
    #[inline]
    pub fn cartesian(&self, fractional: &Vector3D) -> Vector3D {
        return self.cell * fractional;
    }

    /// Periodic boundary conditions distance between the point `u` and the point `v`
    pub fn distance(&self, u: &Vector3D, v: &Vector3D) -> f64 {
        let mut d = v - u;
        self.vector_image(&mut d);
        return d.norm();
    }

    /// Get the angle formed by the points at `r1`, `r2` and `r3` using periodic
    /// boundary conditions.
    pub fn angle(&self, r1: &Vector3D, r2: &Vector3D, r3: &Vector3D) -> f64 {
        let mut r12 = r1 - r2;
        self.vector_image(&mut r12);
        let mut r23 = r3 - r2;
        self.vector_image(&mut r23);

        return f64::acos(r12 * r23 / (r12.norm() * r23.norm()));
    }

    /// Get the angle formed by the points at `r1`, `r2` and `r3` using periodic
    /// boundary conditions and its derivatives.
    pub fn angle_and_derivatives(
        &self,
        r1: &Vector3D,
        r2: &Vector3D,
        r3: &Vector3D,
    ) -> (f64, Vector3D, Vector3D, Vector3D) {
        let mut r12 = r1 - r2;
        self.vector_image(&mut r12);
        let mut r23 = r3 - r2;
        self.vector_image(&mut r23);

        let r12_norm = r12.norm();
        let r23_norm = r23.norm();
        let r12n = r12 / r12_norm;
        let r23n = r23 / r23_norm;

        let cos = r12n * r23n;
        let sin_inv = 1.0 / f64::sqrt(1.0 - cos * cos);

        let d1 = sin_inv * (cos * r12n - r23n) / r12_norm;
        let d3 = sin_inv * (cos * r23n - r12n) / r23_norm;
        let d2 = -(d1 + d3);

        return (f64::acos(cos), d1, d2, d3);
    }


    /// Get the dihedral angle formed by the points at `r1`, `r2`, `r3`, and `r4` using
    /// periodic boundary conditions.
    pub fn dihedral(&self, r1: &Vector3D, r2: &Vector3D, r3: &Vector3D, r4: &Vector3D) -> f64 {
        let mut r12 = r2 - r1;
        self.vector_image(&mut r12);
        let mut r23 = r3 - r2;
        self.vector_image(&mut r23);
        let mut r34 = r4 - r3;
        self.vector_image(&mut r34);

        let u = r12 ^ r23;
        let v = r23 ^ r34;
        return f64::atan2(r23.norm() * v * r12, u * v);
    }

    /// Get the dihedral angle and and its derivatives defined by the points at
    /// `r1`, `r2`, `r3`, and `r4` using periodic boundary conditions.
    pub fn dihedral_and_derivatives(
        &self,
        r1: &Vector3D,
        r2: &Vector3D,
        r3: &Vector3D,
        r4: &Vector3D,
    ) -> (f64, Vector3D, Vector3D, Vector3D, Vector3D) {
        let mut r12 = r2 - r1;
        self.vector_image(&mut r12);
        let mut r23 = r3 - r2;
        self.vector_image(&mut r23);
        let mut r34 = r4 - r3;
        self.vector_image(&mut r34);

        let u = r12 ^ r23;
        let v = r23 ^ r34;
        let u_norm2 = u.norm2();
        let v_norm2 = v.norm2();
        let r23_norm2 = r23.norm2();
        let r23_norm = f64::sqrt(r23_norm2);

        let d1 = (-r23_norm / u_norm2) * u;
        let d4 = ( r23_norm / v_norm2) * v;

        let r23_r34 = r23 * r34;
        let r12_r23 = r12 * r23;

        let d2 = (-r12_r23 / r23_norm2 - 1.0) * d1 + (r23_r34 / r23_norm2) * d4;
        let d3 = (-r23_r34 / r23_norm2 - 1.0) * d4 + (r12_r23 / r23_norm2) * d1;

        let phi = f64::atan2(r23_norm * v * r12, u * v);
        return (phi, d1, d2, d3, d4);
    }
}

/// Get the angles between the vectors `u` and `v`.
fn angle(u: Vector3D, v: Vector3D) -> f64 {
    let un = u.normalized();
    let vn = v.normalized();
    f64::acos(un * vn)
}

#[cfg(test)]
#[allow(clippy::unreadable_literal)]
mod tests {
    use super::*;
    use std::f64;
    use std::f64::consts::PI;
    use crate::Matrix3;

    use approx::{assert_ulps_eq, assert_relative_eq};

    #[test]
    #[should_panic(expected="Cell lengths must be positive")]
    fn negative_cubic() {
        let _ = UnitCell::cubic(-4.0);
    }

    #[test]
    #[should_panic(expected="Cell lengths must be positive")]
    fn negative_ortho() {
        let _ = UnitCell::ortho(3.0, 0.0, -5.0);
    }

    #[test]
    #[should_panic(expected="Cell lengths must be positive")]
    fn negative_triclinic() {
        let _ = UnitCell::triclinic(3.0, 0.0, -5.0, 90.0, 90.0, 90.0);
    }

    #[test]
    fn infinite() {
        let cell = UnitCell::infinite();
        assert_eq!(cell.shape(), CellShape::Infinite);
        assert!(cell.is_infinite());

        assert_eq!(cell.vect_a(), Vector3D::zero());
        assert_eq!(cell.vect_b(), Vector3D::zero());
        assert_eq!(cell.vect_c(), Vector3D::zero());

        assert_eq!(cell.a(), 0.0);
        assert_eq!(cell.b(), 0.0);
        assert_eq!(cell.c(), 0.0);

        assert_eq!(cell.alpha(), 90.0);
        assert_eq!(cell.beta(), 90.0);
        assert_eq!(cell.gamma(), 90.0);

        assert_eq!(cell.volume(), 0.0);
    }

    #[test]
    fn cubic() {
        let cell = UnitCell::cubic(3.0);
        assert_eq!(cell.shape(), CellShape::Orthorhombic);
        assert!(!cell.is_infinite());

        assert_eq!(cell.vect_a(), Vector3D::new(3.0, 0.0, 0.0));
        assert_eq!(cell.vect_b(), Vector3D::new(0.0, 3.0, 0.0));
        assert_eq!(cell.vect_c(), Vector3D::new(0.0, 0.0, 3.0));

        assert_eq!(cell.a(), 3.0);
        assert_eq!(cell.b(), 3.0);
        assert_eq!(cell.c(), 3.0);

        assert_eq!(cell.alpha(), 90.0);
        assert_eq!(cell.beta(), 90.0);
        assert_eq!(cell.gamma(), 90.0);

        assert_eq!(cell.volume(), 3.0 * 3.0 * 3.0);
    }

    #[test]
    fn orthorhombic() {
        let cell = UnitCell::ortho(3.0, 4.0, 5.0);
        assert_eq!(cell.shape(), CellShape::Orthorhombic);
        assert!(!cell.is_infinite());

        assert_eq!(cell.vect_a(), Vector3D::new(3.0, 0.0, 0.0));
        assert_eq!(cell.vect_b(), Vector3D::new(0.0, 4.0, 0.0));
        assert_eq!(cell.vect_c(), Vector3D::new(0.0, 0.0, 5.0));

        assert_eq!(cell.a(), 3.0);
        assert_eq!(cell.b(), 4.0);
        assert_eq!(cell.c(), 5.0);

        assert_eq!(cell.alpha(), 90.0);
        assert_eq!(cell.beta(), 90.0);
        assert_eq!(cell.gamma(), 90.0);

        assert_eq!(cell.volume(), 3.0 * 4.0 * 5.0);
    }

    #[test]
    fn triclinic() {
        let cell = UnitCell::triclinic(3.0, 4.0, 5.0, 80.0, 90.0, 110.0);
        assert_eq!(cell.shape(), CellShape::Triclinic);
        assert!(!cell.is_infinite());

        assert_eq!(cell.vect_a(), Vector3D::new(3.0, 0.0, 0.0));
        assert_eq!(cell.vect_b()[2], 0.0);

        assert_eq!(cell.a(), 3.0);
        assert_eq!(cell.b(), 4.0);
        assert_eq!(cell.c(), 5.0);

        assert_eq!(cell.alpha(), 80.0);
        assert_eq!(cell.beta(), 90.0);
        assert_eq!(cell.gamma(), 110.0);

        assert_relative_eq!(cell.volume(), 55.410529, epsilon = 1e-6);
    }

    #[test]
    fn lengths() {
        let ortho = UnitCell::ortho(3.0, 4.0, 5.0);
        assert_eq!(ortho.lengths(), Vector3D::new(3.0, 4.0, 5.0));

        let triclinic = UnitCell::triclinic(3.0, 4.0, 5.0, 90.0, 90.0, 90.0);
        assert_eq!(triclinic.lengths(), Vector3D::new(3.0, 4.0, 5.0));

        let triclinic = UnitCell::triclinic(3.0, 4.0, 5.0, 90.0, 80.0, 100.0);
        assert_eq!(triclinic.lengths(), Vector3D::new(2.908132319388713, 3.9373265973230853, 4.921658246653857));
    }

    #[test]
    fn scale() {
        let cell = UnitCell::ortho(3.0, 4.0, 5.0);
        let cell = cell.scale(2.0 * Matrix3::one());

        assert_eq!(cell.a(), 6.0);
        assert_eq!(cell.b(), 8.0);
        assert_eq!(cell.c(), 10.0);
    }

    #[test]
    #[should_panic(expected="can not scale infinite cells")]
    fn scale_infinite() {
        let cell = UnitCell::infinite();
        let _ = cell.scale(2.0 * Matrix3::one());
    }

    #[test]
    fn scale_mut() {
        let mut cell = UnitCell::ortho(3.0, 4.0, 5.0);
        cell.scale_mut(2.0 * Matrix3::one());

        assert_eq!(cell.a(), 6.0);
        assert_eq!(cell.b(), 8.0);
        assert_eq!(cell.c(), 10.0);
    }

    #[test]
    #[should_panic(expected="can not scale infinite cells")]
    fn scale_mut_infinite() {
        let mut cell = UnitCell::infinite();
        cell.scale_mut(2.0 * Matrix3::one());
    }

    #[test]
    fn k_vectors() {
        let cell = UnitCell::ortho(3.0, 4.0, 5.0);
        let two_pi_vol = 2.0 * PI / cell.volume();
        let kvec = cell.k_vector([1.0, 1.0, 1.0]);

        assert_eq!(kvec, Vector3D::new(
            4.0 * 5.0 * two_pi_vol,
            3.0 * 5.0 * two_pi_vol,
            3.0 * 4.0 * two_pi_vol,
        ));

        let cell = UnitCell::triclinic(3.0, 4.0, 5.0, 90.0, 90.0, 90.0);
        let kvec = cell.k_vector([1.0, 1.0, 1.0]);

        assert_ulps_eq!(kvec, Vector3D::new(
            4.0 * 5.0 * two_pi_vol,
            3.0 * 5.0 * two_pi_vol,
            3.0 * 4.0 * two_pi_vol,
        ));

        let cell = UnitCell::infinite();
        let kvec = cell.k_vector([1.0, 1.0, 1.0]);
        assert_ulps_eq!(kvec, Vector3D::new(0.0, 0.0, 0.0));
    }

    #[test]
    fn distances() {
        // Orthorhombic unit cell
        let cell = UnitCell::ortho(3.0, 4.0, 5.0);
        let u = &Vector3D::zero();
        let v = &Vector3D::new(1.0, 2.0, 6.0);
        assert_eq!(cell.distance(u, v), f64::sqrt(6.0));

        // Infinite unit cell
        let cell = UnitCell::infinite();
        assert_eq!(cell.distance(u, v), v.norm());

        // Triclinic unit cell
        let cell = UnitCell::triclinic(3.0, 4.0, 5.0, 90.0, 90.0, 90.0);
        assert_eq!(cell.distance(u, v), f64::sqrt(6.0));
    }

    #[test]
    fn wrap_vector() {
        // Cubic unit cell
        let cell = UnitCell::cubic(10.0);
        let mut v = Vector3D::new(9.0, 18.0, -6.0);
        cell.wrap_vector(&mut v);
        assert_eq!(v, Vector3D::new(9.0, 8.0, 4.0));

        // Orthorhombic unit cell
        let cell = UnitCell::ortho(3.0, 4.0, 5.0);
        let mut v = Vector3D::new(1.0, 1.5, 6.0);
        cell.wrap_vector(&mut v);
        assert_eq!(v, Vector3D::new(1.0, 1.5, 1.0));

        // Infinite unit cell
        let cell = UnitCell::infinite();
        let mut v = Vector3D::new(1.0, 1.5, 6.0);
        cell.wrap_vector(&mut v);
        assert_eq!(v, Vector3D::new(1.0, 1.5, 6.0));

        // Triclinic unit cell
        let cell = UnitCell::triclinic(3.0, 4.0, 5.0, 90.0, 90.0, 90.0);
        let mut v = Vector3D::new(1.0, 1.5, 6.0);
        cell.wrap_vector(&mut v);
        let res = Vector3D::new(1.0, 1.5, 1.0);
        assert_ulps_eq!(v[0], res[0], max_ulps = 5);
        assert_ulps_eq!(v[1], res[1], max_ulps = 5);
        assert_ulps_eq!(v[2], res[2], max_ulps = 5);
    }

    #[test]
    fn vector_image() {
        // Cubic unit cell
        let cell = UnitCell::cubic(10.0);
        let mut v = Vector3D::new(9.0, 18.0, -6.0);
        cell.vector_image(&mut v);
        assert_eq!(v, Vector3D::new(-1.0, -2.0, 4.0));

        // Orthorhombic unit cell
        let cell = UnitCell::ortho(3.0, 4.0, 5.0);
        let mut v = Vector3D::new(1.0, 1.5, 6.0);
        cell.vector_image(&mut v);
        assert_eq!(v, Vector3D::new(1.0, 1.5, 1.0));

        // Infinite unit cell
        let cell = UnitCell::infinite();
        let mut v = Vector3D::new(1.0, 1.5, 6.0);
        cell.vector_image(&mut v);
        assert_eq!(v, Vector3D::new(1.0, 1.5, 6.0));

        // Triclinic unit cell
        let cell = UnitCell::triclinic(3.0, 4.0, 5.0, 90.0, 90.0, 90.0);
        let mut v = Vector3D::new(1.0, 1.5, 6.0);
        cell.vector_image(&mut v);
        let res = Vector3D::new(1.0, 1.5, 1.0);
        assert_ulps_eq!(v[0], res[0], max_ulps = 5);
        assert_ulps_eq!(v[1], res[1], max_ulps = 5);
        assert_ulps_eq!(v[2], res[2], max_ulps = 5);
    }

    #[test]
    fn fractional_cartesian() {
        let cell = UnitCell::cubic(5.0);

        assert_eq!(cell.fractional(&Vector3D::new(0.0, 10.0, 4.0)), Vector3D::new(0.0, 2.0, 0.8));
        assert_eq!(cell.cartesian(&Vector3D::new(0.0, 2.0, 0.8)), Vector3D::new(0.0, 10.0, 4.0));

        let cell = UnitCell::triclinic(5.0, 6.0, 3.6, 90.0, 53.0, 77.0);
        let tests = vec![
            Vector3D::new(0.0, 10.0, 4.0),
            Vector3D::new(-5.0, 12.0, 4.9),
        ];

        for test in &tests {
            let transformed = cell.cartesian(&cell.fractional(test));
            assert_ulps_eq!(test, &transformed, epsilon = 1e-15);
        }
    }

    #[test]
    fn angles() {
        let cell = UnitCell::infinite();

        let a = Vector3D::new(1.0, 0.0, 0.0);
        let b = Vector3D::zero();
        let c = Vector3D::new(0.0, 1.0, 0.0);
        assert_eq!(cell.angle(&a, &b, &c), PI / 2.0);

        let a = Vector3D::new(1.0, 0.0, 0.0);
        let b = Vector3D::zero();
        let c = Vector3D::new(f64::cos(1.877), f64::sin(1.877), 0.0);
        assert_eq!(cell.angle(&a, &b, &c), 1.877);
    }

    #[test]
    fn angle_derivatives() {
        const EPS: f64 = 1e-6;
        let cell = UnitCell::infinite();
        let a = Vector3D::new(0.0, 0.02, 0.0);
        let b = Vector3D::new(-0.784729, -0.5548997, 0.0);
        let c = Vector3D::new(0.784729, -0.5548997, 0.0);

        let (angle, d1, d2, d3) = cell.angle_and_derivatives(&a, &b, &c);
        assert_eq!(angle, cell.angle(&a, &b, &c));

        // Check by comparison to finite differences
        for i in 0..3 {
            let mut p = a;
            p[i] += EPS;
            assert_ulps_eq!((cell.angle(&p, &b, &c) - angle) / EPS, d1[i], epsilon = 1e-6);
        }

        for i in 0..3 {
            let mut p = b;
            p[i] += EPS;
            assert_ulps_eq!((cell.angle(&a, &p, &c) - angle) / EPS, d2[i], epsilon = 1e-6);
        }

        for i in 0..3 {
            let mut p = c;
            p[i] += EPS;
            assert_ulps_eq!((cell.angle(&a, &b, &p) - angle) / EPS, d3[i], epsilon = 1e-6);
        }
    }

    #[test]
    fn dihedrals() {
        let cell = UnitCell::infinite();

        let a = Vector3D::zero();
        let b = Vector3D::new(1.0, 0.0, 0.0);
        let c = Vector3D::new(1.0, 1.0, 0.0);
        let d = Vector3D::new(2.0, 1.0, 0.0);
        assert_eq!(cell.dihedral(&a, &b, &c, &d), PI);

        let a = Vector3D::new(1.241, 0.444, 0.349);
        let b = Vector3D::new(-0.011, -0.441, 0.333);
        let c = Vector3D::new(-1.176, 0.296, -0.332);
        let d = Vector3D::new(-1.396, 1.211, 0.219);
        assert_relative_eq!(cell.dihedral(&a, &b, &c, &d), -1.0453789626063168);
    }

    #[test]
    #[allow(clippy::many_single_char_names)]
    fn dihedral_derivatives() {
        const EPS: f64 = 1e-6;
        let cell = UnitCell::infinite();
        let a = Vector3D::new(1.241, 0.444, 0.349);
        let b = Vector3D::new(-0.011, -0.441, 0.333);
        let c = Vector3D::new(-1.176, 0.296, -0.332);
        let d = Vector3D::new(-1.396, 1.211, 0.219);

        let (angle, d1, d2, d3, d4) = cell.dihedral_and_derivatives(&a, &b, &c, &d);
        assert_eq!(angle, cell.dihedral(&a, &b, &c, &d));

        // Check by comparison to finite differences
        for i in 0..3 {
            let mut p = a;
            p[i] += EPS;
            assert_ulps_eq!((cell.dihedral(&p, &b, &c, &d) - angle) / EPS, d1[i], epsilon = 1e-6);
        }

        for i in 0..3 {
            let mut p = b;
            p[i] += EPS;
            assert_ulps_eq!((cell.dihedral(&a, &p, &c, &d) - angle) / EPS, d2[i], epsilon = 1e-6);
        }

        for i in 0..3 {
            let mut p = c;
            p[i] += EPS;
            assert_ulps_eq!((cell.dihedral(&a, &b, &p, &d) - angle) / EPS, d3[i], epsilon = 1e-6);
        }

        for i in 0..3 {
            let mut p = d;
            p[i] += EPS;
            assert_ulps_eq!((cell.dihedral(&a, &b, &c, &p) - angle) / EPS, d4[i], epsilon = 1e-6);
        }
    }
}
