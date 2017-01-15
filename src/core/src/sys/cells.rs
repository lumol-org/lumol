// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

//! Simulations in computational chemistry are often made using periodic
//! boundaries conditions. The `UnitCell` type represents the enclosing box of
//! a simulated system, with some type of periodic condition.
use std::f64::consts::PI;

use types::{Matrix3, Vector3D, Zero};

/// The shape of a cell determine how we will be able to compute the periodic
/// boundaries condition.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CellShape {
    /// Infinite unit cell, with no boundaries
    Infinite,
    /// Orthorombic unit cell, with cuboid shape
    Orthorombic,
    /// Triclinic unit cell, with arbitrary parallelepipedic shape
    Triclinic,
}

/// An UnitCell defines the system physical boundaries.
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

impl Default for UnitCell {
    fn default() -> UnitCell {
        UnitCell::new()
    }
}

impl UnitCell {
    /// Create an infinite unit cell
    pub fn new() -> UnitCell {
        UnitCell{
            cell: Matrix3::zero(),
            inv: Matrix3::zero(),
            shape: CellShape::Infinite
        }
    }
    /// Create an orthorhombic unit cell, with side lengths `a, b, c`.
    pub fn ortho(a: f64, b: f64, c: f64) -> UnitCell {
        assert!(a > 0.0 && b > 0.0 && c > 0.0, "Cell lengths must be positive");
        let cell = Matrix3::new(a, 0.0, 0.0,
                                0.0, b, 0.0,
                                0.0, 0.0, c);
        UnitCell{
            cell: cell,
            inv: cell.inverse(),
            shape: CellShape::Orthorombic
        }
    }
    /// Create a cubic unit cell, with side lengths `length, length, length`.
    pub fn cubic(length: f64) -> UnitCell {
        assert!(length > 0.0, "Cell lengths must be positive");
        let cell = Matrix3::new(length,  0.0  ,  0.0  ,
                                  0.0 , length,  0.0  ,
                                  0.0 ,  0.0  , length);
        UnitCell{
            cell: cell,
            inv: cell.inverse(),
            shape: CellShape::Orthorombic
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
        let c_y = c * (cos_alpha - cos_beta*cos_gamma)/sin_gamma;
        let c_z = f64::sqrt(c*c - c_y*c_y - c_x*c_x);

        let cell = Matrix3::new(a,   b_x, c_x,
                                0.0, b_y, c_y,
                                0.0, 0.0, c_z);

        UnitCell{
            cell: cell,
            inv: cell.inverse(),
            shape: CellShape::Triclinic
        }
    }

    /// Get the cell shape
    #[inline] pub fn shape(&self) -> CellShape {
        self.shape
    }

    /// Get the first vector of the cell
    pub fn vect_a(&self) -> Vector3D {
        let x = self.cell[(0, 0)];
        let y = self.cell[(1, 0)];
        let z = self.cell[(2, 0)];
        Vector3D::new(x, y, z)
    }
    /// Get the first length of the cell (i.e. the norm of the first vector of
    /// the cell)
    pub fn a(&self) -> f64 {
        match self.shape {
            CellShape::Triclinic => self.vect_a().norm(),
            CellShape::Orthorombic | CellShape::Infinite => self.cell[(0, 0)],
        }
    }

    /// Get the second vector of the cell
    pub fn vect_b(&self) -> Vector3D {
        let x = self.cell[(0, 1)];
        let y = self.cell[(1, 1)];
        let z = self.cell[(2, 1)];
        Vector3D::new(x, y, z)
    }

    /// Get the second length of the cell (i.e. the norm of the second vector of
    /// the cell)
    pub fn b(&self) -> f64 {
        match self.shape {
            CellShape::Triclinic => self.vect_b().norm(),
            CellShape::Orthorombic | CellShape::Infinite => self.cell[(1, 1)],
        }
    }

    /// Get the third vector of the cell
    pub fn vect_c(&self) -> Vector3D {
        let x = self.cell[(0, 2)];
        let y = self.cell[(1, 2)];
        let z = self.cell[(2, 2)];
        Vector3D::new(x, y, z)
    }

    /// Get the third length of the cell (i.e. the norm of the third vector of
    /// the cell)
    pub fn c(&self) -> f64 {
        match self.shape {
            CellShape::Triclinic => self.vect_c().norm(),
            CellShape::Orthorombic | CellShape::Infinite => self.cell[(2, 2)],
        }
    }

    /// Get the distances between faces of the unit cell
    pub fn lengths(&self) -> [f64; 3] {
        assert!(self.shape != CellShape::Infinite);

        let (a, b, c) = (self.vect_a(), self.vect_b(), self.vect_c());
        // Plans normal vectors
        let na = (b ^ c).normalized();
        let nb = (c ^ a).normalized();
        let nc = (a ^ b).normalized();

        let x = f64::abs(na[0]*a[0] + na[1]*a[1] + na[2]*a[2]);
        let y = f64::abs(nb[0]*b[0] + nb[1]*b[1] + nb[2]*b[2]);
        let z = f64::abs(nc[0]*c[0] + nc[1]*c[1] + nc[2]*c[2]);

        [x, y, z]
    }

    /// Get the first angle of the cell
    pub fn alpha(&self) -> f64 {
        match self.shape {
            CellShape::Triclinic => {
                let b = self.vect_b();
                let c = self.vect_c();
                angle(b, c).to_degrees()
            },
            CellShape::Orthorombic | CellShape::Infinite => 90.0,
        }
    }

    /// Get the second angle of the cell
    pub fn beta(&self) -> f64 {
        match self.shape {
            CellShape::Triclinic => {
                let a = self.vect_a();
                let c = self.vect_c();
                angle(a, c).to_degrees()
            },
            CellShape::Orthorombic | CellShape::Infinite => 90.0,
        }
    }

    /// Get the third angle of the cell
    pub fn gamma(&self) -> f64 {
        match self.shape {
            CellShape::Triclinic => {
                let a = self.vect_a();
                let b = self.vect_b();
                angle(a, b).to_degrees()
            },
            CellShape::Orthorombic | CellShape::Infinite => 90.0,
        }
    }

    /// Get the volume angle of the cell
    pub fn volume(&self) -> f64 {
        let volume = match self.shape {
            CellShape::Infinite => 0.0,
            CellShape::Orthorombic => self.a()*self.b()*self.c(),
            CellShape::Triclinic => {
                // The volume is the mixed product of the three cell vectors
                let a = self.vect_a();
                let b = self.vect_b();
                let c = self.vect_c();
                a * (b ^ c)
            },
        };
        assert!(volume >= 0.0, "Volume is not positive!");
        return volume;
    }

    /// Scale this unit cell in-place by multiplying the cell matrix by `factor`.
    #[inline] pub fn scale_mut(&mut self, factor: Matrix3) {
        self.cell *= factor;
        self.inv = self.cell.inverse();
    }

    /// Scale this unit cell by multiplying the cell matrix by `s`, and return a
    /// new scaled unit cell
    #[inline] pub fn scale(&self, s: Matrix3) -> UnitCell {
        let cell = s * self.cell;
        UnitCell{cell: cell, inv: cell.inverse(), shape: self.shape}
    }

    /// Get the reciprocal vectors of this unit cell
    pub fn reciprocal_vectors(&self) -> (Vector3D, Vector3D, Vector3D) {
        let volume = self.volume();
        let (a, b, c) = (self.vect_a(), self.vect_b(), self.vect_c());

        let rec_a = (2.0 * PI / volume) * (b ^ c);
        let rec_b = (2.0 * PI / volume) * (c ^ a);
        let rec_c = (2.0 * PI / volume) * (a ^ b);

        return (rec_a, rec_b, rec_c);
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
            CellShape::Orthorombic => {
                vect[0] -= f64::floor(vect[0] / self.a()) * self.a();
                vect[1] -= f64::floor(vect[1] / self.b()) * self.b();
                vect[2] -= f64::floor(vect[2] / self.c()) * self.c();
            },
            CellShape::Triclinic => {
                let mut fractional = self.fractional(vect);
                fractional[0] -= f64::floor(fractional[0]);
                fractional[1] -= f64::floor(fractional[1]);
                fractional[2] -= f64::floor(fractional[2]);
                *vect = self.cartesian(&fractional);
            },
        }
    }

    /// Find the image of a vector in the unit cell, obeying the periodic
    /// boundary conditions. For a cubic cell of side length `L`, this produce a
    /// vector with all components in `[-L/2, L/2)`.
    pub fn vector_image(&self, vect: &mut Vector3D) {
        match self.shape {
            CellShape::Infinite => (),
            CellShape::Orthorombic => {
                vect[0] -= f64::round(vect[0] / self.a()) * self.a();
                vect[1] -= f64::round(vect[1] / self.b()) * self.b();
                vect[2] -= f64::round(vect[2] / self.c()) * self.c();
            },
            CellShape::Triclinic => {
                let mut fractional = self.fractional(vect);
                fractional[0] -= f64::round(fractional[0]);
                fractional[1] -= f64::round(fractional[1]);
                fractional[2] -= f64::round(fractional[2]);
                *vect = self.cartesian(&fractional);
            },
        }
    }

    /// Get the fractional representation of the `v` vector in this cell
    #[inline]
    pub fn fractional(&self, vect: &Vector3D) -> Vector3D {
        return self.inv * vect;
    }

    /// Get the Cartesian representation of the fractional `v` vector in this
    /// cell
    #[inline]
    pub fn cartesian(&self, frac: &Vector3D) -> Vector3D {
        return self.cell * frac;
    }

    /// Periodic boundary conditions distance between the point `u` and the point `v`
    pub fn distance(&self, u: &Vector3D, v: &Vector3D) -> f64 {
        let mut d = *v - *u;
        self.vector_image(&mut d);
        return d.norm();
    }

    /// Get the angle formed by the points at `a`, `b` and `c` using periodic
    /// boundary conditions.
    pub fn angle(&self, a: &Vector3D, b: &Vector3D, c: &Vector3D) -> f64 {
        let mut x = *a - *b;
        self.vector_image(&mut x);
        let mut y = *c - *b;
        self.vector_image(&mut y);

        let xn = x.normalized();
        let yn = y.normalized();
        return f64::acos(xn * yn);
    }

    /// Get the angle formed by the points at `a`, `b` and `c` using periodic
    /// boundary conditions and its derivatives.
    pub fn angle_and_derivatives(&self, a: &Vector3D, b: &Vector3D, c: &Vector3D) -> (f64, Vector3D, Vector3D, Vector3D) {
        let mut x = *a - *b;
        self.vector_image(&mut x);
        let mut y = *c - *b;
        self.vector_image(&mut y);

        let x_norm = x.norm();
        let y_norm = y.norm();
        let xn = x.normalized();
        let yn = y.normalized();

        let cos = xn * yn;
        let sin_inv = 1.0 / f64::sqrt(1.0 - cos*cos);

        let d1 = sin_inv * (cos * xn - yn) / x_norm;
        let d3 = sin_inv * (cos * yn - xn) / y_norm;
        let d2 = - (d1 + d3);

        return (f64::acos(cos), d1, d2, d3);
    }


    /// Get the dihedral angle formed by the points at `a`, `b`, `c`, `d` using
    /// periodic boundary conditions.
    pub fn dihedral(&self, a: &Vector3D, b: &Vector3D, c: &Vector3D, d: &Vector3D) -> f64 {
        let mut r12 = *b - *a;
        self.vector_image(&mut r12);
        let mut r23 = *c - *b;
        self.vector_image(&mut r23);
        let mut r34 = *d - *c;
        self.vector_image(&mut r34);

        let u = r12 ^ r23;
        let v = r23 ^ r34;
        return f64::atan2(r23.norm() * v * r12, u * v);
    }

    /// Get the dihedral angle and and its derivatives defined by the points at
    /// `a`, `b`, `c` and `d` using periodic boundary conditions.
    pub fn dihedral_and_derivatives(&self, a: &Vector3D, b: &Vector3D, c: &Vector3D, d: &Vector3D)
    -> (f64, Vector3D, Vector3D, Vector3D , Vector3D) {
        let mut r12 = *b - *a;
        self.vector_image(&mut r12);
        let mut r23 = *c - *b;
        self.vector_image(&mut r23);
        let mut r34 = *d - *c;
        self.vector_image(&mut r34);

        let x = r12 ^ r23;
        let y = r23 ^ r34;
        let norm2_x = x.norm2();
        let norm2_y = y.norm2();
        let norm2_r23 = r23.norm2();
        let norm_r23 = f64::sqrt(norm2_r23);

        let d1 = (-norm_r23 / norm2_x) * x;
        let d4 = (norm_r23 / norm2_y) * y;

        let r23r34 = r23 * r34;
        let r12r23 = r12 * r23;

        let d2 = (-r12r23 / norm2_r23 - 1.0) * d1 + (r23r34 / norm2_r23) * d4;
        let d3 = (-r23r34 / norm2_r23 - 1.0) * d4 + (r12r23 / norm2_r23) * d1;

        let phi = f64::atan2(r23.norm() * y * r12, x * y);
        return (phi, d1, d2, d3, d4)
    }
}

/// Get the angles between the vectors `u` and `v`.
fn angle(u: Vector3D, v: Vector3D) -> f64 {
    let un = u.normalized();
    let vn = v.normalized();
    f64::acos(un*vn)
}

#[cfg(test)]
mod tests {
    use super::*;
    use types::*;
    use std::f64;
    use std::f64::consts::PI;

    #[test]
    #[should_panic]
    fn negative_cubic() {
        let _ = UnitCell::cubic(-4.0);
    }

    #[test]
    #[should_panic]
    fn negative_ortho() {
        let _ = UnitCell::ortho(3.0, 0.0, -5.0);
    }

    #[test]
    #[should_panic]
    fn negative_triclinic() {
        let _ = UnitCell::triclinic(3.0, 0.0, -5.0, 90.0, 90.0, 90.0);
    }

    #[test]
    fn infinite() {
        let cell = UnitCell::new();
        assert_eq!(cell.shape(), CellShape::Infinite);

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
        assert_eq!(cell.shape(), CellShape::Orthorombic);

        assert_eq!(cell.vect_a(), Vector3D::new(3.0, 0.0, 0.0));
        assert_eq!(cell.vect_b(), Vector3D::new(0.0, 3.0, 0.0));
        assert_eq!(cell.vect_c(), Vector3D::new(0.0, 0.0, 3.0));

        assert_eq!(cell.a(), 3.0);
        assert_eq!(cell.b(), 3.0);
        assert_eq!(cell.c(), 3.0);

        assert_eq!(cell.alpha(), 90.0);
        assert_eq!(cell.beta(), 90.0);
        assert_eq!(cell.gamma(), 90.0);

        assert_eq!(cell.volume(), 3.0*3.0*3.0);
    }

    #[test]
    fn orthorombic() {
        let cell = UnitCell::ortho(3.0, 4.0, 5.0);
        assert_eq!(cell.shape(), CellShape::Orthorombic);

        assert_eq!(cell.vect_a(), Vector3D::new(3.0, 0.0, 0.0));
        assert_eq!(cell.vect_b(), Vector3D::new(0.0, 4.0, 0.0));
        assert_eq!(cell.vect_c(), Vector3D::new(0.0, 0.0, 5.0));

        assert_eq!(cell.a(), 3.0);
        assert_eq!(cell.b(), 4.0);
        assert_eq!(cell.c(), 5.0);

        assert_eq!(cell.alpha(), 90.0);
        assert_eq!(cell.beta(), 90.0);
        assert_eq!(cell.gamma(), 90.0);

        assert_eq!(cell.volume(), 3.0*4.0*5.0);
    }

    #[test]
    fn triclinic() {
        let cell = UnitCell::triclinic(3.0, 4.0, 5.0, 80.0, 90.0, 110.0);
        assert_eq!(cell.shape(), CellShape::Triclinic);

        assert_eq!(cell.vect_a(), Vector3D::new(3.0, 0.0, 0.0));
        assert_eq!(cell.vect_b()[2], 0.0);

        assert_eq!(cell.a(), 3.0);
        assert_eq!(cell.b(), 4.0);
        assert_eq!(cell.c(), 5.0);

        assert_eq!(cell.alpha(), 80.0);
        assert_eq!(cell.beta(), 90.0);
        assert_eq!(cell.gamma(), 110.0);

        assert_approx_eq!(cell.volume(), 55.410529, 1e-6);
    }

    #[test]
    fn lengths() {
        let ortho = UnitCell::ortho(3.0, 4.0, 5.0);
        assert_eq!(ortho.lengths(), [3.0, 4.0, 5.0]);

        let triclinic = UnitCell::triclinic(3.0, 4.0, 5.0, 90.0, 90.0, 90.0);
        assert_eq!(triclinic.lengths(), [3.0, 4.0, 5.0]);

        let triclinic = UnitCell::triclinic(3.0, 4.0, 5.0, 90.0, 80.0, 100.0);
        assert_eq!(triclinic.lengths(), [2.908132319388713, 3.9373265973230853, 4.921658246653857]);
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
    fn scale_mut() {
        let mut cell = UnitCell::ortho(3.0, 4.0, 5.0);
        cell.scale_mut(2.0 * Matrix3::one());

        assert_eq!(cell.a(), 6.0);
        assert_eq!(cell.b(), 8.0);
        assert_eq!(cell.c(), 10.0);
    }

    #[test]
    fn reciprocal_vectors() {
        let cell = UnitCell::ortho(3.0, 4.0, 5.0);
        let v = 3.0 * 4.0 * 5.0;
        let two_pi_vol = 2.0 * PI / v;
        let (rec_a, rec_b, rec_c) = cell.reciprocal_vectors();

        assert_eq!(rec_a, Vector3D::new(4.0 * 5.0 * two_pi_vol, 0.0, 0.0));
        assert_eq!(rec_b, Vector3D::new(0.0, 3.0 * 5.0 * two_pi_vol, 0.0));
        assert_eq!(rec_c, Vector3D::new(0.0, 0.0, 3.0 * 4.0 * two_pi_vol));

        let cell = UnitCell::triclinic(3.0, 4.0, 5.0, 90.0, 90.0, 90.0);
        let (rec_a, rec_b, rec_c) = cell.reciprocal_vectors();

        let delta_a = rec_a - Vector3D::new(4.0 * 5.0 * two_pi_vol, 0.0, 0.0);
        assert_approx_eq!(delta_a.norm(), 0.0);

        let delta_b = rec_b - Vector3D::new(0.0, 3.0 * 5.0 * two_pi_vol, 0.0);
        assert_approx_eq!(delta_b.norm(), 0.0);

        let delta_c = rec_c - Vector3D::new(0.0, 0.0, 3.0 * 4.0 * two_pi_vol);
        assert_approx_eq!(delta_c.norm(), 0.0);
    }

    #[test]
    fn distances() {
        // Orthorhombic unit cell
        let cell = UnitCell::ortho(3.0, 4.0, 5.0);
        let u = &Vector3D::zero();
        let v = &Vector3D::new(1.0, 2.0, 6.0);
        assert_eq!(cell.distance(u, v), f64::sqrt(6.0));

        // Infinite unit cell
        let cell = UnitCell::new();
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
        let cell = UnitCell::new();
        let mut v = Vector3D::new(1.0, 1.5, 6.0);
        cell.wrap_vector(&mut v);
        assert_eq!(v, Vector3D::new(1.0, 1.5, 6.0));

        // Triclinic unit cell
        let cell = UnitCell::triclinic(3.0, 4.0, 5.0, 90.0, 90.0, 90.0);
        let mut v = Vector3D::new(1.0, 1.5, 6.0);
        cell.wrap_vector(&mut v);
        let res = Vector3D::new(1.0, 1.5, 1.0);
        assert_approx_eq!(v[0], res[0]);
        assert_approx_eq!(v[1], res[1]);
        assert_approx_eq!(v[2], res[2]);
    }

    #[test]
    fn test_get_wrapped_vector() {
        // Cubic unit cell
        let cell = UnitCell::cubic(10.0);
        let v = Vector3D::new(9.0, 18.0, -6.0);
        let wv = cell.get_wrapped_vector(&v);
        assert_eq!(wv, Vector3D::new(9.0, 8.0, 4.0));

        // Orthorhombic unit cell
        let cell = UnitCell::ortho(3.0, 4.0, 5.0);
        let v = Vector3D::new(1.0, 1.5, 6.0);
        let wv = cell.get_wrapped_vector(&v);
        assert_eq!(wv, Vector3D::new(1.0, 1.5, 1.0));

        // Infinite unit cell
        let cell = UnitCell::new();
        let v = Vector3D::new(1.0, 1.5, 6.0);
        let wv = cell.get_wrapped_vector(&v);
        assert_eq!(wv, Vector3D::new(1.0, 1.5, 6.0));

        // Triclinic unit cell
        let cell = UnitCell::triclinic(3.0, 4.0, 5.0, 90.0, 90.0, 90.0);
        let v = Vector3D::new(1.0, 1.5, 6.0);
        let wv = cell.get_wrapped_vector(&v);
        let res = Vector3D::new(1.0, 1.5, 1.0);
        assert_approx_eq!(wv[0], res[0]);
        assert_approx_eq!(wv[1], res[1]);
        assert_approx_eq!(wv[2], res[2]);
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
        let cell = UnitCell::new();
        let mut v = Vector3D::new(1.0, 1.5, 6.0);
        cell.vector_image(&mut v);
        assert_eq!(v, Vector3D::new(1.0, 1.5, 6.0));

        // Triclinic unit cell
        let cell = UnitCell::triclinic(3.0, 4.0, 5.0, 90.0, 90.0, 90.0);
        let mut v = Vector3D::new(1.0, 1.5, 6.0);
        cell.vector_image(&mut v);
        let res = Vector3D::new(1.0, 1.5, 1.0);
        assert_approx_eq!(v[0], res[0]);
        assert_approx_eq!(v[1], res[1]);
        assert_approx_eq!(v[2], res[2]);
    }

    #[test]
    fn fractional_cartesian() {
        let cell = UnitCell::cubic(5.0);

        assert_eq!(cell.fractional(&Vector3D::new(0.0, 10.0, 4.0)), Vector3D::new(0.0, 2.0, 0.8));
        assert_eq!(cell.cartesian(&Vector3D::new(0.0, 2.0, 0.8)), Vector3D::new(0.0, 10.0, 4.0));

        let cell = UnitCell::triclinic(5.0, 6.0, 3.6, 90.0, 53.0, 77.0);
        let tests = vec![Vector3D::new(0.0, 10.0, 4.0), Vector3D::new(-5.0, 12.0, 4.9)];

        for test in &tests {
            let transformed = cell.cartesian(&cell.fractional(test));
            for i in 0..3 {
                assert_approx_eq!(test[i], transformed[i], 1e-12);
            }
        }
    }

    #[test]
    fn angles() {
        let cell = UnitCell::new();

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
        let cell = UnitCell::new();
        let a = Vector3D::new(0.0, 0.02, 0.0);
        let b = Vector3D::new(-0.784729, -0.5548997, 0.0);
        let c = Vector3D::new(0.784729, -0.5548997, 0.0);

        let (angle, d1, d2, d3) = cell.angle_and_derivatives(&a, &b, &c);
        assert_eq!(angle, cell.angle(&a, &b, &c));

        // Check by comparison to finite differences
        for i in 0..3 {
            let mut p = a;
            p[i] += EPS;
            assert_approx_eq!((cell.angle(&p, &b, &c) - angle)/EPS, d1[i], EPS);
        }

        for i in 0..3 {
            let mut p = b;
            p[i] += EPS;
            assert_approx_eq!((cell.angle(&a, &p, &c) - angle)/EPS, d2[i], EPS);
        }

        for i in 0..3 {
            let mut p = c;
            p[i] += EPS;
            assert_approx_eq!((cell.angle(&a, &b, &p) - angle)/EPS, d3[i], EPS);
        }
    }

    #[test]
    fn dihedrals() {
        let cell = UnitCell::new();

        let a = Vector3D::zero();
        let b = Vector3D::new(1.0, 0.0, 0.0);
        let c = Vector3D::new(1.0, 1.0, 0.0);
        let d = Vector3D::new(2.0, 1.0, 0.0);
        assert_eq!(cell.dihedral(&a, &b, &c, &d), PI);

        let a = Vector3D::new(1.241, 0.444, 0.349);
        let b = Vector3D::new(-0.011, -0.441, 0.333);
        let c = Vector3D::new(-1.176, 0.296, -0.332);
        let d = Vector3D::new(-1.396, 1.211, 0.219);
        assert_approx_eq!(cell.dihedral(&a, &b, &c, &d), -1.0453789, 1e-6);
    }

    #[test]
    fn dihedral_derivatives() {
        const EPS: f64 = 1e-6;
        let cell = UnitCell::new();
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
            assert_approx_eq!((cell.dihedral(&p, &b, &c, &d) - angle)/EPS, d1[i], EPS);
        }

        for i in 0..3 {
            let mut p = b;
            p[i] += EPS;
            assert_approx_eq!((cell.dihedral(&a, &p, &c, &d) - angle)/EPS, d2[i], EPS);
        }

        for i in 0..3 {
            let mut p = c;
            p[i] += EPS;
            assert_approx_eq!((cell.dihedral(&a, &b, &p, &d) - angle)/EPS, d3[i], EPS);
        }

        for i in 0..3 {
            let mut p = d;
            p[i] += EPS;
            assert_approx_eq!((cell.dihedral(&a, &b, &c, &p) - angle)/EPS, d4[i], EPS);
        }
    }
}
