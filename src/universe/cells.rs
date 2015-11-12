/* Cymbalum, Molecular Simulation in Rust - Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
 */
//! Simulations in computational chemistry are often made using periodic
//! boundaries conditions. The `UnitCell` type represents the enclosing box of
//! a simulated system, with some type of periodic condition.
use std::f64::consts::PI;

use types::{Matrix3, Vector3D};

/// The type of a cell determine how we will be able to compute the periodic
/// boundaries condition.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CellType {
    /// Infinite unit cell, with no boundaries
    Infinite,
    /// Orthorombic unit cell, with cuboide shape
    Orthorombic,
    /// Triclinic unit cell, with arbitrary parallelepipedic shape
    Triclinic,
}

/// The Universe type hold all the data about a system.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct UnitCell {
    data: Matrix3,
    celltype: CellType,
}

impl UnitCell {
    /// Create an infinite unit cell
    pub fn new() -> UnitCell {
        UnitCell{data: Matrix3::zero(), celltype: CellType::Infinite}
    }
    /// Create an orthorombic unit cell, with side lengths `a, b, c`.
    pub fn ortho(a: f64, b: f64, c: f64) -> UnitCell {
        assert!(a > 0.0 && b > 0.0 && c > 0.0, "Cell lengths must be positive");
        UnitCell{data: Matrix3::new(a, 0.0, 0.0,
                                    0.0, b, 0.0,
                                    0.0, 0.0, c),
                 celltype: CellType::Orthorombic}
    }
    /// Create a cubic unit cell, with side lengths `L, L, L`.
    pub fn cubic(L: f64) -> UnitCell {
        assert!(L > 0.0, "Cell lengths must be positive");
        UnitCell{data: Matrix3::new(L, 0.0, 0.0,
                                    0.0, L, 0.0,
                                    0.0, 0.0, L),
                 celltype: CellType::Orthorombic}
    }
    /// Create a triclinic unit cell, with side lengths `a, b, c` and angles
    /// `alpha, beta, gamma`.
    pub fn triclinic(a: f64, b: f64, c: f64, alpha: f64, beta: f64, gamma: f64) -> UnitCell {
        assert!(a > 0.0 && b > 0.0 && c > 0.0, "Cell lengths must be positive");
        let cos_alpha = deg2rad(alpha).cos();
        let cos_beta = deg2rad(beta).cos();
        let (sin_gamma, cos_gamma) = deg2rad(gamma).sin_cos();

        let b_x = b * cos_gamma;
        let b_y = b * sin_gamma;

        let c_x = c * cos_beta;
        let c_y = c * (cos_alpha - cos_beta*cos_gamma)/sin_gamma;
        let c_z = f64::sqrt(c*c - c_y*c_y - c_x*c_x);

        UnitCell{data: Matrix3::new(a,   b_x, c_x,
                                    0.0, b_y, c_y,
                                    0.0, 0.0, c_z),
                 celltype: CellType::Triclinic}
    }

    /// Get the cell type
    #[inline] pub fn celltype(&self) -> CellType {
        self.celltype
    }
    /// Set the cell type to `ctype`
    #[inline] pub fn set_celltype(&mut self, ctype: CellType) {
        self.celltype = ctype;
    }

    /// Get the first vector of the cell
    pub fn vect_a(&self) -> Vector3D {
        let x = self.data[(0, 0)];
        let y = self.data[(1, 0)];
        let z = self.data[(2, 0)];
        Vector3D::new(x, y, z)
    }
    /// Get the first length of the cell (i.e. the norm of the first vector of
    /// the cell)
    pub fn a(&self) -> f64 {
        match self.celltype {
            CellType::Triclinic => self.vect_a().norm(),
            _ => self.data[(0, 0)],
        }
    }

    /// Get the second vector of the cell
    pub fn vect_b(&self) -> Vector3D {
        let x = self.data[(0, 1)];
        let y = self.data[(1, 1)];
        let z = self.data[(2, 1)];
        Vector3D::new(x, y, z)
    }

    /// Get the second length of the cell (i.e. the norm of the second vector of
    /// the cell)
    pub fn b(&self) -> f64 {
        match self.celltype {
            CellType::Triclinic => self.vect_b().norm(),
            _ => self.data[(1, 1)],
        }
    }

    /// Get the third vector of the cell
    pub fn vect_c(&self) -> Vector3D {
        let x = self.data[(0, 2)];
        let y = self.data[(1, 2)];
        let z = self.data[(2, 2)];
        Vector3D::new(x, y, z)
    }

    /// Get the third length of the cell (i.e. the norm of the third vector of
    /// the cell)
    pub fn c(&self) -> f64 {
        match self.celltype {
            CellType::Triclinic => self.vect_c().norm(),
            _ => self.data[(2, 2)],
        }
    }

    /// Get the distances between faces of the unit cell
    pub fn lengths(&self) -> (f64, f64, f64) {
        assert!(self.celltype != CellType::Infinite);

        let (a, b, c) = (self.vect_a(), self.vect_b(), self.vect_c());
        // Plans normal vectors
        let na = (b ^ c).normalized();
        let nb = (c ^ a).normalized();
        let nc = (a ^ b).normalized();

        let x = f64::abs(na[0]*a[0] + na[1]*a[1] + na[2]*a[2]);
        let y = f64::abs(nb[0]*b[0] + nb[1]*b[1] + nb[2]*b[2]);
        let z = f64::abs(nc[0]*c[0] + nc[1]*c[1] + nc[2]*c[2]);

        return (x, y, z);
    }

    /// Get the first angle of the cell
    pub fn alpha(&self) -> f64 {
        match self.celltype {
            CellType::Triclinic => {
                let b = self.vect_b();
                let c = self.vect_c();
                rad2deg(angle(b, c))
            },
            _ => 90.0,
        }
    }

    /// Get the second angle of the cell
    pub fn beta(&self) -> f64 {
        match self.celltype {
            CellType::Triclinic => {
                let a = self.vect_a();
                let c = self.vect_c();
                rad2deg(angle(a, c))
            },
            _ => 90.0,
        }
    }

    /// Get the third angle of the cell
    pub fn gamma(&self) -> f64 {
        match self.celltype {
            CellType::Triclinic => {
                let a = self.vect_a();
                let b = self.vect_b();
                rad2deg(angle(a, b))
            },
            _ => 90.0,
        }
    }

    /// Get the volume angle of the cell
    pub fn volume(&self) -> f64 {
        let V = match self.celltype {
            CellType::Infinite => 0.0,
            CellType::Orthorombic => self.a()*self.b()*self.c(),
            CellType::Triclinic => {
                // The volume is the mixed product of the three cell vectors
                let a = self.vect_a();
                let b = self.vect_b();
                let c = self.vect_c();
                a * (b ^ c)
            },
        };
        assert!(V >= 0.0, "Volume is not positive!");
        return V;
    }

    /// Scale this unit cell in-place by multiplying the cell matrix by `s`.
    #[inline] pub fn scale_mut(&mut self, s: Matrix3) {
        self.data = s * self.data;
    }

    /// Scale this unit cell by multiplying the cell matrix by `s`, and return a
    /// new scaled unit cell
    #[inline] pub fn scale(&self, s: Matrix3) -> UnitCell {
        UnitCell{data: s * self.data, celltype: self.celltype}
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
    /// Wrap a vector in the unit cell, obeying the periodic boundary conditions
    pub fn wrap_vector(&self, vect: &mut Vector3D) {
        match self.celltype {
            CellType::Infinite => (),
            CellType::Orthorombic => {
                vect.x = vect.x - f64::round(vect.x/self.a())*self.a();
                vect.y = vect.y - f64::round(vect.y/self.b())*self.b();
                vect.z = vect.z - f64::round(vect.z/self.c())*self.c();
            },
            CellType::Triclinic => {
                let mut fractional = self.fractional(vect);
                fractional.x = fractional.x - f64::round(fractional.x);
                fractional.y = fractional.y - f64::round(fractional.y);
                fractional.z = fractional.z - f64::round(fractional.z);
                *vect = self.cartesian(&fractional);
            },
        }
    }

    /// Get the fractional representation of the `v` vector in this cell
    #[inline]
    pub fn fractional(&self, v: &Vector3D) -> Vector3D {
        let inv = self.data.inverse();
        let vect = v.clone();
        return inv * vect;
    }

    /// Get the cartesian representation of the fractional `v` vector in this cell
    #[inline]
    pub fn cartesian(&self, v: &Vector3D) -> Vector3D {
        let frac = v.clone();
        return self.data * frac;
    }

    /// Periodic boundary conditions distance between the point `u` and the point `v`
    pub fn distance(&self, u: &Vector3D, v: &Vector3D) -> f64 {
        let mut d = *v - *u;
        self.wrap_vector(&mut d);
        return d.norm();
    }

    /// Get the angle formed by the points at `a`, `b` and `c` using periodic
    /// boundary conditions.
    pub fn angle(&self, a: &Vector3D, b: &Vector3D, c: &Vector3D) -> f64 {
        let mut x = *a - *b;
        self.wrap_vector(&mut x);
        let mut y = *c - *b;
        self.wrap_vector(&mut y);

        let xn = x.normalized();
        let yn = y.normalized();
        return f64::acos(xn * yn);
    }

    /// Get the angle formed by the points at `a`, `b` and `c` using periodic
    /// boundary conditions and its derivatives.
    pub fn angle_and_derivatives(&self, a: &Vector3D, b: &Vector3D, c: &Vector3D) -> (f64, Vector3D, Vector3D, Vector3D) {
        let mut x = *a - *b;
        self.wrap_vector(&mut x);
        let mut y = *c - *b;
        self.wrap_vector(&mut y);

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
        self.wrap_vector(&mut r12);
        let mut r23 = *c - *b;
        self.wrap_vector(&mut r23);
        let mut r34 = *d - *c;
        self.wrap_vector(&mut r34);

        let u = r12 ^ r23;
        let v = r23 ^ r34;
        return f64::atan2(r23.norm() * v * r12, u * v);
    }

    /// Get the dihedral angle formed by the points at `u`, `v` and `w` using
    /// periodic boundary conditions and its derivatives.
    pub fn dihedral_and_derivatives(&self, a: &Vector3D, b: &Vector3D, c: &Vector3D, d: &Vector3D)
    -> (f64, Vector3D, Vector3D, Vector3D , Vector3D) {
        let mut r12 = *b - *a;
        self.wrap_vector(&mut r12);
        let mut r23 = *c - *b;
        self.wrap_vector(&mut r23);
        let mut r34 = *d - *c;
        self.wrap_vector(&mut r34);

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

/// Convert `x` from degrees to radians
#[inline] fn deg2rad(x: f64) -> f64 {
    x * PI / 180.0
}
/// Convert `x` from radians to degrees
#[inline] fn rad2deg(x: f64) -> f64 {
    x * 180.0 / PI
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
        UnitCell::cubic(-4.0);
    }

    #[test]
    #[should_panic]
    fn negative_ortho() {
        UnitCell::ortho(3.0, 0.0, -5.0);
    }

    #[test]
    #[should_panic]
    fn negative_triclinic() {
        UnitCell::triclinic(3.0, 0.0, -5.0, 90.0, 90.0, 90.0);
    }

    #[test]
    fn infinite() {
        let cell = UnitCell::new();
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
        assert_eq!(ortho.lengths(), (3.0, 4.0, 5.0));

        let triclinic = UnitCell::triclinic(3.0, 4.0, 5.0, 90.0, 90.0, 90.0);
        assert_eq!(triclinic.lengths(), (3.0, 4.0, 5.0));

        let triclinic = UnitCell::triclinic(3.0, 4.0, 5.0, 90.0, 80.0, 100.0);
        assert_eq!(triclinic.lengths(), (2.908132319388713, 3.9373265973230853, 4.921658246653857));
    }

    #[test]
    fn scale() {
        let cell = UnitCell::ortho(3.0, 4.0, 5.0);
        let A = 2.0 * Matrix3::one();
        let cell = cell.scale(A);

        assert_eq!(cell.a(), 6.0);
        assert_eq!(cell.b(), 8.0);
        assert_eq!(cell.c(), 10.0);
    }

    #[test]
    fn scale_mut() {
        let mut cell = UnitCell::ortho(3.0, 4.0, 5.0);
        let A = 2.0 * Matrix3::one();
        cell.scale_mut(A);

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
        // Orthorombic unit cell
        let cell = UnitCell::ortho(3.0, 4.0, 5.0);
        let u = &Vector3D::new(0.0, 0.0, 0.0);
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
        // Orthorombic unit cell
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
        assert_approx_eq!(v.x, res.x);
        assert_approx_eq!(v.y, res.y);
        assert_approx_eq!(v.z, res.z);
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
        let b = Vector3D::new(0.0, 0.0, 0.0);
        let c = Vector3D::new(0.0, 1.0, 0.0);
        assert_eq!(cell.angle(&a, &b, &c), PI / 2.0);

        let a = Vector3D::new(1.0, 0.0, 0.0);
        let b = Vector3D::new(0.0, 0.0, 0.0);
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
            let mut p = a.clone();
            p[i] += EPS;
            assert_approx_eq!((cell.angle(&p, &b, &c) - angle)/EPS, d1[i], EPS);
        }

        for i in 0..3 {
            let mut p = b.clone();
            p[i] += EPS;
            assert_approx_eq!((cell.angle(&a, &p, &c) - angle)/EPS, d2[i], EPS);
        }

        for i in 0..3 {
            let mut p = c.clone();
            p[i] += EPS;
            assert_approx_eq!((cell.angle(&a, &b, &p) - angle)/EPS, d3[i], EPS);
        }
    }

    #[test]
    fn dihedrals() {
        let cell = UnitCell::new();

        let a = Vector3D::new(0.0, 0.0, 0.0);
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
            let mut p = a.clone();
            p[i] += EPS;
            assert_approx_eq!((cell.dihedral(&p, &b, &c, &d) - angle)/EPS, d1[i], EPS);
        }

        for i in 0..3 {
            let mut p = b.clone();
            p[i] += EPS;
            assert_approx_eq!((cell.dihedral(&a, &p, &c, &d) - angle)/EPS, d2[i], EPS);
        }

        for i in 0..3 {
            let mut p = c.clone();
            p[i] += EPS;
            assert_approx_eq!((cell.dihedral(&a, &b, &p, &d) - angle)/EPS, d3[i], EPS);
        }

        for i in 0..3 {
            let mut p = d.clone();
            p[i] += EPS;
            assert_approx_eq!((cell.dihedral(&a, &b, &c, &p) - angle)/EPS, d4[i], EPS);
        }
    }
}
