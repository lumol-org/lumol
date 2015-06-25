/*
 * Cymbalum, Molecular Simulation in Rust
 * Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
*/
use std::f64::consts::PI;
use std::ops::Mul;
use std::convert::Into;
use ::types::*;

/// The type of a cell determine how we will be able to compute the periodic
/// boundaries condition.
#[derive(Clone, Copy)]
pub enum CellType {
    /// Infinite unit cell, with no boundaries
    INFINITE,
    /// Orthorombic unit cell, with cuboide shape
    ORTHOROMBIC,
    /// Triclinic unit cell, with arbitrary parallelepipedic shape
    TRICLINIC,
}

/// The Universe type hold all the data about a system.
#[derive(Clone, Copy)]
pub struct UnitCell {
    data: Matrix3,
    celltype: CellType,
}

impl UnitCell {
    /// Create an infinite unit cell
    pub fn new() -> UnitCell {
        UnitCell{data: Matrix3::zero(), celltype: CellType::INFINITE}
    }
    /// Create an orthorombic unit cell
    pub fn ortho(a: f64, b: f64, c: f64) -> UnitCell {
        UnitCell{data: Matrix3::new(a, 0.0, 0.0,
                                    0.0, b, 0.0,
                                    0.0, 0.0, c),
                 celltype: CellType::ORTHOROMBIC}
    }
    /// Create a cubic unit cell
    pub fn cubic(L: f64) -> UnitCell {
        UnitCell{data: Matrix3::new(L, 0.0, 0.0,
                                    0.0, L, 0.0,
                                    0.0, 0.0, L),
                 celltype: CellType::ORTHOROMBIC}
    }
    /// Create a triclinic unit cell
    pub fn triclinic(a: f64, b: f64, c: f64, alpha: f64, beta: f64, gamma: f64) -> UnitCell {
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
                 celltype: CellType::TRICLINIC}
    }

    /// Get the cell type
    pub fn celltype(&self) -> CellType {
        self.celltype
    }
    /// Set the cell type
    pub fn set_celltype(&mut self, ctype: CellType) {
        self.celltype = ctype;
    }

    /// Get the first vector of the cell
    pub fn vect_a(&self) -> Vector3D {
        let x = self.data[(0, 0)];
        let y = self.data[(1, 0)];
        let z = self.data[(2, 0)];
        Vector3D::new(x, y, z)
    }
    /// Get the first length of the cell
    pub fn a(&self) -> f64 {
        match self.celltype {
            CellType::TRICLINIC => self.vect_a().norm(),
            _ => self.data[(0, 0)],
        }
    }

    /// Get the second length of the cell
    pub fn vect_b(&self) -> Vector3D {
        let x = self.data[(0, 1)];
        let y = self.data[(1, 1)];
        let z = self.data[(2, 1)];
        Vector3D::new(x, y, z)
    }
    /// Get the second length of the cell
    pub fn b(&self) -> f64 {
        match self.celltype {
            CellType::TRICLINIC => self.vect_b().norm(),
            _ => self.data[(1, 1)],
        }
    }

    /// Get the third length of the cell
    pub fn vect_c(&self) -> Vector3D {
        let x = self.data[(0, 2)];
        let y = self.data[(1, 2)];
        let z = self.data[(2, 2)];
        Vector3D::new(x, y, z)
    }
    /// Get the second length of the cell
    pub fn c(&self) -> f64 {
        match self.celltype {
            CellType::TRICLINIC => self.vect_c().norm(),
            _ => self.data[(2, 2)],
        }
    }

    /// Get the first angle of the cell
    pub fn alpha(&self) -> f64 {
        match self.celltype {
            CellType::TRICLINIC => {
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
            CellType::TRICLINIC => {
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
            CellType::TRICLINIC => {
                let a = self.vect_a();
                let b = self.vect_b();
                rad2deg(angle(a, b))
            },
            _ => 90.0,
        }
    }

    /// Get the volume angle of the cell
    pub fn volume(&self) -> f64 {
        match self.celltype {
            CellType::INFINITE => 0.0,
            CellType::ORTHOROMBIC => self.a()*self.b()*self.c(),
            CellType::TRICLINIC => {
                // The volume is the mixed product of the three cell vectors
                let a = self.vect_a();
                let b = self.vect_b();
                let c = self.vect_c();
                a * (b ^ c)
            },
        }
    }

    /// Scale this unit cell in-place by multiplying the H matrix by `s`.
    pub fn scale_mut(&mut self, s: Matrix3) {
        self.data = s * self.data;
    }

    /// Scale this unit cell by multiplying the H matrix by `s`, and return a new unit cell
    pub fn scale(&self, s: Matrix3) -> UnitCell {
        UnitCell{data: s * self.data, celltype: self.celltype}
    }

    /// Periodic boundary conditions distance between the point `u` and the point `v`
    pub fn distance(&self, u: &Vector3D, v: &Vector3D) -> f64 {
        let d = *v - *u;
        match self.celltype {
            CellType::INFINITE => d.norm(),
            CellType::ORTHOROMBIC => {
                let x = d.x - f64::round(d.x/self.a())*self.a();
                let y = d.y - f64::round(d.y/self.b())*self.b();
                let z = d.z - f64::round(d.z/self.c())*self.c();
                return f64::sqrt(x*x + y*y + z*z);
            },
            CellType::TRICLINIC => {
                let inv = self.data.inverse();
                let s = inv * d;
                let x = s.x - f64::round(s.x);
                let y = s.y - f64::round(s.y);
                let z = s.z - f64::round(s.z);
                let res = Vector3D::new(x, y, z);
                return (self.data*res).norm();
            },
        }
    }
}

/// Convert `x` from degrees to radians
fn deg2rad(x: f64) -> f64 {
    x * PI / 180.0
}
/// Convert `x` from radians to degrees
fn rad2deg(x: f64) -> f64 {
    x * 180.0 / PI
}

/// Get the angles between the vectors `u` and `v`.
fn angle(u: Vector3D, v: Vector3D) -> f64 {
    let un = u.normalize();
    let vn = v.normalize();
    f64::acos(un*vn)
}

#[cfg(test)]
mod tests {
    use super::*;
    use ::types::*;

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
}
