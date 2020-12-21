// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

//! Complex type
use std::f64;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use num_traits::{One, Zero};

#[derive(Debug, Clone, Copy, PartialEq, Default)]
/// Complex number, with double precision real and imaginary parts.
///
/// `Complex` implements all the usual arithmetic operations:
///
/// ```
/// # use lumol_core::types::Complex;
///
/// let w = Complex::cartesian(-1.0, 0.5);
/// let z = Complex::cartesian(4.0, 2.0);
///
/// // Addition
/// let c = w + z;
/// assert_eq!(c, Complex::cartesian(3.0, 2.5));
///
/// // Subtraction
/// let c = w - z;
/// assert_eq!(c, Complex::cartesian(-5.0, -1.5));
///
/// // Multiplication
/// let c = w * z;
/// assert_eq!(c, Complex::cartesian(-5.0, 0.0));
///
/// let c = 42.0 * w;
/// assert_eq!(c, Complex::cartesian(-42.0, 21.0));
///
/// // Division
/// let c = z / 2.0;
/// assert_eq!(c, Complex::cartesian(2.0, 1.0));
/// ```
pub struct Complex {
    /// Real part of the complex
    real: f64,
    /// Imaginary part of the complex
    imag: f64,
}

impl Complex {
    /// Create a new `Complex` from a norm `r` and a phase `phi` in radians.
    ///
    /// # Examples
    ///
    /// ```
    /// # use lumol_core::types::Complex;
    /// # use std::f64;
    /// let z = Complex::polar(3.0, f64::consts::PI);
    /// assert_eq!(z.norm(), 3.0);
    /// ```
    pub fn polar(r: f64, phi: f64) -> Complex {
        Complex {
            real: r * f64::cos(phi),
            imag: r * f64::sin(phi),
        }
    }

    /// Create a complex from Cartesian coordinates
    ///
    /// # Examples
    ///
    /// ```
    /// # use lumol_core::types::Complex;
    /// let z = Complex::cartesian(3.0, -2.0);
    /// assert_eq!(z.real(), 3.0);
    /// assert_eq!(z.imag(), -2.0);
    /// ```
    pub fn cartesian(x: f64, y: f64) -> Complex {
        Complex { real: x, imag: y }
    }

    /// Create a new `Complex` with both cartesian coordinate set to 0.
    ///
    /// # Examples
    ///
    /// ```
    /// # use lumol_core::types::Complex;
    /// let z = Complex::zero();
    /// assert_eq!(z.norm(), 0.0);
    /// ```
    pub fn zero() -> Complex {
        <Complex as Zero>::zero()
    }

    /// Get the real part of the complex
    /// # Examples
    /// ```
    /// # use lumol_core::types::Complex;
    /// let z = Complex::cartesian(3.0, -2.0);
    /// assert_eq!(z.real(), 3.0);
    /// ```
    #[inline]
    pub fn real(&self) -> f64 {
        self.real
    }

    /// Get the imaginary part of the complex
    /// # Examples
    /// ```
    /// # use lumol_core::types::Complex;
    /// let z = Complex::cartesian(3.0, -2.0);
    /// assert_eq!(z.imag(), -2.0);
    /// ```
    #[inline]
    pub fn imag(&self) -> f64 {
        self.imag
    }

    /// Get the phase of the complex in the [-π, π) interval
    /// # Examples
    /// ```
    /// # use lumol_core::types::Complex;
    /// let z = Complex::polar(2.0, 0.3);
    /// assert_eq!(z.phase(), 0.3);
    /// ```
    #[inline]
    pub fn phase(&self) -> f64 {
        f64::atan2(self.imag, self.real)
    }

    /// Get the norm of the complex
    /// # Examples
    /// ```
    /// # use lumol_core::types::Complex;
    /// # use std::f64;
    /// let z = Complex::polar(2.0, 0.3);
    /// assert_eq!(z.norm(), 2.0);
    ///
    /// let z = Complex::cartesian(2.0, 1.0);
    /// assert_eq!(z.norm(), f64::sqrt(5.0));
    /// ```
    #[inline]
    pub fn norm(&self) -> f64 {
        f64::sqrt(self.norm2())
    }

    /// Get the square of the norm if this complex
    /// # Examples
    /// ```
    /// # use lumol_core::types::Complex;
    /// let z = Complex::cartesian(2.0, 1.0);
    /// assert_eq!(z.norm2(), 5.0);
    /// ```
    #[inline]
    pub fn norm2(&self) -> f64 {
        self.real * self.real + self.imag * self.imag
    }

    /// Get the conjugate of the complex
    /// # Examples
    /// ```
    /// # use lumol_core::types::Complex;
    /// let z = Complex::cartesian(2.0, 1.0);
    /// assert_eq!(z.conj(), Complex::cartesian(2.0, -1.0));
    /// ```
    #[inline]
    pub fn conj(&self) -> Complex {
        Complex {
            real: self.real,
            imag: -self.imag,
        }
    }
}

impl Add<Complex> for Complex {
    type Output = Complex;

    fn add(self, other: Complex) -> Complex {
        Complex {
            real: self.real + other.real,
            imag: self.imag + other.imag,
        }
    }
}

impl AddAssign<Complex> for Complex {
    fn add_assign(&mut self, other: Complex) {
        self.real += other.real;
        self.imag += other.imag;
    }
}

impl Sub<Complex> for Complex {
    type Output = Complex;

    fn sub(self, other: Complex) -> Complex {
        Complex {
            real: self.real - other.real,
            imag: self.imag - other.imag,
        }
    }
}

impl SubAssign<Complex> for Complex {
    fn sub_assign(&mut self, other: Complex) {
        self.real -= other.real;
        self.imag -= other.imag;
    }
}

impl Neg for Complex {
    type Output = Complex;

    fn neg(self) -> Complex {
        Complex {
            real: -self.real,
            imag: -self.imag,
        }
    }
}

impl Mul<Complex> for Complex {
    type Output = Complex;

    fn mul(self, other: Complex) -> Complex {
        Complex {
            real: self.real * other.real - self.imag * other.imag,
            imag: self.real * other.imag + self.imag * other.real,
        }
    }
}

impl Mul<f64> for Complex {
    type Output = Complex;

    fn mul(self, other: f64) -> Complex {
        Complex {
            real: self.real * other,
            imag: self.imag * other,
        }
    }
}

impl Mul<Complex> for f64 {
    type Output = Complex;

    fn mul(self, other: Complex) -> Complex {
        Complex {
            real: self * other.real,
            imag: self * other.imag,
        }
    }
}

impl MulAssign<Complex> for Complex {
    fn mul_assign(&mut self, other: Complex) {
        let real = self.real * other.real - self.imag * other.imag;
        let imag = self.real * other.imag + self.imag * other.real;

        self.real = real;
        self.imag = imag;
    }
}

impl MulAssign<f64> for Complex {
    fn mul_assign(&mut self, other: f64) {
        self.real *= other;
        self.imag *= other;
    }
}

impl Div<Complex> for Complex {
    type Output = Complex;

    fn div(self, other: Complex) -> Complex {
        let norm = other.norm2();
        let real = self.real * other.real + self.imag * other.imag;
        let imag = -self.real * other.imag + self.imag * other.real;

        Complex {
            real: real / norm,
            imag: imag / norm,
        }
    }
}

impl Div<f64> for Complex {
    type Output = Complex;

    fn div(self, other: f64) -> Complex {
        Complex {
            real: self.real / other,
            imag: self.imag / other,
        }
    }
}

impl DivAssign<Complex> for Complex {
    fn div_assign(&mut self, other: Complex) {
        let norm = other.norm2();
        let real = self.real * other.real + self.imag * other.imag;
        let imag = -self.real * other.imag + self.imag * other.real;

        self.real = real / norm;
        self.imag = imag / norm;
    }
}

impl DivAssign<f64> for Complex {
    fn div_assign(&mut self, other: f64) {
        self.real /= other;
        self.imag /= other;
    }
}

impl Zero for Complex {
    fn zero() -> Complex {
        Complex::cartesian(0.0, 0.0)
    }

    fn is_zero(&self) -> bool {
        self.norm2() == 0.0
    }
}

impl One for Complex {
    fn one() -> Complex {
        Complex::cartesian(1.0, 0.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts;
    use approx::assert_ulps_eq;

    #[test]
    fn norm() {
        let c = Complex::polar(3.0, 5.0);
        assert_ulps_eq!(c.norm(), 3.0);

        let c = Complex::polar(-3.0, 0.0);
        assert_ulps_eq!(c.norm(), 3.0);
    }

    #[test]
    fn phase() {
        for &phase in &[-consts::PI, -3.1, -1.5, 0.0, 0.1, 2.0, 3.1] {
            let c = Complex::polar(1.0, phase);
            assert_ulps_eq!(c.phase(), phase);
        }
        let c = Complex::polar(1.0, -8.0);
        assert_ulps_eq!(c.phase(), -8.0 + 2.0 * consts::PI);

        let c = Complex::polar(1.0, 12.0);
        assert_ulps_eq!(c.phase(), 12.0 - 4.0 * consts::PI, max_ulps = 10);

        let c = Complex::polar(1.0, consts::PI);
        assert_eq!(c.phase(), consts::PI);
        let c = Complex::polar(1.0, -consts::PI);
        assert_eq!(c.phase(), -consts::PI);

        let c = Complex::polar(3.0, 0.0);
        assert_eq!(c.phase(), 0.0);

        let c = Complex::polar(-3.0, 0.0);
        assert_eq!(c.phase(), -consts::PI);
    }

    #[test]
    fn conj() {
        let c = Complex::cartesian(3.0, 5.0);
        assert_eq!(c.conj(), Complex::cartesian(3.0, -5.0));

        let c = Complex::cartesian(3.0, -5.0);
        assert_eq!(c.conj(), Complex::cartesian(3.0, 5.0));

        let c = Complex::cartesian(3.0, 0.0);
        assert_eq!(c.conj(), Complex::cartesian(3.0, 0.0));
    }

    #[test]
    fn cartesian() {
        let c = Complex::polar(1.0, 0.0);
        assert_ulps_eq!(c.real(), 1.0);
        assert_ulps_eq!(c.imag(), 0.0);

        let c = Complex::polar(1.0, consts::PI);
        assert_ulps_eq!(c.real(), -1.0);
        assert_ulps_eq!(c.imag(), 0.0);

        let c = Complex::polar(1.0, consts::FRAC_PI_2);
        assert_ulps_eq!(c.real(), 0.0);
        assert_ulps_eq!(c.imag(), 1.0);

        let c = Complex::polar(1.0, consts::FRAC_PI_4);
        assert_ulps_eq!(c.real(), consts::FRAC_1_SQRT_2);
        assert_ulps_eq!(c.imag(), consts::FRAC_1_SQRT_2);

        let c = Complex::cartesian(consts::FRAC_1_SQRT_2, consts::FRAC_1_SQRT_2);
        assert_ulps_eq!(c.norm(), 1.0);
        assert_ulps_eq!(c.phase(), consts::FRAC_PI_4);
    }

    #[test]
    fn add() {
        let mut a = Complex::polar(2.0, 0.2);
        let b = Complex::polar(1.0, 0.5);
        let c = a + b;

        assert_eq!(c.real(), a.real() + b.real());
        assert_eq!(c.imag(), a.imag() + b.imag());

        a += b;
        assert_eq!(a, c);
    }

    #[test]
    fn sub() {
        let mut a = Complex::polar(2.0, 0.2);
        let b = Complex::polar(1.0, 0.5);
        let c = a - b;

        assert_eq!(c.real(), a.real() - b.real());
        assert_eq!(c.imag(), a.imag() - b.imag());

        a -= b;
        assert_eq!(a, c);
    }

    #[test]
    fn neg() {
        let a = Complex::polar(2.0, 0.2);
        let c = -a;

        assert_ulps_eq!(c.norm(), a.norm());
        assert_ulps_eq!(c.phase(), a.phase() - consts::PI);

        assert_ulps_eq!(c.real(), -a.real());
        assert_ulps_eq!(c.imag(), -a.imag());
    }

    #[test]
    fn mul() {
        let mut a = Complex::polar(2.0, 1.4);
        let b = Complex::polar(1.0, 0.5);
        let c = a * b;

        assert_eq!(c.norm(), a.norm() * b.norm());
        assert_eq!(c.phase(), a.phase() + b.phase());

        let c = 3.0 * a;
        assert_eq!(c.norm(), 3.0 * a.norm());
        assert_eq!(c.phase(), a.phase());

        let c = a * 3.0;
        assert_eq!(c.norm(), 3.0 * a.norm());
        assert_eq!(c.phase(), a.phase());

        let c = -2.0 * a;
        assert_eq!(c.norm(), 2.0 * a.norm());
        assert_ulps_eq!(c.phase(), a.phase() - consts::PI);

        let c = a * b;
        a *= b;
        assert_eq!(a, c);

        let c = 3.0 * a;
        a *= 3.0;
        assert_eq!(a, c);
    }

    #[test]
    fn div() {
        let mut a = Complex::polar(2.0, 0.2);
        let b = Complex::polar(1.0, 0.5);
        let c = a / b;

        assert_eq!(c.norm(), a.norm() / b.norm());
        assert_eq!(c.phase(), a.phase() - b.phase());

        let c = a / 3.0;
        assert_eq!(c.norm(), a.norm() / 3.0);
        assert_ulps_eq!(c.phase(), a.phase());

        let c = a / (-2.0);
        assert_eq!(c.norm(), a.norm() / 2.0);
        assert_ulps_eq!(c.phase(), a.phase() - consts::PI);

        let c = a / b;
        a /= b;
        assert_eq!(a, c);

        let c = a / 3.0;
        a /= 3.0;
        assert_eq!(a, c);
    }
}
