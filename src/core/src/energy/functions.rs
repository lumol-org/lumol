// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 Lumol's contributors — BSD license

use energy::Potential;
use energy::{PairPotential, BondPotential, AnglePotential, DihedralPotential};

/// No-op potential.
///
/// The `NullPotential` always returns 0 as energy and force. It should be used
/// to indicate that there is no potential interaction for a given set of
/// particles.
///
/// # Examples
///
/// ```
/// use lumol::energy::Potential;
/// use lumol::energy::NullPotential;
///
/// let potential = NullPotential;
/// assert_eq!(potential.energy(0.1), 0.0);
/// assert_eq!(potential.energy(100000.0), 0.0);
///
/// assert_eq!(potential.force(0.1), 0.0);
/// assert_eq!(potential.force(100000.0), 0.0);
/// ```
#[derive(Clone, Copy)]
pub struct NullPotential;
impl Potential for NullPotential {
    #[inline] fn energy(&self, _: f64) -> f64 {0.0}
    #[inline] fn force(&self, _: f64) -> f64 {0.0}
}

impl PairPotential for NullPotential {
    fn tail_energy(&self, _: f64) -> f64 {0.0}
    fn tail_virial(&self, _: f64) -> f64 {0.0}
}
impl BondPotential for NullPotential {}
impl AnglePotential for NullPotential {}
impl DihedralPotential for NullPotential {}

/// Lennard-Jones potential.
///
/// The following expression of Lennard-Jones potential is used: `V(r) = 4 *
/// epsilon * ((sigma/r)^12 - (sigma/r)^6)` where `sigma` is the Lennard-Jones
/// distance constant, and `epsilon` the Lennard-Jones energetic constant.
///
/// # Examples
///
/// ```
/// use lumol::energy::Potential;
/// use lumol::energy::LennardJones;
///
/// let potential = LennardJones{sigma: 2.0, epsilon: 10.0};
/// assert_eq!(potential.energy(2.0), 0.0);
/// assert_eq!(potential.energy(3.0), -3.203365942785746);
///
/// assert_eq!(potential.force(2.0), 120.0);
/// ```
#[derive(Clone, Copy)]
pub struct LennardJones {
    /// Distance constant of the Lennard-Jones potential
    pub sigma: f64,
    /// Energy constant of the Lennard-Jones potential
    pub epsilon: f64,
}

impl Potential for LennardJones {
    #[inline]
    fn energy(&self, r: f64) -> f64 {
        let s6 = f64::powi(self.sigma/r, 6);
        4.0 * self.epsilon * (f64::powi(s6, 2) - s6)
    }

    #[inline]
    fn force(&self, r: f64) -> f64 {
        let s6 = f64::powi(self.sigma/r, 6);
        -24.0 * self.epsilon * (s6 - 2.0 * f64::powi(s6, 2)) / r
    }
}

impl PairPotential for LennardJones {
    fn tail_energy(&self, cutoff: f64) -> f64 {
        let s3 = self.sigma * self.sigma * self.sigma;
        let rc3 = cutoff * cutoff * cutoff;
        let s9 = s3 * s3 * s3;
        let rc9 = rc3 * rc3 * rc3;
        return 4.0 / 3.0 * self.epsilon * s3 * (1.0 / 3.0 * s9 / rc9 - s3 / rc3);
    }

    fn tail_virial(&self, cutoff: f64) -> f64 {
        let s3 = self.sigma * self.sigma * self.sigma;
        let rc3 = cutoff * cutoff * cutoff;
        let s9 = s3 * s3 * s3;
        let rc9 = rc3 * rc3 * rc3;
        return 8.0 * self.epsilon * s3 * (2.0 / 3.0 * s9 / rc9 - s3 / rc3);
    }
}

/// Harmonic potential.
///
/// The following energy expression is used: `V(x) = 1/2 * k * (x - x0)^2` where
/// `x0` is the distance equilibrium, and `k` the elastic constant.
///
/// # Examples
///
/// ```
/// use lumol::energy::Potential;
/// use lumol::energy::Harmonic;
///
/// let potential = Harmonic{x0: 2.0, k: 100.0};
/// assert_eq!(potential.energy(2.0), 0.0);
/// assert_eq!(potential.energy(3.0), 50.0);
///
/// assert_eq!(potential.force(2.0), 0.0);
/// assert_eq!(potential.force(1.5), 50.0);
/// ```
#[derive(Clone, Copy)]
pub struct Harmonic {
    /// Spring constant
    pub k: f64,
    /// Equilibrium value
    pub x0: f64,
}

impl Potential for Harmonic {
    #[inline]
    fn energy(&self, x: f64) -> f64 {
        let dx = x - self.x0;
        0.5 * self.k * dx * dx
    }

    #[inline]
    fn force(&self, x: f64) -> f64 {
        self.k * (self.x0 - x)
    }
}

impl PairPotential for Harmonic {
    // These two functions should return infinity, as the Harmonic potential
    // does not goes to zero at infinity. We use 0 instead to ignore the tail
    // contribution to the energy/virial.
    fn tail_energy(&self, _: f64) -> f64 {0.0}
    fn tail_virial(&self, _: f64) -> f64 {0.0}
}

impl BondPotential for Harmonic {}
impl AnglePotential for Harmonic {}
impl DihedralPotential for Harmonic {}

/// Cosine harmonic potential.
///
/// The following potential expression is used: `V(x) = 1/2 * k * (cos(x) -
/// cos(x0))^2` where `x0` is the distance equilibrium, and `k` the elastic
/// constant.
///
/// # Examples
///
/// ```
/// use lumol::energy::Potential;
/// use lumol::energy::CosineHarmonic;
///
/// let potential = CosineHarmonic::new(/* k */ 100.0, /* x0 */ 2.0);
/// assert_eq!(potential.energy(2.0), 0.0);
/// assert_eq!(potential.energy(3.0), 16.464942078100552);
///
/// assert_eq!(potential.force(2.0), 0.0);
/// ```
#[derive(Clone, Copy)]
pub struct CosineHarmonic {
    /// Spring constant
    k: f64,
    /// Cosine of the equilibrium value
    cos_x0: f64,
}

impl CosineHarmonic {
    /// Create a new `CosineHarmonic` potentials, with elastic constant of `k`
    /// and equilibrium value of `x0`
    pub fn new(k: f64, x0:f64) -> CosineHarmonic {
        CosineHarmonic{k: k, cos_x0: f64::cos(x0)}
    }
}

impl Potential for CosineHarmonic {
    #[inline]
    fn energy(&self, x: f64) -> f64 {
        let dr = f64::cos(x) - self.cos_x0;
        0.5 * self.k * dr * dr
    }

    #[inline]
    fn force(&self, x: f64) -> f64 {
        self.k * (f64::cos(x) - self.cos_x0) * f64::sin(x)
    }
}

impl AnglePotential for CosineHarmonic {}
impl DihedralPotential for CosineHarmonic {}

/// Torsion potential.
///
/// This potential is intended for use with dihedral angles, using a custom
/// periodicity and multiple minima.
///
/// The following potential expression is used: `V(x) = k * (1 + cos(n * x -
/// delta))` where `k` is the force constant, `n` the periodicity of the
/// potential, and `delta` the equilibrium angle.
///
/// # Examples
///
/// ```
/// use lumol::energy::Potential;
/// use lumol::energy::Torsion;
/// use std::f64::consts::PI;
///
/// let potential = Torsion{delta: PI / 2.0, k: 10.0, n: 3};
/// assert_eq!(potential.energy(PI / 2.0), 0.0);
/// assert_eq!(potential.energy(PI / 3.0), 10.0);
///
/// assert!(potential.force(PI / 2.0).abs() < 1e-12);
/// ```
#[derive(Clone, Copy)]
pub struct Torsion {
    /// Force constant
    pub k: f64,
    /// Equilibrium value
    pub delta: f64,
    /// Multiplicity of the potential
    pub n: usize,
}

impl Potential for Torsion {
    #[inline]
    fn energy(&self, phi: f64) -> f64 {
        let n = self.n as f64;
        let cos = f64::cos(n * phi - self.delta);
        self.k * (1.0 + cos)
    }

    #[inline]
    fn force(&self, phi: f64) -> f64 {
        let n = self.n as f64;
        let sin = f64::sin(n * phi - self.delta);
        self.k * n * sin
    }
}

impl DihedralPotential for Torsion {}


/// Born-Mayer-Huggins potential.
///
/// The following potential expression is used: `V(x) = A * exp((sigma - r) /
/// rho) - C/r^6 + D/r^8`; where `A`, `C` and `D` are energetic constants;
/// `sigma` and `rho` are length parameters.
///
/// # Examples
///
/// ```
/// use lumol::energy::Potential;
/// use lumol::energy::BornMayerHuggins;
///
/// let potential = BornMayerHuggins{a: 2.0, c: 1.0, d: 0.5, sigma: 1.5, rho: 5.3};
/// assert_eq!(potential.energy(2.2), 1.7446409593340713);
/// assert_eq!(potential.force(2.2), 0.30992873382584607);
/// ```
#[derive(Clone, Copy)]
pub struct BornMayerHuggins {
    /// Exponential term energetic constant
    pub a: f64,
    /// `1/r^6` term energetic constant
    pub c: f64,
    /// `1/r^8` term energetic constant
    pub d: f64,
    /// Sphere diameter length constant
    pub sigma: f64,
    /// Width of the exponential term length constant
    pub rho: f64,
}

impl Potential for BornMayerHuggins {
    #[inline]
    fn energy(&self, r: f64) -> f64 {
        let r2 = r * r;
        let r6 = r2 * r2 * r2;
        let exp = f64::exp((self.sigma - r) / self.rho);
        self.a * exp - self.c / r6 + self.d / (r6 * r2)
    }

    #[inline]
    fn force(&self, r: f64) -> f64 {
        let r2 = r * r;
        let r7 = r2 * r2 * r2 * r;
        let exp = f64::exp((self.sigma - r) / self.rho);
        self.a / self.rho * exp - 6.0 * self.c / r7 + 8.0 * self.d / (r7 * r2)
    }
}

impl PairPotential for BornMayerHuggins {
    fn tail_energy(&self, rc: f64) -> f64 {
        let rc2 = rc * rc;
        let rc3 = rc2 * rc;
        let exp = f64::exp((self.sigma - rc) / self.rho);
        let factor = rc2 - 2.0 * rc * self.rho + 2.0 * self.rho * self.rho;
        self.a * self.rho * exp * factor - self.c / (3.0 * rc3) + self.d / (5.0 * rc2 * rc3)
    }

    fn tail_virial(&self, rc: f64) -> f64 {
        let rc2 = rc * rc;
        let rc3 = rc2 * rc;
        let exp = f64::exp((self.sigma - rc) / self.rho);
        let factor = rc3 + 3.0 * rc2 * self.rho + 6.0 * rc * self.rho * self.rho + 6.0 * self.rho * self.rho * self.rho;
        self.a * exp * factor - 20.0 * self.c / rc3 + 8.0 * self.d / (5.0 * rc2 * rc3)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use energy::{Potential, PairPotential};
    const EPS: f64 = 1e-9;

    #[test]
    fn null() {
        let null = NullPotential;
        assert_eq!(null.energy(2.0), 0.0);
        assert_eq!(null.energy(2.5), 0.0);

        assert_eq!(null.force(2.0), 0.0);
        assert_eq!(null.force(2.5), 0.0);

        assert_eq!(null.tail_energy(1.0), 0.0);
        assert_eq!(null.tail_virial(1.0), 0.0);

        let e0 = null.energy(2.0);
        let e1 = null.energy(2.0 + EPS);
        assert_approx_eq!((e0 - e1)/EPS, null.force(2.0), EPS);
    }

    #[test]
    fn lj() {
        let lj = LennardJones{epsilon: 0.8, sigma: 2.0};
        assert_eq!(lj.energy(2.0), 0.0);
        assert_eq!(lj.energy(2.5), -0.6189584744448002);

        assert_eq!(lj.tail_energy(1.0), 1388.0888888888887);
        assert_eq!(lj.tail_energy(2.0), -5.688888888888889);
        assert_eq!(lj.tail_energy(14.42), -0.022767318648783084);

        assert_eq!(lj.tail_virial(1.0), 17066.666666666668);
        assert_eq!(lj.tail_virial(2.0), -17.06666666666667);
        assert_eq!(lj.tail_virial(14.42), -0.1366035877536718);

        assert_approx_eq!(lj.force(f64::powf(2.0, 1.0/6.0) * 2.0), 0.0);
        assert_approx_eq!(lj.force(2.5), -0.95773475733504);

        let e0 = lj.energy(4.0);
        let e1 = lj.energy(4.0 + EPS);
        assert_approx_eq!((e0 - e1)/EPS, lj.force(4.0), 1e-6);
    }

    #[test]
    fn harmonic() {
        let harm = Harmonic{k: 50.0, x0: 2.0};
        assert_eq!(harm.energy(2.0), 0.0);
        assert_eq!(harm.energy(2.5), 6.25);

        assert_eq!(harm.force(2.0), 0.0);
        assert_eq!(harm.force(2.5), -25.0);

        assert_eq!(harm.tail_energy(1.0), 0.0);
        assert_eq!(harm.tail_virial(1.0), 0.0);

        let e0 = harm.energy(2.1);
        let e1 = harm.energy(2.1 + EPS);
        assert_approx_eq!((e0 - e1)/EPS, harm.force(2.1), 1e-6);
    }

    #[test]
    fn cosine_harmonic() {
        let harm = CosineHarmonic::new(50.0, 2.0);
        assert_eq!(harm.energy(2.0), 0.0);
        let dcos = f64::cos(2.5) - f64::cos(2.0);
        assert_eq!(harm.energy(2.5), 0.5 * 50.0 * dcos * dcos);

        assert_eq!(harm.force(2.0), 0.0);
        let dcos = f64::cos(2.5) - f64::cos(2.0);
        assert_eq!(harm.force(2.5), 50.0 * dcos * f64::sin(2.5));

        let e0 = harm.energy(2.3);
        let e1 = harm.energy(2.3 + EPS);
        assert_approx_eq!((e0 - e1)/EPS, harm.force(2.3), 1e-6);
    }

    #[test]
    fn torsion() {
        let torsion = Torsion{k: 5.0, n: 3, delta: 3.0};
        assert_eq!(torsion.energy(1.0), 10.0);
        let energy = 5.0 * (1.0 + f64::cos(3.0 * 1.1 - 3.0));
        assert_eq!(torsion.energy(1.1), energy);

        assert_eq!(torsion.force(1.0), 0.0);

        let e0 = torsion.energy(4.0);
        let e1 = torsion.energy(4.0 + EPS);
        assert_approx_eq!((e0 - e1)/EPS, torsion.force(4.0), 1e-6);
    }

    #[test]
    fn born() {
        let born = BornMayerHuggins{a: 2.0, c: 1.0, d: 0.5, sigma: 2.0, rho: 2.0};

        // Comparing to externally computed values
        assert_eq!(born.energy(2.0), 1.986328125);
        assert_eq!(born.force(2.0), 0.9609375);

        assert_eq!(born.tail_energy(10.0), 4.981521444402363);
        assert_eq!(born.tail_virial(10.0), 69.13986044386026);

        let e0 = born.energy(4.0);
        let e1 = born.energy(4.0 + EPS);
        assert_approx_eq!((e0 - e1)/EPS, born.force(4.0), 1e-6);
    }
}
