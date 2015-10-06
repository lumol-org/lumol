/* Cymbalum, Molecular Simulation in Rust - Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
 */

//! Pair potentials computation

use super::PotentialFunction;

/// A `Computation` is a way to compute a potential.
pub trait Computation {
    /// Compute the compute_energy value at `r`
    fn compute_energy(&self, r: f64) -> f64;

    /// Compute the compute_force value at `r`
    fn compute_force(&self, r: f64) -> f64;
}

/******************************************************************************/
/// Direct computation of a potential, delegating to the underlying
/// `PotentialFunction` all the computations.
pub struct DirectComputation{
    /// Potential to compute
    potential: Box<PotentialFunction>
}

impl DirectComputation {
    /// Create a new DirectComputation for the given `potential`.
    pub fn new(potential: Box<PotentialFunction>) -> DirectComputation {
        DirectComputation{potential: potential}
    }
}

impl Computation for DirectComputation {
    #[inline]
    fn compute_energy(&self, r: f64) -> f64 {
        return self.potential.energy(r);
    }

    #[inline]
    fn compute_force(&self, r: f64) -> f64 {
        return self.potential.force(r);
    }
}

/******************************************************************************/
/// Direct computation of the potential with a cutoff applied. The compute_energy is
/// shifted to ensure `E(rc) = 0`, where `rc` is the cutoff distance.
pub struct CutoffComputation {
    /// Potential to compute
    potential: Box<PotentialFunction>,
    /// Cutoff distance
    cutoff: f64,
    /// Energy at cutoff
    delta: f64,
}

impl CutoffComputation {
    /// Create a new `CutoffComputation` for `potential`, with cutoff distance
    /// of `cutoff`.
    pub fn new(potential: Box<PotentialFunction>, cutoff: f64) -> CutoffComputation {
        let delta = potential.energy(cutoff);
        CutoffComputation{potential: potential, cutoff: cutoff, delta: delta}
    }
}

impl Computation for CutoffComputation {
    #[inline]
    fn compute_energy(&self, r: f64) -> f64 {
        if r > self.cutoff {
            return 0.0;
        } else {
            return self.potential.energy(r) - self.delta;
        }
    }

    #[inline]
    fn compute_force(&self, r: f64) -> f64 {
        if r > self.cutoff {
            return 0.0;
        } else {
            return self.potential.force(r);
        }
    }
}

/******************************************************************************/
/// Computation of a potential using tabulated values. This can be faster than
/// direct computation for smooth potentials, but either uses more memory or is
/// less precise than direct computation. Values are tabulated in the `[0, max)`
/// range, and a cutoff is applied after `max`. Energy is shifted before the
/// cutoff to ensure that `E(max) = 0`
pub struct TableComputation {
    // TODO: use genericity over static values here if it ever comes out
    // see https://internals.rust-lang.org/t/pre-rfc-genericity-over-static-values/1538/19

    /// Number of tabulated values
    N: usize,
    /// Step for tabulated value. compute_energy[i]/compute_force[i] contains compute_energy/compute_force at
    /// r = i * delta
    delta: f64,
    /// Tabulated potential
    compute_energy: Vec<f64>,
    /// Tabulated compute_force
    compute_force: Vec<f64>,
}


impl TableComputation {
    /// Create a new `TableComputation` for `potential`, with `N` points and a
    /// maximum value of `max`.
    pub fn new(potential: Box<PotentialFunction>, N: usize, max:f64) -> TableComputation {
        let delta = max/(N as f64);
        let dE = potential.energy(max);
        let mut compute_energy = Vec::with_capacity(N);
        let mut compute_force = Vec::with_capacity(N);
        for i in 0..N {
            let pos = i as f64 * delta;
            compute_energy.push(potential.energy(pos) - dE);
            compute_force.push(potential.force(pos));
        }
        TableComputation{N: N, delta: delta, compute_energy: compute_energy, compute_force: compute_force}
    }
}

impl Computation for TableComputation {
    fn compute_energy(&self, r: f64) -> f64 {
        let bin = f64::floor(r / self.delta) as usize;
        if bin < self.N - 1 {
            let dx = r - (bin as f64)*self.delta;
            let slope = (self.compute_energy[bin + 1] - self.compute_energy[bin])/self.delta;
            return self.compute_energy[bin] + dx*slope;
        } else {
            return 0.0;
        }
    }

    fn compute_force(&self, r: f64) -> f64 {
        let bin = f64::floor(r / self.delta) as usize;
        if bin < self.N - 1 {
            let dx = r - (bin as f64)*self.delta;
            let slope = (self.compute_force[bin + 1] - self.compute_force[bin])/self.delta;
            return self.compute_force[bin] + dx*slope;
        } else {
            return 0.0;
        }
    }
}

/******************************************************************************/

#[cfg(test)]
mod test {
    use super::*;
    use ::potentials::Harmonic;


    #[test]
    fn direct() {
        let direct = DirectComputation::new(Box::new(Harmonic{k: 50.0, x0: 2.0}));

        assert_eq!(direct.compute_force(2.5), -25.0);
        assert_eq!(direct.compute_energy(2.5), 6.25);
    }

    #[test]
    fn cutoff() {
        let cutoff = CutoffComputation::new(Box::new(Harmonic{k: 50.0, x0: 2.0}), 4.0);

        assert_eq!(cutoff.compute_force(2.5), -25.0);
        assert_eq!(cutoff.compute_energy(2.5), -93.75);

        assert_eq!(cutoff.compute_force(4.1), 0.0);
        assert_eq!(cutoff.compute_energy(4.1), 0.0);
    }

    #[test]
    fn table() {
        let table = TableComputation::new(Box::new(Harmonic{k: 50.0, x0: 2.0}), 1000, 4.0);

        assert_eq!(table.compute_force(2.5), -25.0);
        assert_eq!(table.compute_energy(2.5), -93.75);

        assert_eq!(table.compute_force(4.1), 0.0);
        assert_eq!(table.compute_energy(4.1), 0.0);
    }
}
