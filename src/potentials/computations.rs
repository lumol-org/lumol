/*
 * Cymbalum, Molecular Simulation in Rust
 * Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
*/

//! Pair potentials computation

use super::PotentialFunction;

/// A `PotentialComputation` is a way to compute a potential.
pub trait PotentialComputation {
    /// Compute the energy value at `r`
    fn energy(&self, r: f64) -> f64;

    /// Compute the force value at `r`
    fn force(&self, r: f64) -> f64;
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

impl PotentialComputation for DirectComputation {
    #[inline]
    fn energy(&self, r: f64) -> f64 {
        return self.potential.energy(r);
    }

    #[inline]
    fn force(&self, r: f64) -> f64 {
        return self.potential.force(r);
    }
}

/******************************************************************************/
/// Direct computation of the potential with a cutoff applied. The energy is
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

impl PotentialComputation for CutoffComputation {
    #[inline]
    fn energy(&self, r: f64) -> f64 {
        if r > self.cutoff {
            return 0.0;
        } else {
            return self.potential.energy(r) - self.delta;
        }
    }

    #[inline]
    fn force(&self, r: f64) -> f64 {
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
    /// Step for tabulated value. energy[i]/force[i] contains energy/force at
    /// r = i * delta
    delta: f64,
    /// Tabulated potential
    energy: Vec<f64>,
    /// Tabulated force
    force: Vec<f64>,
}


impl TableComputation {
    /// Create a new `TableComputation` for `potential`, with `N` points and a
    /// maximum value of `max`.
    pub fn new(potential: Box<PotentialFunction>, N: usize, max:f64) -> TableComputation {
        let delta = max/(N as f64);
        let dE = potential.energy(max);
        let mut energy = Vec::with_capacity(N);
        let mut force = Vec::with_capacity(N);
        for i in 0..N {
            let pos = i as f64 * delta;
            energy.push(potential.energy(pos) - dE);
            force.push(potential.force(pos));
        }
        TableComputation{N: N, delta: delta, energy: energy, force: force}
    }
}

impl PotentialComputation for TableComputation {
    fn energy(&self, r: f64) -> f64 {
        let bin = f64::floor(r / self.delta) as usize;
        if bin < self.N - 1 {
            let dx = r - (bin as f64)*self.delta;
            let slope = (self.energy[bin + 1] - self.energy[bin])/self.delta;
            return self.energy[bin] + dx*slope;
        } else {
            return 0.0;
        }
    }

    fn force(&self, r: f64) -> f64 {
        let bin = f64::floor(r / self.delta) as usize;
        if bin < self.N - 1 {
            let dx = r - (bin as f64)*self.delta;
            let slope = (self.force[bin + 1] - self.force[bin])/self.delta;
            return self.force[bin] + dx*slope;
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
        let direct = DirectComputation::new(Box::new(Harmonic{k: 50.0, r0: 2.0}));

        assert_eq!(direct.force(2.5), -25.0);
        assert_eq!(direct.energy(2.5), 6.25);
    }

    #[test]
    fn cutoff() {
        let cutoff = CutoffComputation::new(Box::new(Harmonic{k: 50.0, r0: 2.0}), 4.0);

        assert_eq!(cutoff.force(2.5), -25.0);
        assert_eq!(cutoff.energy(2.5), -93.75);

        assert_eq!(cutoff.force(4.1), 0.0);
        assert_eq!(cutoff.energy(4.1), 0.0);
    }

    #[test]
    fn table() {
        let table = TableComputation::new(Box::new(Harmonic{k: 50.0, r0: 2.0}), 1000, 4.0);

        assert_eq!(table.force(2.5), -25.0);
        assert_eq!(table.energy(2.5), -93.75);

        assert_eq!(table.force(4.1), 0.0);
        assert_eq!(table.energy(4.1), 0.0);
    }
}
