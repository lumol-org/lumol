// Cymbalum, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

//! A potential computation is a way of computing a potential given its
//! expression (represented by a `PotentialFunction`). The same potential can be
//! computed either direcly, or using a cutoff, or by a table interpolation, ...
use super::{PotentialFunction, PairPotential};

/// A `Computation` is a way to compute a potential.
pub trait Computation: Sync + Send {
    // / Kind of potential we can apply this computation to
    // type Potential: PotentialFunction + ?Sized;

    /// Compute the energy value at `r`
    fn compute_energy(&self, r: f64) -> f64;

    /// Compute the force value at `r`
    fn compute_force(&self, r: f64) -> f64;
}

impl<P: Computation + Clone + 'static> PotentialFunction for P {
    #[inline] fn energy(&self, r:f64) -> f64 {
        self.compute_energy(r)
    }

    #[inline] fn force(&self, r:f64) -> f64 {
        self.compute_force(r)
    }
}

/******************************************************************************/
/// Direct computation of the potential with a cutoff applied. The computed
/// energy is shifted to ensure `E(rc) = 0`, where `rc` is the cutoff distance.
#[derive(Clone)]
pub struct CutoffComputation {
    /// Potential to compute
    potential: Box<PairPotential>,
    /// Cutoff distance
    cutoff: f64,
    /// Energy at cutoff
    delta: f64,
}

impl CutoffComputation {
    /// Create a new `CutoffComputation` for `potential`, with cutoff distance
    /// of `cutoff`.
    pub fn new(potential: Box<PairPotential>, cutoff: f64) -> CutoffComputation {
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

impl PairPotential for CutoffComputation {}

/******************************************************************************/
/// Computation of a potential using tabulated values. This can be faster than
/// direct computation for smooth potentials, but will either uses more memory
/// or is less precise than direct computation. Values are tabulated in the `[0,
/// max)` range, and a cutoff is applied after `max`. Energy is shifted before
/// the cutoff to ensure that `E(max) = 0`
#[derive(Clone)]
pub struct TableComputation {
    // TODO: use genericity over static values here if it ever comes out
    // see https://internals.rust-lang.org/t/pre-rfc-genericity-over-static-values/1538/19

    /// Number of tabulated values
    size: usize,
    /// Step for tabulated value. compute_energy[i]/compute_force[i] contains compute_energy/compute_force at
    /// r = i * delta
    delta: f64,
    /// Tabulated potential
    compute_energy: Vec<f64>,
    /// Tabulated compute_force
    compute_force: Vec<f64>,
}


impl TableComputation {
    /// Create a new `TableComputation` for `potential`, with `size` points and
    /// a maximum value of `max`.
    pub fn new(potential: Box<PairPotential>, size: usize, max:f64) -> TableComputation {
        let delta = max/(size as f64);
        let energy_shift = potential.energy(max);
        let mut compute_energy = Vec::with_capacity(size);
        let mut compute_force = Vec::with_capacity(size);
        for i in 0..size {
            let pos = i as f64 * delta;
            compute_energy.push(potential.energy(pos) - energy_shift);
            compute_force.push(potential.force(pos));
        }
        TableComputation{size: size, delta: delta, compute_energy: compute_energy, compute_force: compute_force}
    }
}

impl Computation for TableComputation {
    fn compute_energy(&self, r: f64) -> f64 {
        let bin = f64::floor(r / self.delta) as usize;
        if bin < self.size - 1 {
            let dx = r - (bin as f64)*self.delta;
            let slope = (self.compute_energy[bin + 1] - self.compute_energy[bin])/self.delta;
            return self.compute_energy[bin] + dx*slope;
        } else {
            return 0.0;
        }
    }

    fn compute_force(&self, r: f64) -> f64 {
        let bin = f64::floor(r / self.delta) as usize;
        if bin < self.size - 1 {
            let dx = r - (bin as f64)*self.delta;
            let slope = (self.compute_force[bin + 1] - self.compute_force[bin])/self.delta;
            return self.compute_force[bin] + dx*slope;
        } else {
            return 0.0;
        }
    }
}

impl PairPotential for TableComputation {}

/******************************************************************************/

#[cfg(test)]
mod test {
    use super::*;
    use potentials::Harmonic;

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
