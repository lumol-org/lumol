// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

use energy::{Potential, PairPotential};

/// Alternative energy and forces computation.
///
/// The `Computation` trait represent an alternative way to compute a given
/// potential. For example using interpolation on a table or on a grid, from a
/// Fourrier deecompostion, *etc.*
pub trait Computation: Sync + Send {
    /// Compute the energy value at `r`
    fn compute_energy(&self, r: f64) -> f64;

    /// Compute the force value at `r`
    fn compute_force(&self, r: f64) -> f64;
}

impl<P: Computation + Clone + 'static> Potential for P {
    #[inline] fn energy(&self, r:f64) -> f64 {
        self.compute_energy(r)
    }

    #[inline] fn force(&self, r:f64) -> f64 {
        self.compute_force(r)
    }
}

/// Computation of a potential using tabulated values.
///
/// This can be faster than direct computation for smooth potentials, but will
/// uses more memory and be less precise than direct computation. Values are
/// tabulated in the `[0, max)` range, and a cutoff is applied after `max`.
/// Energy is shifted to ensure `E(max) = 0`
#[derive(Clone)]
pub struct TableComputation {
    /// Step for tabulated value. `energy_table[i]`/`force_table[i]` contains
    /// energy/force at `r = i * delta`
    delta: f64,
    /// Tabulated potential
    energy_table: Vec<f64>,
    /// Tabulated compute_force
    force_table: Vec<f64>,
}


impl TableComputation {
    /// Create a new `TableComputation` for `potential`, with `size` points and
    /// a maximum value of `max`.
    pub fn new(potential: &PairPotential, size: usize, max:f64) -> TableComputation {
        let delta = max / (size as f64);
        let energy_shift = potential.energy(max);
        let mut energy_table = Vec::with_capacity(size);
        let mut force_table = Vec::with_capacity(size);
        for i in 0..size {
            let pos = i as f64 * delta;
            energy_table.push(potential.energy(pos) - energy_shift);
            force_table.push(potential.force(pos));
        }
        TableComputation {
            delta: delta,
            energy_table: energy_table,
            force_table: force_table
        }
    }
}

impl Computation for TableComputation {
    fn compute_energy(&self, r: f64) -> f64 {
        debug_assert!(self.energy_table.len() == self.force_table.len());
        let bin = f64::floor(r / self.delta) as usize;
        if bin < self.energy_table.len() - 1 {
            let dx = r - (bin as f64) * self.delta;
            let slope = (self.energy_table[bin + 1] - self.energy_table[bin]) / self.delta;
            return self.energy_table[bin] + dx * slope;
        } else {
            return 0.0;
        }
    }

    fn compute_force(&self, r: f64) -> f64 {
        debug_assert!(self.energy_table.len() == self.force_table.len());
        let bin = f64::floor(r / self.delta) as usize;
        if bin < self.force_table.len() - 1 {
            let dx = r - (bin as f64) * self.delta;
            let slope = (self.force_table[bin + 1] - self.force_table[bin]) / self.delta;
            return self.force_table[bin] + dx * slope;
        } else {
            return 0.0;
        }
    }
}

impl PairPotential for TableComputation {}

#[cfg(test)]
mod test {
    use super::*;
    use energy::Harmonic;

    #[test]
    fn table() {
        let table = TableComputation::new(&Harmonic{k: 50.0, x0: 2.0}, 1000, 4.0);

        assert_eq!(table.compute_energy(2.5), -93.75);
        assert_eq!(table.compute_force(2.5), -25.0);

        // Check that the table is defined up to the cutoff value
        let delta = 4.0 / 1000.0;
        assert_eq!(table.compute_energy(4.0 - 2.0 * delta), -0.7984000000000009);
        assert_eq!(table.compute_force(4.0 - 2.0 * delta), -99.6);

        assert_eq!(table.compute_energy(4.0 - delta), 0.0);
        assert_eq!(table.compute_force(4.0 - delta), 0.0);

        assert_eq!(table.compute_energy(4.0), 0.0);
        assert_eq!(table.compute_force(4.0), 0.0);

        assert_eq!(table.compute_energy(4.1), 0.0);
        assert_eq!(table.compute_force(4.1), 0.0);
    }
}
