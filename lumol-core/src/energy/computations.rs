// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

use log_once::warn_once;

use crate::{PairPotential, Potential};

/// Alternative energy and forces computation.
///
/// The `Computation` trait represent an alternative way to compute a given
/// potential. For example using interpolation on a table or on a grid, from a
/// Fourier decomposition, *etc.*
///
/// # Examples
///
/// ```
/// # use lumol_core::energy::Potential;
/// use lumol_core::energy::Computation;
/// use lumol_core::energy::Harmonic;
///
/// /// This is just a thin wrapper logging every time the `energy/force`
/// /// methods are called.
/// #[derive(Clone)]
/// struct LoggingComputation<T: Potential>(T);
///
/// impl<T: Potential> Computation for LoggingComputation<T> {
///     fn compute_energy(&self, r: f64) -> f64 {
///         println!("Called energy");
///         self.0.energy(r)
///     }
///
///     fn compute_force(&self, r: f64) -> f64 {
///         println!("Called force");
///         self.0.force(r)
///     }
/// }
///
/// let potential = Harmonic{x0: 0.5, k: 4.2};
/// let computation = LoggingComputation(potential.clone());
///
/// assert_eq!(computation.energy(1.0), potential.energy(1.0));
/// assert_eq!(computation.force(2.0), potential.force(2.0));
/// ```
pub trait Computation: Sync + Send {
    /// Compute the energy value at `r`
    fn compute_energy(&self, r: f64) -> f64;

    /// Compute the force value at `r`
    fn compute_force(&self, r: f64) -> f64;
}

impl<P: Computation + Clone + 'static> Potential for P {
    #[inline]
    fn energy(&self, r: f64) -> f64 {
        self.compute_energy(r)
    }

    #[inline]
    fn force(&self, r: f64) -> f64 {
        self.compute_force(r)
    }
}

/// Computation of a potential using tabulated values.
///
/// This can be faster than direct computation for smooth potentials, but will
/// uses more memory and be less precise than direct computation. Values are
/// tabulated in the `[0, max)` range, and a cutoff is applied after `max`.
#[derive(Clone)]
pub struct TableComputation {
    /// Step for tabulated value. `energy_table[i]`/`force_table[i]` contains
    /// energy/force at `r = i * delta`
    delta: f64,
    /// Cutoff distance
    cutoff: f64,
    /// Tabulated potential
    energy_table: Vec<f64>,
    /// Tabulated compute_force
    force_table: Vec<f64>,
    /// Initial potential, kept around for tail corrections
    potential: Box<dyn PairPotential>,
}


impl TableComputation {
    /// Create a new `TableComputation` for `potential`, with `size` points and
    /// a maximum value of `max`.
    ///
    /// # Examples
    ///
    /// ```
    /// # use lumol_core::energy::Potential;
    /// use lumol_core::energy::TableComputation;
    /// use lumol_core::energy::Harmonic;
    ///
    /// let potential = Box::new(Harmonic{x0: 0.5, k: 4.2});
    /// let table = TableComputation::new(potential, 1000, 2.0);
    ///
    /// assert_eq!(table.energy(1.0), 0.525);
    /// assert_eq!(table.energy(3.0), 0.0);
    /// ```
    pub fn new(potential: Box<dyn PairPotential>, size: usize, max: f64) -> TableComputation {
        let delta = max / (size as f64);
        let mut energy_table = Vec::with_capacity(size);
        let mut force_table = Vec::with_capacity(size);
        for i in 0..size {
            let r = i as f64 * delta;
            energy_table.push(potential.energy(r));
            force_table.push(potential.force(r));
        }

        TableComputation {
            delta: delta,
            cutoff: max,
            energy_table: energy_table,
            force_table: force_table,
            potential: potential,
        }
    }
}

impl Computation for TableComputation {
    fn compute_energy(&self, r: f64) -> f64 {
        debug_assert_eq!(self.energy_table.len(), self.force_table.len());
        let bin = f64::floor(r / self.delta) as usize;
        if bin < self.energy_table.len() - 1 {
            let dx = r - (bin as f64) * self.delta;
            let slope = (self.energy_table[bin + 1] - self.energy_table[bin]) / self.delta;
            return self.energy_table[bin] + dx * slope;
        }

        return 0.0;
    }

    fn compute_force(&self, r: f64) -> f64 {
        debug_assert_eq!(self.energy_table.len(), self.force_table.len());
        let bin = f64::floor(r / self.delta) as usize;
        if bin < self.force_table.len() - 1 {
            let dx = r - (bin as f64) * self.delta;
            let slope = (self.force_table[bin + 1] - self.force_table[bin]) / self.delta;
            return self.force_table[bin] + dx * slope;
        }

        return 0.0;
    }
}

impl PairPotential for TableComputation {
    fn tail_energy(&self, cutoff: f64) -> f64 {
        if cutoff > self.cutoff {
            warn_once!(
                "Cutoff in table computation ({}) is smaller than the \
                 pair interaction cutoff ({}) when computing tail correction. \
                 This may lead to wrong values for energy.",
                cutoff,
                self.cutoff
            );
        }
        return self.potential.tail_energy(cutoff);
    }

    fn tail_virial(&self, cutoff: f64) -> f64 {
        if cutoff > self.cutoff {
            warn_once!(
                "Cutoff in table computation ({}) is smaller than the \
                 pair interaction cutoff ({}) when computing tail correction. \
                 This may lead to wrong values for pressure.",
                cutoff,
                self.cutoff
            );
        }
        return self.potential.tail_virial(cutoff);
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::{Harmonic, LennardJones};
    use crate::PairPotential;

    #[test]
    fn table() {
        let table = TableComputation::new(Box::new(Harmonic { k: 50.0, x0: 2.0 }), 1000, 4.0);

        assert_eq!(table.compute_energy(2.5), 6.25);
        assert_eq!(table.compute_force(2.5), -25.0);

        // Check that the table is defined up to the cutoff value
        let delta = 4.0 / 1000.0;
        assert_eq!(table.compute_energy(4.0 - 2.0 * delta), 99.2016);
        assert_eq!(table.compute_force(4.0 - 2.0 * delta), -99.6);

        assert_eq!(table.compute_energy(4.0 - delta), 0.0);
        assert_eq!(table.compute_force(4.0 - delta), 0.0);

        assert_eq!(table.compute_energy(4.0), 0.0);
        assert_eq!(table.compute_force(4.0), 0.0);

        assert_eq!(table.compute_energy(4.1), 0.0);
        assert_eq!(table.compute_force(4.1), 0.0);


        let lj = LennardJones {
            epsilon: 50.0,
            sigma: 2.0,
        };
        let table = TableComputation::new(Box::new(lj), 1000, 4.0);
        assert_eq!(table.tail_energy(5.0), lj.tail_energy(5.0));
        assert_eq!(table.tail_virial(5.0), lj.tail_virial(5.0));
    }
}
