// Cymbalum, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

use rand::distributions::{Sample, Range};
use rand::Rng;

use std::f64;
use std::mem;

use super::MCMove;

use types::{Matrix3, One};
use sys::{System, Configuration, EnergyCache};

/// Monte Carlo move that changes the size of the simulation cell
pub struct Resize {
    /// Delta for translation of the box length
    delta: f64,
    /// Sampling range for volume scaling
    range: Range<f64>,
    /// Configuration before applying changes to the simulation cell
    previous: Configuration,
    /// target pressure
    pressure: f64,
    /// largest cutoff diameter of potentials in `Interactions`
    maximum_cutoff: Option<f64>,
}

impl Resize {
    /// Create a new `Resize` move, with target pressure `pressure` and maximum
    /// displacement of `delta`.
    pub fn new(pressure: f64, delta: f64) -> Resize {
        assert!(delta > 0.0, "delta must be positive in Resize move");
        Resize {
            delta: delta,
            range: Range::new(-delta, delta),
            previous: Configuration::new(),
            pressure: pressure,
            maximum_cutoff: None,
        }
    }
}

impl MCMove for Resize {
    fn describe(&self) -> &str {
        "resizing of the cell"
    }

    fn setup(&mut self, system: &System) {
        // check if the cell is infinite
        if system.cell.is_infinite() {
            fatal_error!("Cannot use `Resize` move with infinite simulation cell.")
        }

        // Get the largest cutoff of all intermolecular interactions in the
        // system.
        self.maximum_cutoff = system.maximum_cutoff()
    }

    fn prepare(&mut self, system: &mut System, rng: &mut Box<Rng>) -> bool {
        let delta = self.range.sample(rng);

        // Store the previous configuration
        self.previous = (**system).clone();

        let volume = system.volume();
        let scaling_factor = f64::cbrt((volume + delta) / volume);
        // Change the simulation cell
        system.cell.scale_mut(Matrix3::one() * scaling_factor);
        // Check the radius of the smallest inscribed sphere and compare to the
        // cut off distance.
        // Abort simulation when box gets smaller than twice the cutoff radius.
        if let Some(maximum_cutoff) = self.maximum_cutoff {
            if system.cell.lengths().iter().any(|&d| 0.5 * d <= maximum_cutoff) {
                fatal_error!("Tried to decrease the cell size but new size \
                              conflicts with the cut off radius. Increase the number of \
                              particles to get rid of this problem.");
            }
        };

        for (mi, molecule) in self.previous.molecules().iter().enumerate() {
            // We don't want to change the intramolecular distances so we
            // compute the translation vector of the center-of-mass (com) of a
            // molecule and apply it to all its particles. Note that to do
            // this, the com of a molecule *always* has to reside inside the
            // simulation cell.
            let old_com = self.previous.molecule_com(mi);
            let frac_com = self.previous.cell.fractional(&old_com);
            let delta_com = system.cell.cartesian(&frac_com) - old_com;
            for position in &mut system.particles_mut().position[molecule.iter()] {
                *position += delta_com;
            }
        }
        true
    }

    fn cost(&self, system: &System, beta: f64, cache: &mut EnergyCache) -> f64 {
        let delta_energy = cache.move_all_rigid_molecules_cost(system);
        let new_volume = system.volume();
        let old_volume = self.previous.cell.volume();
        let delta_volume = new_volume - old_volume;
        // Build and return the cost function.
        beta * (delta_energy + self.pressure * delta_volume) -
        (system.molecules().len() as f64) * f64::ln(new_volume / old_volume)
    }

    fn apply(&mut self, _: &mut System) {
        // Nothing to do.
    }

    fn restore(&mut self, system: &mut System) {
        // Exchange configurations
        mem::swap(&mut **system, &mut self.previous)
    }

    fn update_amplitude(&mut self, scaling_factor: Option<f64>) {
        if let Some(s) = scaling_factor {
            self.delta *= s;
            self.range = Range::new(-self.delta, self.delta);
        }
    }
}
