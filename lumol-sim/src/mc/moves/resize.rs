// Cymbalum, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

use rand::RngCore;
use rand_distr::{Uniform, Distribution};

use std::f64;
use std::mem;

use super::{MCDegreeOfFreedom, MCMove};

use lumol_core::{Configuration, EnergyCache, System, Matrix3};

/// Monte Carlo move that changes the size of the simulation cell
pub struct Resize {
    /// Delta for translation of the box length
    delta: f64,
    /// Sampling range for volume scaling
    range: Uniform<f64>,
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
            range: Uniform::new(-delta, delta),
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

    fn degrees_of_freedom(&self) -> MCDegreeOfFreedom {
        MCDegreeOfFreedom::AllMolecules
    }

    fn setup(&mut self, system: &System) {
        // check if the cell is infinite
        assert!(!system.cell.is_infinite(), "cannot use `Resize` move with infinite simulation cell");

        // Get the largest cutoff of all intermolecular interactions in the
        // system.
        self.maximum_cutoff = system.maximum_cutoff();
    }

    #[allow(clippy::manual_assert)]
    fn prepare(&mut self, system: &mut System, rng: &mut dyn RngCore) -> bool {
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
                panic!(
                    "Tried to decrease the cell size in Monte Carlo Resize move \
                     but the new size is smaller than the interactions cut off \
                     radius. You can try to increase the cell size or the number \
                     of particles."
                );
            }
        };

        let cell = system.cell;
        for mut molecule in system.molecules_mut() {
            // We don't want to change the intramolecular distances so we
            // compute the translation vector of the center-of-mass (com) of a
            // molecule and apply it to all its particles. Note that to do
            // this, the com of a molecule *always* has to reside inside the
            // simulation cell.
            let old_com = molecule.as_ref().center_of_mass();
            let frac_com = self.previous.cell.fractional(&old_com);
            let delta_com = cell.cartesian(&frac_com) - old_com;
            for position in molecule.particles_mut().position.iter_mut() {
                *position += delta_com;
            }
        }
        true
    }

    fn cost(&self, system: &System, beta: f64, cache: &mut EnergyCache) -> f64 {
        let delta_energy = cache.move_all_molecules_cost(system);
        let new_volume = system.volume();
        let old_volume = self.previous.cell.volume();
        let delta_volume = new_volume - old_volume;
        // Build and return the cost function.
        beta * (delta_energy + self.pressure * delta_volume)
            - (system.molecules().count() as f64) * f64::ln(new_volume / old_volume)
    }

    fn apply(&mut self, _: &mut System) {
        // Nothing to do.
    }

    fn restore(&mut self, system: &mut System) {
        // Exchange configurations
        mem::swap(&mut **system, &mut self.previous);
    }

    fn update_amplitude(&mut self, scaling_factor: Option<f64>) {
        if let Some(s) = scaling_factor {
            self.delta *= s;
            self.range = Uniform::new(-self.delta, self.delta);
        }
    }
}
