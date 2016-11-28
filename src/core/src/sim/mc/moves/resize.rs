// Cymbalum, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

use rand::distributions::{Sample, Range};
use rand::Rng;

use std::f64;
use std::mem;

use super::MCMove;

use types::{Matrix3, One};
use sys::{System, EnergyCache};

/// Monte-Carlo move that changes the size of the simulation cell
pub struct Resize {
    /// Delta for translation of the box length
    delta: f64,
    /// Sampling range for volume scaling
    range: Range<f64>,
    /// System before applying changes to the simulation box
    old_system: System,
    /// target pressure
    pressure: f64,
    /// largest cut off diameter of `PairPotentials`
    rc_max: f64,
}

impl Resize {
    /// Create a new `Resize` move, with target pressure `pressure` and maximum
    /// displacement of `delta`.
    pub fn new(pressure: f64, delta: f64) -> Resize {
        assert!(delta > 0.0, "delta must be positive in Resize move");
        // TODO: set rc_max from the largest cut off in PairPotentials
        Resize {
            delta: 0.,
            range: Range::new(-delta, delta),
            old_system: System::new(),
            pressure: pressure,
            rc_max: 100.0, // arbitrarily large value
        }
    }
}

impl MCMove for Resize {
    fn describe(&self) -> &str {
        "resizing of the cell"
    }

    fn prepare(&mut self, system: &mut System, rng: &mut Box<Rng>) -> bool {
        self.delta = self.range.sample(rng);
        // clone the current (old) system
        self.old_system = system.clone();
        let volume = system.volume();
        let scaling_factor = f64::cbrt((volume + self.delta) / volume);
        // change the simulation cell
        system.cell_mut().scale_mut(Matrix3::one() * scaling_factor);
        let new_cell = system.cell().clone();
        // check the radius of the smallest inscribed sphere and compare to the 
        // cut off distance.
        // abort simulation when box gets too small
        if new_cell.lengths_as_vec().iter().any(|&d| 0.5 * d < self.rc_max) {
            fatal_error!(
                "Tried to decrease the cell size but \
                new size conflicts with the cut off radius. \
                You could try to increase the number of particles to get rid \
                of this problem."
            )
        }
        for particle in system.iter_mut() {
            // when changing the size (or shape) of the cell, fractional
            // coordinates do not change.
            // compute fractional coordinates using the old cell
            let fractional = self.old_system.cell().fractional(&particle.position);
            // compute coordinates after resizing the cell
            particle.position = new_cell.cartesian(&fractional);
        }
        true
    }

    fn cost(&self, system: &System, beta: f64, cache: &mut EnergyCache) -> f64 {
        cache.unused();
        // get system energy before change: this is stored in `self.old_system`
        let old_energy = self.old_system.potential_energy();
        // get
        let new_energy = system.potential_energy();
        let delta_energy = new_energy - old_energy;

        let new_volume = system.volume();
        let old_volume = self.old_system.volume();
        let delta_volume = new_volume - old_volume;

        // return the cost function
        beta * (delta_energy + self.pressure * delta_volume)
            - (system.size() as f64) * f64::ln(new_volume / old_volume)
    }

    fn apply(&mut self, _: &mut System) {
        // Nothing to do
    }

    fn restore(&mut self, system: &mut System) {
        mem::swap(system, &mut self.old_system);
    }

    fn update_amplitude(&mut self, scaling_factor: Option<f64>) {
        if let Some(s) = scaling_factor {
            self.delta *= s;
            self.range = Range::new(-self.delta, self.delta);
        }
    }
}
