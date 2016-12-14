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
        // TODO: set rc_max from the largest cut off in PairPotentials.
        // For the time being, we set rc_max to zero.
        // That way, the simulation cell can get smaller than 2*rc_max
        // which can lead to wrong energies.
        Resize {
            delta: delta,
            range: Range::new(-delta, delta),
            old_system: System::new(),
            pressure: pressure,
            rc_max: 0.0,
        }
    }
}

impl MCMove for Resize {
    fn describe(&self) -> &str {
        "resizing of the cell"
    }

    fn setup(&mut self, system: &System) {
        // get the largest cutoff of all intermolecular interactions in the
        // system.
        // TODO: include elecetrostatic interactions
        self.rc_max = 0f64;
        // iterate over all (intermolecular) pair interactions to extract rc
        for rc in system.interactions()
            .all_pairs()
            .iter()
            .map(|i| i.get_cutoff()) {
            if rc > self.rc_max {
                self.rc_max = rc
            }
        }
    }

    fn prepare(&mut self, system: &mut System, rng: &mut Box<Rng>) -> bool {
        let delta = self.range.sample(rng);
        // clone the current (old) system
        self.old_system = system.clone();
        let volume = system.volume();
        let scaling_factor = f64::cbrt((volume + delta) / volume);
        // change the simulation cell
        system.cell_mut().scale_mut(Matrix3::one() * scaling_factor);
        let new_cell = system.cell().clone();
        // check the radius of the smallest inscribed sphere and compare to the
        // cut off distance.
        // abort simulation when box gets smaller than twice the cutoff radius
        if new_cell.lengths().iter().any(|&d| 0.5 * d <= self.rc_max) {
            fatal_error!("Tried to decrease the cell size but new size 
                conflicts with the cut off radius. \
                Increase the number of particles to get rid of this problem.")
        }
        // loop over all molecules in the system
        // we don't want to change the intramolecular distances
        // so we compute the translation vector of the center of mass
        // of a molecule and apply it to all particles
        
        // TODO: check if system.size == system.molecules().len
        // if that is the case, skip com computation since it is a
        // system without molecules
        for (mi, molecule) in self.old_system
            .molecules()
            .iter()
            .enumerate() {
            let com_old = self.old_system.molecule_com(mi);
            let com_frac = self.old_system
                .cell()
                .fractional(&com_old);
            // compute translation vector
            let delta_com = new_cell.cartesian(&com_frac) - com_old;
            // loop over all particles (indices) in the molecule
            for pi in molecule.iter() {
                system[pi].position += delta_com;
            }
        }
        true
    }

    fn cost(&self, system: &System, beta: f64, cache: &mut EnergyCache) -> f64 {
        let delta_energy = cache.resize_cell_cost(system);
        let new_volume = system.volume();
        let old_volume = self.old_system.volume();
        let delta_volume = new_volume - old_volume;

        // return the cost function
        beta * (delta_energy + self.pressure * delta_volume) -
        (system.size() as f64) * f64::ln(new_volume / old_volume)
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
