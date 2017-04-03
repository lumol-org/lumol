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
    /// System after applying changes to the simulation cell
    new_system: System,
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
        Resize {
            delta: delta,
            range: Range::new(-delta, delta),
            new_system: System::new(),
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
        // check if the cell is infinite
        if system.cell().is_infinite() {
            fatal_error!("Cannot use `Resize` move with infinite simulation cell.")
        }

        // Get the largest cutoff of all intermolecular interactions in the
        // system.

        // Go through global interactions
        let rc_glob = system.interactions()
                            .globals()
                            .iter()
                            .map(|i| i.cutoff())
                            .filter_map(|rc| rc)
                            .fold(f64::NAN, f64::max);

        // Pair interactions
        let rc_pairs = system.interactions()
                            .all_pairs()
                            .iter()
                            .map(|i| i.cutoff())
                            .fold(f64::NAN, f64::max);

        self.rc_max = f64::max(rc_glob, rc_pairs)
    }

    fn prepare(&mut self, system: &mut System, rng: &mut Box<Rng>) -> bool {
        let delta = self.range.sample(rng);

        // Copy the system: the proposed state will be stored here
        // TODO: we only need to store positions and the cell; all
        // other information stays the same.
        self.new_system = system.clone();

        let volume = system.volume();
        let scaling_factor = f64::cbrt((volume + delta) / volume);
        // Change the simulation cell
        self.new_system.cell_mut().scale_mut(Matrix3::one() * scaling_factor);
        // Check the radius of the smallest inscribed sphere and compare to the
        // cut off distance.
        // Abort simulation when box gets smaller than twice the cutoff radius.
        if self.new_system.cell()
            .lengths()
            .iter()
            .any(|&d| 0.5 * d <= self.rc_max) {
            fatal_error!("Tried to decrease the cell size but new size
                conflicts with the cut off radius. \
                Increase the number of particles to get rid of this problem.")
        }
        // Loop over all molecules in the system.
        // We don't want to change the intramolecular distances
        // so we compute the translation vector of the center-of-mass
        // (com) of a molecule and apply it to all its particles.
        // Note that to do this, the com of a molecule *always* has
        // to reside inside the simulation cell.

        // TODO: Check if system.size == system.molecules().len
        // if that is the case, skip com computation since it is a
        // system without molecules.
        for (mi, molecule) in system
            .molecules()
            .iter()
            .enumerate() {
            let old_com = system.molecule_com(mi);
            let frac_com = system.cell().fractional(&old_com);
            // compute translation vector
            let delta_com =
                self.new_system.cell().cartesian(&frac_com) - old_com;
            // loop over all particles (indices) in the molecule
            for pi in molecule.iter() {
                self.new_system[pi].position += delta_com;
            }
        }
        true
    }

    fn cost(&self, system: &System, beta: f64, cache: &mut EnergyCache) -> f64 {
        let delta_energy = cache.move_all_rigid_molecules_cost(&self.new_system);
        let new_volume = self.new_system.volume();
        let old_volume = system.volume();
        let delta_volume = new_volume - old_volume;
        // Build and return the cost function.
        beta * (delta_energy + self.pressure * delta_volume) -
        (system.molecules().len() as f64) * f64::ln(new_volume / old_volume)
    }

    fn apply(&mut self, system: &mut System) {
        // Exchange systems: This will effectively update
        // new positions and cell.
        mem::swap(system, &mut self.new_system)
    }

    fn restore(&mut self, _: &mut System) {
        // Do nothing.
    }

    fn update_amplitude(&mut self, scaling_factor: Option<f64>) {
        if let Some(s) = scaling_factor {
            self.delta *= s;
            self.range = Range::new(-self.delta, self.delta);
        }
    }
}
