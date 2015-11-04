/* Cymbalum, Molecular Simulation in Rust - Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
 */
//! A Monte-Carlo simulation is based on a set of applicables moves.
//!
//! For example, NVT Monte-Carlo will use the `Translate` move only for particles
//! and add the `Rotate` moves for molecules. NPT Monte-Carlo will add the
//! `VolumeResize` move.
//!
//! In all this module, beta refers to the Boltzmann factor 1/(kB T)
extern crate rand;
use self::rand::distributions::{Sample, Range};
use self::rand::Rng;

use std::usize;

use types::Vector3D;
use universe::Universe;

/// The `MCMove` trait correspond to the set of methods used in Monte-Carlo
/// simulations.
pub trait MCMove {
    /// Give a short description of this move
    fn describe(&self) -> &str;

    /// Prepare the move, by selecting the particles to move, and the parameters
    /// of the move. The `rng` random number generator should be used to
    /// generate the parameters of the move.
    fn prepare(&mut self, universe: &mut Universe, rng: &mut Box<Rng>);

    /// Get the cost of performing this move, as the exponential factor. For
    /// simple NVT simulations, this cost is the energetic difference over
    /// `beta`.
    ///
    /// The cost must be dimmensionless, and will be placed in an exponential.
    fn cost(&self, universe: &Universe, beta: f64) -> f64;

    /// Effectivelly apply the move, if it has not already been done in
    /// `prepare`.
    fn apply(&mut self, universe: &mut Universe);

    /// Restore the universe to it's initial state, if it has been changed in
    /// `prepare`.
    fn restore(&mut self, universe: &mut Universe);
}


/******************************************************************************/

/// Monte-Carlo move for translating a molecule
pub struct Translate {
    /// Index of the molecule to translate
    molecule: usize,
    /// Translation vector
    delta: Vector3D,
    /// Translation range
    delta_range: Range<f64>,
    /// Potential energy before the move
    e_before: f64,
}

impl Translate {
    /// Create a new `Translate` move, with maximum displacement of `dr`.
    pub fn new(dr: f64) -> Translate {
        assert!(dr > 0.0, "dr must be positive in Translate move");
        let dr = dr / f64::sqrt(3.0);
        Translate {
            molecule: usize::MAX,
            delta: Vector3D::new(0.0, 0.0, 0.0),
            delta_range: Range::new(-dr, dr),
            e_before: 0.0,
        }
    }
}

impl Default for Translate {
    fn default() -> Translate {
        Translate::new(1.0)
    }
}

impl MCMove for Translate {
    fn describe(&self) -> &str {
        "molecular translation"
    }

    fn prepare(&mut self, universe: &mut Universe, rng: &mut Box<Rng>) {
        self.delta.x = self.delta_range.sample(rng);
        self.delta.y = self.delta_range.sample(rng);
        self.delta.z = self.delta_range.sample(rng);

        let nmols = universe.molecules().len();
        self.molecule = rng.gen_range(0, nmols);

        self.e_before = universe.potential_energy();
        for i in universe.molecules()[self.molecule].clone() {
            universe[i].position = universe[i].position + self.delta;
        }
    }

    fn cost(&self, universe: &Universe, beta: f64) -> f64 {
        let e_after = universe.potential_energy();
        return (e_after - self.e_before)/beta;
    }

    fn apply(&mut self, _: &mut Universe) {
        // Nothing to do
    }

    fn restore(&mut self, universe: &mut Universe) {
        for i in universe.molecules()[self.molecule].clone() {
            universe[i].position = universe[i].position - self.delta;
        }
    }
}
