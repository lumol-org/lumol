/* Cymbalum, Molecular Simulation in Rust - Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
 */
extern crate rand;
use self::rand::distributions::{Sample, Range};
use self::rand::Rng;

use std::usize;

use super::MCMove;
use super::select_molecule;

use types::Vector3D;
use universe::Universe;

/// Monte-Carlo move for translating a molecule
pub struct Translate {
    /// Type of molecule to translate. `None` means all molecules.
    moltype: Option<u64>,
    /// Index of the molecule to translate
    molid: usize,
    /// Translation vector
    delta: Vector3D,
    /// Translation range for random number generation
    delta_range: Range<f64>,
    /// Potential energy before the move
    e_before: f64,
}

impl Translate {
    /// Create a new `Translate` move, with maximum displacement of `dr`,
    /// translating all the molecules in the system.
    pub fn new(dr: f64) -> Translate {
        Translate::create(dr, None)
    }

    /// Create a new `Translate` move, with maximum displacement of `dr`,
    /// translating only molecules with `moltype` type.
    pub fn with_moltype(dr: f64, moltype: u64) -> Translate {
        Translate::create(dr, Some(moltype))
    }

    /// Factorizing the constructors
    fn create(dr: f64, moltype: Option<u64>) -> Translate {
        assert!(dr > 0.0, "dr must be positive in Translate move");
        let dr = dr / f64::sqrt(3.0);
        Translate {
            moltype: moltype,
            molid: usize::MAX,
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

    fn prepare(&mut self, universe: &mut Universe, rng: &mut Box<Rng>) -> bool {
        if let Some(id) = select_molecule(universe, self.moltype, rng) {
            self.molid = id;
        } else {
            warn!("Can not translate molecule: no molecule of this type in the universe.");
            return false;
        }

        self.delta.x = self.delta_range.sample(rng);
        self.delta.y = self.delta_range.sample(rng);
        self.delta.z = self.delta_range.sample(rng);

        self.e_before = universe.potential_energy();
        for i in universe.molecule(self.molid) {
            universe[i].position = universe[i].position + self.delta;
        }
        return true;
    }

    fn cost(&self, universe: &Universe, beta: f64) -> f64 {
        let e_after = universe.potential_energy();
        return (e_after - self.e_before)/beta;
    }

    fn apply(&mut self, _: &mut Universe) {
        // Nothing to do
    }

    fn restore(&mut self, universe: &mut Universe) {
        for i in universe.molecule(self.molid) {
            universe[i].position = universe[i].position - self.delta;
        }
    }
}
