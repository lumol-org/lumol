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
use universe::{Universe, EnergyCache};

/// Monte-Carlo move for translating a molecule
pub struct Translate {
    /// Type of molecule to translate. `None` means all molecules.
    moltype: Option<u64>,
    /// Index of the molecule to translate
    molid: usize,
    /// New positions of the atom in the translated molecule
    newpos: Vec<Vector3D>,
    /// Translation vector
    delta: Vector3D,
    /// Translation range for random number generation
    delta_range: Range<f64>,
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
            newpos: Vec::new(),
            delta: Vector3D::new(0.0, 0.0, 0.0),
            delta_range: Range::new(-dr, dr),
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

        self.newpos.clear();
        for i in universe.molecule(self.molid) {
            self.newpos.push(universe[i].position + self.delta);
        }
        return true;
    }

    fn cost(&self, universe: &Universe, beta: f64, cache: &mut EnergyCache) -> f64 {
        let idxes = universe.molecule(self.molid).iter().collect::<Vec<_>>();
        let cost = cache.move_particles_cost(universe, idxes, &self.newpos);
        return cost/beta;
    }

    fn apply(&mut self, universe: &mut Universe) {
        for (i, pi) in universe.molecule(self.molid).iter().enumerate() {
            universe[pi].position = self.newpos[i];
        }
    }

    fn restore(&mut self, _: &mut Universe) {
        // Nothing to do
    }
}
