// Cymbalum, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

extern crate rand;
use self::rand::distributions::{Sample, Range};
use self::rand::Rng;

use std::usize;

use super::MCMove;
use super::select_molecule;

use types::Vector3D;
use system::{System, EnergyCache};

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

    fn prepare(&mut self, system: &mut System, rng: &mut Box<Rng>) -> bool {
        if let Some(id) = select_molecule(system, self.moltype, rng) {
            self.molid = id;
        } else {
            warn!("Can not translate molecule: no molecule of this type in the system.");
            return false;
        }

        self.delta.x = self.delta_range.sample(rng);
        self.delta.y = self.delta_range.sample(rng);
        self.delta.z = self.delta_range.sample(rng);

        self.newpos.clear();
        for i in system.molecule(self.molid) {
            self.newpos.push(system[i].position + self.delta);
        }
        return true;
    }

    fn cost(&self, system: &System, beta: f64, cache: &mut EnergyCache) -> f64 {
        let idxes = system.molecule(self.molid).iter().collect::<Vec<_>>();
        let cost = cache.move_particles_cost(system, idxes, &self.newpos);
        return cost/beta;
    }

    fn apply(&mut self, system: &mut System) {
        for (i, pi) in system.molecule(self.molid).iter().enumerate() {
            system[pi].position = self.newpos[i];
        }
    }

    fn restore(&mut self, _: &mut System) {
        // Nothing to do
    }
}
