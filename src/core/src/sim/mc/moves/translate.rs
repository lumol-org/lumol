// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

use rand::distributions::{Sample, Range};
use rand::Rng;

use std::usize;
use std::f64;

use super::MCMove;
use super::select_molecule;

use types::Vector3D;
use sys::{System, EnergyCache};

/// Monte Carlo move for translating a molecule
pub struct Translate {
    /// Type of molecule to translate. `None` means all molecules.
    moltype: Option<u64>,
    /// Index of the molecule to translate
    molid: usize,
    /// New positions of the atom in the translated molecule
    newpos: Vec<Vector3D>,
    /// Maximum displacement value
    delta: f64,
    /// The maximum value must not exceed this value, if set
    maximum_cutoff: Option<f64>,
    /// Translation range for random number generation
    range: Range<f64>,
}

impl Translate {
    /// Create a new `Translate` move, with maximum displacement of `delta`.
    /// Translating all the molecules in the system.
    pub fn new(delta: f64) -> Translate {
        Translate::create(delta, None)
    }

    /// Create a new `Translate` move, with maximum displacement of `delta`.
    /// Translating only molecules with `moltype` type.
    pub fn with_moltype(delta: f64, moltype: u64) -> Translate {
        Translate::create(delta, Some(moltype))
    }

    /// Factorizing the constructors
    fn create(delta: f64, moltype: Option<u64>) -> Translate {
        assert!(delta > 0.0, "delta must be positive in Translate move");
        let delta = delta / f64::sqrt(3.0);
        Translate {
            moltype: moltype,
            molid: usize::MAX,
            newpos: Vec::new(),
            delta: delta,
            maximum_cutoff: None,
            range: Range::new(-delta, delta),
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

    fn setup(&mut self, system: &System) {
        // Limit the displacement range to the maximum cutoff
        self.maximum_cutoff = system.maximum_cutoff();
        if let Some(max) = self.maximum_cutoff {
            if self.delta > max {
                warn!("Changing the maximal displacement for Translate, because the interactions \
                       cutoff is too low.");
                self.delta = max
            }
        }
    }

    fn prepare(&mut self, system: &mut System, rng: &mut Box<Rng>) -> bool {
        if let Some(id) = select_molecule(system, self.moltype, rng) {
            self.molid = id;
        } else {
            warn!("Can not translate molecule: no molecule of this type in the system.");
            return false;
        }

        // Create random displacement vector.
        let delta = Vector3D::new(self.range.sample(rng),
                                  self.range.sample(rng),
                                  self.range.sample(rng));

        // Generate displaced coordinates
        // Note that this may move a particles' center-of-mass (com) out of
        // the cell. If the move is accepted, we have to wrap the com such
        // that it lies inside the cell.
        let indexes = system.molecule(self.molid).iter();
        self.newpos = system.particles().position[indexes].to_vec();
        for newpos in &mut self.newpos {
            *newpos += delta;
        }
        return true;
    }

    fn cost(&self, system: &System, beta: f64, cache: &mut EnergyCache) -> f64 {
        let idxes = system.molecule(self.molid).iter().collect::<Vec<_>>();
        let cost = cache.move_particles_cost(system, idxes, &self.newpos);
        return cost * beta;
    }

    fn apply(&mut self, system: &mut System) {
        {
            // Update positions.
            let indexes = system.molecule(self.molid).iter();
            let positions = &mut system.particles_mut().position[indexes];
            for (position, newpos) in izip!(positions, &self.newpos) {
                *position = *newpos;
            }
        }
        // Move molecule such that its center-of-mass is inside the simulation
        // cell. Note that particles of the molecule may still be outside the
        // cell, but that is not important.
        system.wrap_molecule(self.molid)
    }

    fn restore(&mut self, _: &mut System) {
        // Nothing to do.
    }

    fn update_amplitude(&mut self, scaling_factor: Option<f64>) {
        if let Some(s) = scaling_factor {
            if let Some(max) = self.maximum_cutoff {
                if (self.delta * s) > max {
                    warn_once!("Tried to increase the maximum amplitude for translations to more \
                                than the maximum cutoff -- ignoring.");
                    return;
                }
            }

            self.delta *= s;
            self.range = Range::new(-self.delta, self.delta);
        };
    }
}
