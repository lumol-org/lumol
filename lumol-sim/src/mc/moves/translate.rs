// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

use rand::RngCore;
use rand_distr::{Uniform, Distribution};

use std::collections::BTreeSet;
use std::f64;
use std::usize;

use log::warn;
use log_once::warn_once;

use soa_derive::soa_zip;

use super::{MCDegreeOfFreedom, MCMove};
use super::select_molecule;

use lumol_core::{EnergyCache, System, MoleculeHash, Vector3D};

/// Monte Carlo move for translating a molecule
pub struct Translate {
    /// Hash of molecule to translate. `None` means all molecules.
    hash: Option<MoleculeHash>,
    /// Index of the molecule to translate
    molid: usize,
    /// New positions of the atom in the translated molecule
    newpos: Vec<Vector3D>,
    /// Maximum displacement value
    delta: f64,
    /// The maximum value must not exceed this value, if set
    maximum_cutoff: Option<f64>,
    /// Translation range for random number generation
    range: Uniform<f64>,
}

impl Translate {
    /// Create a new `Translate` move, with maximum displacement of `delta`.
    /// This move will apply to the molecules with the given `hash`, or all
    /// molecules if `hash` is `None`.
    pub fn new<H: Into<Option<MoleculeHash>>>(delta: f64, hash: H) -> Translate {
        assert!(delta > 0.0, "delta must be positive in Translate move");
        let delta = delta / f64::sqrt(3.0);
        Translate {
            hash: hash.into(),
            molid: usize::max_value(),
            newpos: Vec::new(),
            delta: delta,
            maximum_cutoff: None,
            range: Uniform::new(-delta, delta),
        }
    }
}

impl MCMove for Translate {
    fn describe(&self) -> &str {
        "molecular translation"
    }

    fn degrees_of_freedom(&self) -> MCDegreeOfFreedom {
        match self.hash {
            Some(hash) => {
                let mut all = BTreeSet::new();
                let _ = all.insert(hash);
                MCDegreeOfFreedom::Molecules(all)
            }
            None => MCDegreeOfFreedom::AllMolecules,
        }
    }

    fn setup(&mut self, system: &System) {
        // Limit the displacement range to the maximum cutoff
        self.maximum_cutoff = system.maximum_cutoff();
        if let Some(max) = self.maximum_cutoff {
            if self.delta > max {
                warn!(
                    "Changing the maximal displacement for Translate, \
                     because the interactions cutoff is too low."
                );
                self.delta = max;
            }
        }
    }

    fn prepare(&mut self, system: &mut System, rng: &mut dyn RngCore) -> bool {
        if let Some(id) = select_molecule(system, self.hash, rng) {
            self.molid = id;
        } else {
            warn!("Can not translate molecule: no molecule of this type in the system.");
            return false;
        }

        // Create random displacement vector.
        let delta = Vector3D::new(
            self.range.sample(rng),
            self.range.sample(rng),
            self.range.sample(rng)
        );

        // Generate displaced coordinates
        // Note that this may move a particles' center-of-mass (com) out of
        // the cell. If the move is accepted, we have to wrap the com such
        // that it lies inside the cell.
        self.newpos = system.molecule(self.molid).particles().position.to_vec();
        for newpos in &mut self.newpos {
            *newpos += delta;
        }
        return true;
    }

    fn cost(&self, system: &System, beta: f64, cache: &mut EnergyCache) -> f64 {
        return beta * cache.move_molecule_cost(system, self.molid, &self.newpos);
    }

    fn apply(&mut self, system: &mut System) {
        let cell = system.cell;
        let mut molecule = system.molecule_mut(self.molid);
        for (position, newpos) in soa_zip!(molecule.particles_mut(), [mut position], &self.newpos) {
            *position = *newpos;
        }

        // Move molecule such that its center-of-mass is inside the simulation
        // cell. Note that particles of the molecule may still be outside the
        // cell, but that is not important.
        molecule.wrap(&cell);
    }

    fn restore(&mut self, _: &mut System) {
        // Nothing to do.
    }

    fn update_amplitude(&mut self, scaling_factor: Option<f64>) {
        if let Some(s) = scaling_factor {
            if let Some(max) = self.maximum_cutoff {
                if (self.delta * s) > max {
                    warn_once!(
                        "Tried to increase the maximum amplitude for translations \
                         to more than the maximum cutoff -- ignoring."
                    );
                    return;
                }
            }

            self.delta *= s;
            self.range = Uniform::new(-self.delta, self.delta);
        };
    }
}
