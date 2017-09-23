// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

use rand::distributions::{Sample, Range};
use rand::Rng;

use std::usize;

use super::{MCMove, select_molecule};
use super::configurational_bias::*;

use consts::K_BOLTZMANN;
use types::Vector3D;
use sys::{System, Molecule, EnergyCache};

/// Monte-Carlo move for rebuilding a molecule
///
/// TODO:
/// - dual cutoff scheme
/// - cutoff between molecules
/// - early rejection scheme
/// - handle first particle differently
pub struct Rebuild {
    /// Type of molecule to Rebuild. `None` means all molecules.
    moltype: Option<u64>,
    /// Index of the molecule to rebuild
    molid: usize,
    /// New positions of the rebuilt molecule
    newpos: Vec<Vector3D>,
    /// Rosenbluth weight of the new configuration
    new_weight: f64,
    /// Rosenbluth weight of the old configuration
    old_weight: f64,
    /// number of configurational bias steps
    cb_steps: u32,
    // /// early rejection tolerance, hardcode?
    //early_rejection_tolerance: f64,
}

impl Rebuild {
    /// Create a new `Rebuild` move, using `k`configurational bias steps.
    pub fn new(k: u32) -> Rebuild {
        Rebuild::create(k, None)
    }

    /// Create a new `Rebuild` move, using `k`configurational bias steps.
    /// Rebuilding only molecules of type `moltype`.
    pub fn with_moltype(k: u32, moltype: u64) -> Rebuild {
        Rebuild::create(k, Some(moltype))
    }

    /// Factorizing the constructors
    fn create(cb_steps: u32, moltype: Option<u64>) -> Rebuild {
        Rebuild {
            moltype: None,
            molid: usize::MAX,
            newpos: Vec::new(),
            new_weight: 1.0,
            old_weight: 1.0,
            cb_steps: cb_steps,
            // early_rejection_tolerance: 1.0e-200,
        }
    }
}

impl Default for Rebuild {
    fn default() -> Rebuild {
        Rebuild::new(1)
    }
}

impl MCMove for Rebuild {
    fn describe(&self) -> &str {
        "rebuilding a molecule"
    }

    fn setup(&mut self, _: &System) {
    }

    fn prepare(&mut self, system: &mut System, rng: &mut Box<Rng>) -> bool {
        if let Some(id) = select_molecule(system, self.moltype, rng) {
            self.molid = id;
        } else {
            warn!("Cannot rebuild molecule: no molecule of this type in the system.");
            return false;
        }

        // clear
        self.newpos.clear();
        self.new_weight = 1.0;
        self.old_weight = 1.0;

        let mut trial_weights: Vec<f64> = Vec::with_capacity(self.cb_steps as usize);
        let mut trial_boltzmann_factor: Vec<f64> = Vec::with_capacity(self.cb_steps as usize);
        let mut trialpos: Vec<Vector3D> = Vec::with_capacity(self.cb_steps as usize);

        let beta = 1.0 / (K_BOLTZMANN * system.temperature());
        let molecule = system.molecule(self.molid);
        let mut oldpos: Vec<Vector3D> = Vec::new();
        for pi in molecule.iter() {
            oldpos.push(system[pi].position);
        }

        // loop over each atom (segment) of the molecule
        // .iter() returns the index of the atom to use to index into system
        // pid == particle index
        for (i, pid) in molecule.iter().enumerate() {
            trial_weights.clear();
            trial_boltzmann_factor.clear();
            trialpos.clear();

            // loop over trials
            // this is how you'd find it in literature
            // maybe it is better to lift loop into `trial_position` and computation of the Boltzmann factor?
            // that would save us all (but one) lookup in trial_position
            // i.e write:
            // fn trial_positions(system, pid, self.newpos, beta, rng) -> &Vec<Vector3D>
            // where the resulting array would contain cb_steps elements
            for step in 0..self.cb_steps {
                // new weight
                let new_trial = trial_position(system, molecule, &self.newpos, beta, rng);
                // energy between the trial particle and ALL other particles
                // (excluding intramolecular energy wrt. the growing molecule)
                let energy = trial_non_covalent_energy(system, &self.newpos, &new_trial, &molecule);
                trialpos.push(new_trial);
                trial_boltzmann_factor.push(f64::exp(-beta * energy));
            }
            let w = trial_boltzmann_factor.iter().sum();
            for factor in trial_boltzmann_factor.iter() {
                trial_weights.push(factor / w)
            }
            self.new_weight *= w; // do we need an array here or could we just build the product here?
            self.newpos.push(select_position(&trialpos, &trial_weights, rng));

            // compute weight for OLD configuration
            // the zero'th weight is the old position
            trial_weights.clear();
            trial_boltzmann_factor.clear();
            trialpos.clear();

            let energy = trial_non_covalent_energy(system, &self.newpos, &oldpos[i], &molecule);
            trial_boltzmann_factor.push(f64::exp(-beta * energy));
            trialpos.push(oldpos[i]);
            for step in 1..self.cb_steps {
                // new weight
                let new_trial = trial_position(system, &molecule, &oldpos[0..i], beta, rng);;
                let energy = trial_non_covalent_energy(system, &oldpos[0..i], &new_trial, &molecule);
                trial_boltzmann_factor.push(f64::exp(-beta * energy));
            }
            let w = trial_boltzmann_factor.iter().sum();
            for factor in trial_boltzmann_factor.iter() {
                trial_weights.push(factor / w)
            }
            self.old_weight *= w
        }
        true
    }

    fn cost(&self, system: &System, beta: f64, cache: &mut EnergyCache) -> f64 {
        let idxes = system.molecule(self.molid).iter().collect::<Vec<_>>();
        // update cache
        let _ = cache.move_particles_cost(system, idxes, &self.newpos);
        return beta * (self.new_weight / self.old_weight);
    }

    fn apply(&mut self, system: &mut System) {
        for (i, pi) in system.molecule(self.molid).iter().enumerate() {
            system[pi].position = self.newpos[i];
        }
    }

    fn restore(&mut self, _: &mut System) {
        // Nothing to do.
    }

    fn update_amplitude(&mut self, _: Option<f64>) {
        // We could check the scaling and increase the number of steps
        // if the scaling is larger than one or decrease if its lower.
    }
}
