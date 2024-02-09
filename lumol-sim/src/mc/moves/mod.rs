// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

//! A Monte Carlo simulation is based on a set of applicable moves.
//!
//! For example, NVT Monte Carlo will use the `Translate` move only for particles
//! and add the `Rotate` moves for molecules. NPT Monte Carlo will add the
//! `VolumeResize` move.
//!
//! In all this module, beta refers to the Boltzmann factor 1/(kB T)
use rand::{RngCore, Rng, seq::SliceRandom};
use std::collections::BTreeSet;
use lumol_core::{EnergyCache, System, MoleculeHash};

/// Possible degrees of freedom simulated by a given Monte Carlo move
#[derive(Clone, PartialEq, Debug)]
pub enum MCDegreeOfFreedom {
    /// All molecules are simulated
    AllMolecules,
    /// All molecules with a molecule type in the `BTreeSet` are simulated
    Molecules(BTreeSet<MoleculeHash>),
    /// All the particles are simulated
    Particles,
}

impl MCDegreeOfFreedom {
    /// Combine the degrees of freedom represented by this `MCDegreeOfFreedom`
    /// and `other`
    pub fn combine(self, other: MCDegreeOfFreedom) -> MCDegreeOfFreedom {
        use super::MCDegreeOfFreedom as DOF;
        match (self, other) {
            (DOF::Particles, _) | (_, DOF::Particles) => DOF::Particles,
            (DOF::AllMolecules, _) | (_, DOF::AllMolecules) => DOF::AllMolecules,
            (DOF::Molecules(mut all), DOF::Molecules(others)) => {
                all.extend(others);
                DOF::Molecules(all)
            }
        }
    }
}

/// The `MCMove` trait correspond to the set of methods used in Monte Carlo
/// simulations.
pub trait MCMove {
    /// Give a short description of this move
    fn describe(&self) -> &str;

    /// Set up move before simulation is run
    fn setup(&mut self, system: &System);

    /// Get the number of degrees of freedom simulated by this move
    fn degrees_of_freedom(&self) -> MCDegreeOfFreedom;

    /// Prepare the move by selecting the particles to move, and the parameters
    /// of the move. The `rng` random number generator should be used to
    /// generate the parameters of the move.
    ///
    /// This function should return true is we can perform the move, and false
    /// otherwise.
    fn prepare(&mut self, system: &mut System, rng: &mut dyn RngCore) -> bool;

    /// Get the cost of performing this move on `system`. For example in
    /// simple NVT simulations, this cost is the energetic difference between
    /// the new and the old state times beta. The cost must be dimmensionless.
    ///
    /// Note that the cost will be placed in an exponential with a negative sign.
    /// For NVT using the Metropolis criterion:
    /// cost = beta*(U_new - U_old) -> P_acc = min[1, exp(-cost)].
    ///
    /// The `cache` should be used to compute the cost, or the
    /// `cache.unused` function should be used to ensure that the cache is
    /// updated as needed after this move.
    fn cost(&self, system: &System, beta: f64, cache: &mut EnergyCache) -> f64;

    /// Apply the move, if it has not already been done in `prepare`.
    fn apply(&mut self, system: &mut System);

    /// Restore the system to it's initial state if it has been changed in
    /// `prepare`.
    fn restore(&mut self, system: &mut System);

    /// Update the sample range for displacements.
    fn update_amplitude(&mut self, scaling_factor: Option<f64>);
}

/// Select a random molecule in the system using `rng` as random number
/// generator. If `hash` is `None`, any molecule can be chosen. If `hash` is
/// `Some(hash)`, then a molecule with matching hash is selected.
///
/// This function returns `None` if no matching molecule was found, and
/// `Some(molid)` with `molid` the index of the molecule if a molecule was
/// selected.
fn select_molecule(system: &System, hash: Option<MoleculeHash>, rng: &mut dyn RngCore) -> Option<usize> {
    if let Some(hash) = hash {
        // Pick a random molecule with matching molecule type
        let mols = system.molecules()
            .enumerate()
            .filter(|(_, m)| m.hash() == hash)
            .map(|(i, _)| i)
            .collect::<Vec<_>>();

        mols.choose(rng).copied()
    } else {
        let molecules_count = system.molecules().count();
        if molecules_count == 0 {
            None
        } else {
            Some(rng.gen_range(0..molecules_count))
        }
    }
}

mod translate;
pub use self::translate::Translate;

mod rotate;
pub use self::rotate::Rotate;

mod resize;
pub use self::resize::Resize;
