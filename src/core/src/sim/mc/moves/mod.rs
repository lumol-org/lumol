// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

//! A Monte-Carlo simulation is based on a set of applicable moves.
//!
//! For example, NVT Monte-Carlo will use the `Translate` move only for particles
//! and add the `Rotate` moves for molecules. NPT Monte-Carlo will add the
//! `VolumeResize` move.
//!
//! In all this module, beta refers to the Boltzmann factor 1/(kB T)
use rand::Rng;

use sys::{System, EnergyCache};

/// The `MCMove` trait correspond to the set of methods used in Monte-Carlo
/// simulations.
pub trait MCMove {
    /// Give a short description of this move
    fn describe(&self) -> &str;

    /// Prepare the move, by selecting the particles to move, and the parameters
    /// of the move. The `rng` random number generator should be used to
    /// generate the parameters of the move.
    ///
    /// This function should return true is we can perform the move, and false
    /// otherwise.
    fn prepare(&mut self, system: &mut System, rng: &mut Box<Rng>) -> bool;

    /// Get the cost of performing this move on `system`. For example in
    /// simple NVT simulations, this cost is the energetic difference over
    /// `beta`. The cost must be dimensionless, and will be placed in an
    /// exponential. The `cache` should be used to compute the cost, or the
    /// `cache.unused` function should be used to ensure that the cache is
    /// updated  as needed after this move.
    fn cost(&self, system: &System, beta: f64, cache: &mut EnergyCache) -> f64;

    /// Apply the move, if it has not already been done in `prepare`.
    fn apply(&mut self, system: &mut System);

    /// Restore the system to it's initial state, if it has been changed in
    /// `prepare`.
    fn restore(&mut self, system: &mut System);
}

/// Select a random molecule in the system using `rng` as random number
/// generator. If `moltype` is `None`, any molecule can be chosen. If `moltype`
/// is `Some(molecule_type)`, then a molecule with matching type is selected.
///
/// This function returns `None` if no matching molecule was found, and
/// `Some(molid)` with `molid` the index of the molecule if a molecule was
/// selected.
fn select_molecule(system: &System, moltype: Option<u64>, rng: &mut Box<Rng>) -> Option<usize> {
    if let Some(moltype) = moltype {
        // Pick a random molecule with matching moltype
        let mols = system.molecules_with_moltype(moltype);
        return rng.choose(&mols).cloned();
    } else {
        let nmols = system.molecules().len();
        if nmols == 0 {
            return None;
        } else {
            return Some(rng.gen_range(0, nmols));
        }
    }
}

mod translate;
pub use self::translate::Translate;

mod rotate;
pub use self::rotate::Rotate;
