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
use self::rand::Rng;

use universe::Universe;

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
    fn prepare(&mut self, universe: &mut Universe, rng: &mut Box<Rng>) -> bool;

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

/// Select a random molecule in the universe using `rng` as random number
/// generator. If `moltype` is `None`, any molecule can be choosen. If `moltype`
/// is `Some(molecule_type)`, then a molecule with matching type is selected.
///
/// This function returns `None` if no matching molecule was found, and
/// `Some(molid)` with `molid` the index of the molecule if a molecule was
/// selected.
fn select_molecule(universe: &Universe, moltype: Option<u64>, rng: &mut Box<Rng>) -> Option<usize> {
    if let Some(moltype) = moltype {
        // Pick a random molecule with matching moltype
        let mols = universe.molecules_with_moltype(moltype);
        return rng.choose(&mols).map(|id| *id);
    } else {
        let nmols = universe.molecules().len();
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
