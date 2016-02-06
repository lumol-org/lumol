/* Cymbalum, Molecular Simulation in Rust - Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
 */
//! Global potential are potentials acting on the whole system at once
//!
//! They can be coulombic potentials, or external provided potential function
//! for example.
use universe::Universe;
use types::{Matrix3, Vector3D};
use super::PairRestriction;

/// The `GlobalPotential` trait represent a potential acting on the whole
/// universe at once.
pub trait GlobalPotential: BoxCloneGlobal {
    /// Compute the energetic contribution of this potential
    fn energy(&self, universe: &Universe) -> f64;
    /// Compute the force contribution of this potential
    fn forces(&self, universe: &Universe) -> Vec<Vector3D>;
    /// Compute the virial contribution of this potential
    fn virial(&self, universe: &Universe) -> Matrix3;
}

/// Electrostatic potential solver should implement the `CoulombicPotential`
/// trait.
pub trait CoulombicPotential : GlobalPotential + BoxCloneCoulombic {
    /// Set the restriction scheme to use to `restriction`. All future call to
    /// `energy`, `force` or `virial` should use this restriction.
    fn set_restriction(&mut self, restriction: PairRestriction);
}

impl_box_clone!(GlobalPotential, BoxCloneGlobal, box_clone_gobal);
impl_box_clone!(CoulombicPotential, BoxCloneCoulombic, box_clone_coulombic);

mod wolf;
pub use self::wolf::Wolf;

mod ewald;
pub use self::ewald::Ewald;
