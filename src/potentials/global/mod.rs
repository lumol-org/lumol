// Cymbalum, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

//! Global potential are potentials acting on the whole system at once
//!
//! They can be coulombic potentials, or external provided potential function
//! for example.
use system::System;
use types::{Matrix3, Vector3D};
use super::PairRestriction;

/// The `GlobalPotential` trait represent a potential acting on the whole
/// system at once.
pub trait GlobalPotential: GlobalCache + BoxCloneGlobal {
    /// Compute the energetic contribution of this potential
    fn energy(&self, system: &System) -> f64;
    /// Compute the force contribution of this potential
    fn forces(&self, system: &System) -> Vec<Vector3D>;
    /// Compute the virial contribution of this potential
    fn virial(&self, system: &System) -> Matrix3;
}

/// Energetic cache for global potentials. This trait have an opt-in default
/// implementation that you can use by implementing the `DefaultGlobalCache`
/// trait.
///
/// # Example
///
/// ```rust
/// struct MyGlobalPotential;
///
/// impl GlobalPotential for MyGlobalPotential {
///    ... // Implementation
/// }
///
/// // Use the default cache
/// impl DefaultGlobalCache for MyGlobalPotential {}
/// ```
pub trait GlobalCache {
    /// Get the cost of moving the `system` particles whose indexes are in
    /// `idxes` to `newpos`
    fn move_particles_cost(&self, system: &System, idxes: &[usize], newpos: &[Vector3D]) -> f64;
    /// Update the cache after a call to a `*_cost` function.
    fn update(&mut self);
}

/// Marker trait for the default implementation of GlobalCache. This default
/// implementation is slow, as it create a copy of the system
pub trait DefaultGlobalCache {}

impl<T> GlobalCache for T where T: DefaultGlobalCache + GlobalPotential {
    fn move_particles_cost(&self, system: &System, idxes: &[usize], newpos: &[Vector3D]) -> f64 {
        let mut system = system.clone();
        let old_e = self.energy(&system);
        for (i, &pi) in idxes.iter().enumerate() {
            system[pi].position = newpos[i];
        }
        let new_e = self.energy(&system);
        return new_e - old_e;
    }

    fn update(&mut self) {
        // Nothing to do
    }
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
