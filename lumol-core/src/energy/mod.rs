// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

//! Interaction potentials for energy and forces computations
//!
//! # Potentials
//!
//! In classical simulations, the total energy of a system is separated as a
//! sum containing terms from the non-bonded pairs, the bonds, the angles, the
//! dihedral angles in the system; and the electrostatic interactions.
//!
//! Potentials are used to compute the interaction energy in a system. They are
//! represented by a [`Potential`][Potential] trait, used to compute the force
//! and the energy of interaction. In order to add a new potential to lumol,
//! one has to implement the Potential trait, and then indicate how the
//! potential can be used. This is done by implementing one or more of the
//! potentials marker traits:
//!
//! - [`PairPotential`][PairPotential] for non-bonded two body interactions;
//! - [`BondPotential`][BondPotential] for covalent bonds interactions;
//! - [`AnglePotential`][AnglePotential] for covalent angles interactions;
//! - [`DihedralPotential`][DihedralPotential] for covalent dihedral angles
//!   interactions.
//!
//! ```
//! use lumol_core::energy::{Potential, PairPotential, DihedralPotential};
//!
//! #[derive(Clone)]
//! struct OnePotential;
//!
//! // OnePotential is a potential
//! impl Potential for OnePotential {
//!     fn energy(&self, _: f64) -> f64 {1.0}
//!     fn force(&self, _: f64) -> f64 {0.0}
//! }
//!
//! // It can be used for pair and dihedral potentials, but not for angles or
//! // bonds.
//! impl PairPotential for OnePotential {
//!     /* some code omitted */
//!     # fn tail_energy(&self, cutoff: f64) -> f64 {0.0}
//!     # fn tail_virial(&self, cutoff: f64) -> f64 {0.0}
//! }
//! impl DihedralPotential for OnePotential {}
//! ```
//!
//! # Global and Coulombic potentials
//!
//! Global potentials are potentials that have an effect on the whole system at
//! once. They are defined by implementing the [`GlobalPotential`]
//! [GlobalPotential] trait. [`CoulombicPotential`][CoulombicPotential] are a
//! specific version of global potentials used to compute electrostatic
//! interactions.
//!
//! [Potential]: trait.Potential.html
//! [PairPotential]: trait.PairPotential.html
//! [BondPotential]: trait.BondPotential.html
//! [AnglePotential]: trait.AnglePotential.html
//! [DihedralPotential]: trait.DihedralPotential.html
//! [GlobalPotential]: trait.GlobalPotential.html
//! [CoulombicPotential]: trait.CoulombicPotential.html
use crate::{Matrix3, Vector3D};

/// A potential for force and energy computations.
///
/// A potential is defined with two functions that takes a single scalar
/// variable and return the corresponding energy or norm of the force. The
/// scalar variable will be the distance for pair potentials, the angle for
/// angles or dihedral angles potentials, *etc.*
///
/// # Example
///
/// ```
/// # use std::f64;
/// use lumol_core::energy::Potential;
///
/// /// An hard sphere potential
/// #[derive(Clone)]
/// struct HardSphere {
///     /// Sphere radius
///     pub r: f64
/// }
///
/// impl Potential for HardSphere {
///     fn energy(&self, x: f64) -> f64 {
///         if x < self.r {
///             f64::INFINITY
///         } else {
///             0.0
///         }
///     }
///
///     fn force(&self, x: f64) -> f64 {
///         0.0
///     }
/// }
/// ```
pub trait Potential: Sync + Send {
    /// Get the energy corresponding to the variable `x`
    fn energy(&self, x: f64) -> f64;
    /// Get the force norm corresponding to the variable `x`
    fn force(&self, x: f64) -> f64;
}

/// Marker trait for potentials that can be used for non-bonded two body
/// interactions.
///
/// # Example
///
/// ```
/// use lumol_core::energy::{Potential, PairPotential};
///
/// // A no-op potential
/// #[derive(Clone)]
/// struct Null;
///
/// impl Potential for Null {
///     fn energy(&self, x: f64) -> f64 {0.0}
///     fn force(&self, x: f64) -> f64 {0.0}
/// }
///
/// // By implementing this trait, we can use the Null potential for pair
/// // interactions
/// impl PairPotential for Null {
///     fn tail_energy(&self, cutoff: f64) -> f64 {
///         return 0.0;
///     }
///
///     fn tail_virial(&self, cutoff: f64) -> f64 {
///         return 0.0;
///     }
/// }
/// ```
pub trait PairPotential: Potential + BoxClonePair {
    /// Compute the virial contribution corresponding to the distance `r`
    /// between the particles.
    fn virial(&self, r: &Vector3D) -> Matrix3 {
        let fact = self.force(r.norm());
        let rn = r.normalized();
        let force = fact * rn;
        force.tensorial(r)
    }

    /// Compute the tail correction to the energy for the given cutoff.
    ///
    /// Calling `V(r)` the `Potential::energy(r)` function corresponding to this
    /// potential, this function should return the integral from `cutoff` to
    /// infinity of `r^2 V(r)`: `\int_{cutoff}^\infty r^2 V(r) dr`.
    ///
    /// If this integral does not converge for the current potential, this
    /// function should then return 0 to disable tail corrections.
    fn tail_energy(&self, cutoff: f64) -> f64;

    /// Compute the tail correction to the virial for the given cutoff.
    ///
    /// Calling `f(r)` the `Potential::force(r)` function corresponding to this
    /// potential, this function should return the integral from `cutoff` to
    /// infinity of `f(r) r^3`: `\int_{cutoff}^\infty r^3 f(r) dr`.
    ///
    /// If this integral does not converge for the current potential, this
    /// function should then return 0.0 to disable tail corrections.
    fn tail_virial(&self, cutoff: f64) -> f64;
}
impl_box_clone!(PairPotential, BoxClonePair, box_clone_pair);

/// Marker trait for potentials that can be used for molecular bonds.
///
/// # Example
///
/// ```
/// use lumol_core::energy::{Potential, BondPotential};
///
/// // A no-op potential
/// #[derive(Clone)]
/// struct Null;
///
/// impl Potential for Null {
///     fn energy(&self, x: f64) -> f64 {0.0}
///     fn force(&self, x: f64) -> f64 {0.0}
/// }
///
/// // Now we can use the Null potential for bonds
/// impl BondPotential for Null {}
/// ```
pub trait BondPotential: Potential + BoxCloneBond {
    /// Compute the virial contribution corresponding to the distance `r`
    /// between the particles.
    fn virial(&self, r: &Vector3D) -> Matrix3 {
        let fact = self.force(r.norm());
        let rn = r.normalized();
        let force = fact * rn;
        force.tensorial(r)
    }
}
impl_box_clone!(BondPotential, BoxCloneBond, box_clone_bond);

/// Marker trait for potentials that can be used for molecular angles.
///
/// # Example
///
/// ```
/// use lumol_core::energy::{Potential, AnglePotential};
///
/// // A no-op potential
/// #[derive(Clone)]
/// struct Null;
///
/// impl Potential for Null {
///     fn energy(&self, x: f64) -> f64 {0.0}
///     fn force(&self, x: f64) -> f64 {0.0}
/// }
///
/// // Now we can use the Null potential for angles
/// impl AnglePotential for Null {}
/// ```
pub trait AnglePotential: Potential + BoxCloneAngle {}
impl_box_clone!(AnglePotential, BoxCloneAngle, box_clone_angle);

/// Marker trait for potentials that can be used for molecular dihedral angles.
///
/// # Example
///
/// ```
/// use lumol_core::energy::{Potential, DihedralPotential};
///
/// // A no-op potential
/// #[derive(Clone)]
/// struct Null;
///
/// impl Potential for Null {
///     fn energy(&self, x: f64) -> f64 {0.0}
///     fn force(&self, x: f64) -> f64 {0.0}
/// }
///
/// // Now we can use the Null potential for dihedral angles
/// impl DihedralPotential for Null {}
/// ```
pub trait DihedralPotential: Potential + BoxCloneDihedral {}
impl_box_clone!(DihedralPotential, BoxCloneDihedral, box_clone_dihedral);

mod functions;
pub use self::functions::{BornMayerHuggins, Buckingham, Gaussian, Morse, Torsion};
pub use self::functions::{CosineHarmonic, Harmonic, LennardJones, NullPotential};
pub use self::functions::Mie;

mod computations;
pub use self::computations::{Computation, TableComputation};

mod restrictions;
pub use self::restrictions::{PairRestriction, RestrictionInfo, BondPath};

mod global;
pub use self::global::{CoulombicPotential, GlobalCache, GlobalPotential};
pub use self::global::{Ewald, SharedEwald, Wolf};

mod pairs;
pub use self::pairs::PairInteraction;
