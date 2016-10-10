// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

//! Interaction potentials for energy and forces computations
//!
//! # Potential functions
//!
//! A potential function is defined as a parametric function giving the energy
//! and the force acting on a specific particle. These function are represented
//! by the [`Potential`][Potential] trait, and allow to compute
//! the energy or force norm corresponding to a scalar variable. This scalar
//! variable will be the distance for non-bonded pair interactions; the angle
//! for covalent angles interactions, *etc.*
//!
//! Additional marker traits are used to specify wheter a given potential can
//! be used for a specific interaction type. These marker traits are:
//!
//! - [`PairPotential`][PairPotential] for bonded and non-bonded two body
//!   interactions;
//! - [`AnglePotential`][AnglePotential] for covalent angles interactions;
//! - [`DihedralPotential`][DihedralPotential] for covalent dihedral angles
//!   interactions.
//!
//! To be able to use a specific potential in one of these contexts, the only
//! thing to do is to implement the corresponding trait.
//!
//! ```
//! # use lumol::potentials::{Potential, PairPotential, DihedralPotential};
//! #[derive(Clone)]
//! struct OnePotential;
//!
//! // OnePotential is a potential function
//! impl Potential for OnePotential {
//!     fn energy(&self, _: f64) -> f64 {1.0}
//!     fn force(&self, _: f64) -> f64 {0.0}
//! }
//!
//! // It can be used for pair and dihedral potentials, but not for angles
//! impl PairPotential for OnePotential {}
//! impl DihedralPotential for OnePotential {}
//! ```
//!
//! # Global potentials
//!
//! Global potentials are potentials that have an effect on the whole system at
//! once. They are defined by implementing the [`GlobalPotential`][GlobalPotential]
//! trait and giving methods for energy, forces and virial contributions.
//!
//! [`CoulombicPotential`][CoulombicPotential] are a specific version of global
//! potentials, that are used to compute charges-charges interactions.
//!
//! [Potential]: trait.Potential.html
//! [PairPotential]: trait.PairPotential.html
//! [AnglePotential]: trait.AnglePotential.html
//! [DihedralPotential]: trait.DihedralPotential.html
//! [GlobalPotential]: trait.GlobalPotential.html
//! [CoulombicPotential]: trait.CoulombicPotential.html
use types::{Matrix3, Vector3D};

/// A set of two parametric functions which takes a single scalar variable and
/// return the corresponding energy or norm of the force.
///
/// The scalar variable will be the distance for pair potentials, the angle for
/// angles or dihedral angles potentials, *etc.*
pub trait Potential : Sync + Send + BoxClonePotential {
    /// Get the energy corresponding to the variable `x`
    fn energy(&self, x: f64) -> f64;
    /// Get the force norm corresponding to the variable `x`
    fn force(&self, x: f64) -> f64;
}

impl_box_clone!(Potential, BoxClonePotential, box_clone_potential);

/// Computation of virial contribution for a potential. The provided fucntion
/// only apply to two-body virial contributions.
pub trait Virial: Potential {
    /// Compute the virial contribution corresponding to the distance `r`
    /// between the particles
    fn virial(&self, r: &Vector3D) -> Matrix3 {
        let fact = self.force(r.norm());
        let rn = r.normalized();
        let force = fact * rn;
        force.tensorial(r)
    }
}

impl<T: Potential> Virial for T {}

/// Potential that can be used for non-bonded two body interactions
pub trait PairPotential : Virial + BoxClonePair {}
impl_box_clone!(PairPotential, BoxClonePair, box_clone_pair);

/// Potential that can be used for molecular bonds
pub trait BondPotential : Virial + BoxCloneBond {}
impl_box_clone!(BondPotential, BoxCloneBond, box_clone_bond);

/// Potential that can be used for molecular angles.
pub trait AnglePotential : Potential + BoxCloneAngle {}
impl_box_clone!(AnglePotential, BoxCloneAngle, box_clone_angle);

/// Potential that can be used for molecular dihedral angles.
pub trait DihedralPotential : Potential + BoxCloneDihedral {}
impl_box_clone!(DihedralPotential, BoxCloneDihedral, box_clone_dihedral);

mod functions;
pub use self::functions::{NullPotential, LennardJones, Harmonic, CosineHarmonic};
pub use self::functions::Torsion;

mod computations;
pub use self::computations::{Computation, TableComputation};

mod restrictions;
pub use self::restrictions::{PairRestriction, RestrictionInfo};

mod global;
pub use self::global::{GlobalPotential, CoulombicPotential};
pub use self::global::{Wolf, Ewald};

mod pairs;
pub use self::pairs::PairInteraction;
