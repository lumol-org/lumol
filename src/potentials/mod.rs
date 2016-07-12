// Cymbalum, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

//! Interaction potentials for energy and forces computations
//!
//! # Potential functions
//!
//! A potential function is defined as a parametric function giving the energy
//! and the force acting on a specific particle. These function are represented
//! by the [`PotentialFunction`][PotentialFunction] trait, and allow to compute
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
//! # use cymbalum::potentials::{PotentialFunction, PairPotential, DihedralPotential};
//! #[derive(Clone)]
//! struct OnePotential;
//!
//! // OnePotential is a potential function
//! impl PotentialFunction for OnePotential {
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
//! [PotentialFunction]: trait.PotentialFunction.html
//! [PairPotential]: trait.PairPotential.html
//! [AnglePotential]: trait.AnglePotential.html
//! [DihedralPotential]: trait.DihedralPotential.html
//! [GlobalPotential]: trait.GlobalPotential.html
//! [CoulombicPotential]: trait.CoulombicPotential.html

mod functions;
pub use self::functions::PotentialFunction;
pub use self::functions::{PairPotential, AnglePotential, DihedralPotential};
pub use self::functions::{NullPotential, LennardJones, Harmonic, CosineHarmonic};
pub use self::functions::Torsion;

mod computations;
pub use self::computations::Computation;
pub use self::computations::{TableComputation, CutoffComputation};

mod restrictions;
pub use self::restrictions::{PairRestriction, RestrictionInfo};

mod global;
pub use self::global::{GlobalPotential, CoulombicPotential};
pub use self::global::{Wolf, Ewald};
