// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

//! Molecular dynamics related algorithms.
//!
//! Each class of algorithms is represented by a `trait`, and each single method
//! is a `struct` implementing said `trait`.
//!
//! The main entry point is the
//! [`MolecularDynamics`](struct.MolecularDynamics.html) struct, implementing
//! the [`Propagator`](../trait.Propagator.html) trait.
//!
//! # Integrator
//!
//! The core of molecular dynamics simulation happens in an
//! [`Integrator`](trait.Integrator.html). The implemented integrators are the
//! following:
//!
//! - [`Verlet`](struct.Verlet.html): simple Verlet integrator;
//! - [`VelocityVerlet`](struct.VelocityVerlet.html): simple velocity-Verlet
//!   integrator;
//! - [`LeapFrog`](struct.LeapFrog.html): Leap-Frog integrator;
//! - [`BerendsenBarostat`](struct.BerendsenBarostat.html): isotropic Berendsen
//!   barostat coupled to a velocity-Verlet integrator;
//! - [`AnisoBerendsenBarostat`](struct.AnisoBerendsenBarostat.html) anisotropic
//!   Berendsen barostat coupled to a velocity-Verlet integrator
//!
//! # Themostats
//!
//! [`Thermostat`](trait.Thermostat.html) are algorihtms used to fix the
//! temperature during a simulation. The following implementation are available:
//!
//! - [`RescaleThermostat`](struct.RescaleThermostat.html): basic themostat
//!   rescaling the velocities of all atoms. This is highly unphysical but can
//!   be usefull for equilibration;
//! - [`CSVRThermostat`](struct.CSVRThermostat.html): Canonical Sampling through
//!   Velocities Rescaling is a well-behaved thermostating algorithm generating
//!   the expected canonical ensemble distribution of states.
//! - [`BerendsenThermostat`](struct.BerendsenThermostat.html): berendsen or
//!   weak-coupling thermostat;
//!
//! # Control
//!
//! [`Control`](trait.Control.html) algorihtms group any algorithm modifying the
//! system to enforce some type of invariant.
//!
//! - [`RemoveTranslation`](struct.RemoveTranslation.html): remove the global
//!   translational momentum of a system;
//! - [`RemoveRotation`](struct.RemoveRotation.html): remove the global
//!   rotational momentum of a system;
//! - [`Rewrap`](struct.Rewrap.html): wrap all atoms from a system inside the
//!   unit cell;

mod integrators;
pub use self::integrators::Integrator;

pub use self::integrators::AnisoBerendsenBarostat;
pub use self::integrators::BerendsenBarostat;
pub use self::integrators::LeapFrog;
pub use self::integrators::VelocityVerlet;
pub use self::integrators::Verlet;

mod controls;
pub use self::controls::Control;

pub use self::controls::RemoveRotation;
pub use self::controls::RemoveTranslation;
pub use self::controls::Rewrap;

mod thermostats;
pub use self::thermostats::Thermostat;

pub use self::thermostats::RescaleThermostat;
pub use self::thermostats::BerendsenThermostat;
pub use self::thermostats::CSVRThermostat;

mod molecular_dynamics;
pub use self::molecular_dynamics::MolecularDynamics;
