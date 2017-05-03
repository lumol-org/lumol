// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 Lumol's contributors — BSD license

//! The `system` module provide a way to store data about a simulated system.

mod config;
pub use self::config::{PeriodicTable, ElementData};
pub use self::config::{Particle, ParticleKind};
pub use self::config::Composition;
pub use self::config::{UnitCell, CellShape};
pub use self::config::{Bond, Angle, Dihedral};
pub use self::config::Connectivity;
pub use self::config::{CONNECT_12, CONNECT_13, CONNECT_14, CONNECT_FAR};
pub use self::config::Molecule;
pub use self::config::molecule_type;
