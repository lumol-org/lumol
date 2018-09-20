// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

//! Types and traits for representing simulation algorithms
mod propagator;
pub use self::propagator::DegreesOfFreedom;
pub use self::propagator::Propagator;
pub use self::propagator::TemperatureStrategy;

pub mod md;
pub mod mc;
pub mod min;

mod simulations;
pub use self::mc::MonteCarlo;
pub use self::md::MolecularDynamics;
pub use self::min::Minimization;
pub use self::simulations::Simulation;

mod utils;
pub use self::utils::Alternator;
