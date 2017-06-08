// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

//! Saving properties of a system during a simulation

use sys::System;

/// The `Output` trait define the interface for all the quantities outputted by
/// the simulation during the run. An Output can be a text or a binary data
/// file, an image, a text log, …
pub trait Output {
    /// Function called once at the beginning of the simulation, which allow
    /// for some setup of the output if needed.
    fn setup(&mut self, _: &System) {}

    /// Write the output from the system.
    fn write(&mut self, system: &System);

    /// Function called once at the end of the simulation.
    fn finish(&mut self, _: &System) {}
}

mod tests;
mod cell;
mod energy;
mod custom;
mod properties;
mod trajectory;

pub use self::cell::CellOutput;
pub use self::energy::EnergyOutput;
pub use self::custom::{CustomOutput, CustomOutputError};
pub use self::properties::PropertiesOutput;
pub use self::trajectory::TrajectoryOutput;
