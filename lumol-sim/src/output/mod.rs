// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

//! Saving properties of a system during a simulation

use lumol_core::System;

/// The `Output` trait defines the interface for all the quantities outputted by
/// the simulation during the run. An Output can be a text or a binary data
/// file, an image, a text log, …
pub trait Output {
    /// Function called once at the beginning of the simulation, which allows
    /// for some setup of the output if needed.
    fn setup(&mut self, _: &System) {}

    /// Write the output from the system.
    fn write(&mut self, system: &System);

    /// Function called once at the end of the simulation.
    fn finish(&mut self, _: &System) {}
}

mod tests;

macro_rules! writeln_or_log {
    ($this: expr, $($args: expr),* $(,)*) => (
        if let Err(err) = writeln!(&mut $this.file, $($args,)*) {
            error!("could not write to file '{}': {}", $this.path.display(), err);
            return;
        }
    );
}

mod cell;
pub use self::cell::CellOutput;

mod stress;
pub use self::stress::StressOutput;

mod energy;
pub use self::energy::EnergyOutput;

mod custom;
pub use self::custom::{CustomOutput, CustomOutputError};

mod forces;
pub use self::forces::ForcesOutput;

mod properties;
pub use self::properties::PropertiesOutput;

mod trajectory;
pub use self::trajectory::TrajectoryOutput;
