// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

use std::io::prelude::*;
use std::io;
use std::fs::File;
use std::path::{Path, PathBuf};

use super::Output;
use utils;
use sys::System;

/// The `EnergyOutput` writes the energy of the system to a text file, organized
/// as: `PotentialEnergy     KineticEnergy     TotalEnergy`.
pub struct EnergyOutput {
    file: File,
    path: PathBuf
}

impl EnergyOutput {
    /// Create a new `EnergyOutput` writing to `filename`. The file is replaced
    /// if it already exists.
    pub fn new<P: AsRef<Path>>(filename: P) -> Result<EnergyOutput, io::Error> {
        Ok(EnergyOutput{
            file: try!(File::create(filename.as_ref())),
            path: filename.as_ref().to_owned(),
        })
    }
}

impl Output for EnergyOutput {
    fn setup(&mut self, _: &System) {
        if let Err(err) = writeln!(&mut self.file, "# Energy of the simulation (kJ/mol)") {
            fatal_error!("Could not write to file '{}': {}", self.path.display(), err);
        }
        if let Err(err) = writeln!(&mut self.file, "# Step Potential Kinetic Total") {
            fatal_error!("Could not write to file '{}': {}", self.path.display(), err);
        }
    }

    fn write(&mut self, system: &System) {
        let potential = utils::unit_to(system.potential_energy(), "kJ/mol");
        let kinetic = utils::unit_to(system.kinetic_energy(), "kJ/mol");
        let total = utils::unit_to(system.total_energy(), "kJ/mol");
        if let Err(err) = writeln!(&mut self.file, "{} {} {} {}", system.step(), potential, kinetic, total) {
            error!("Could not write to file '{}': {}", self.path.display(), err);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::tests::test_output;

    #[test]
    fn energy() {
        test_output(|path| {
            Box::new(EnergyOutput::new(path).unwrap())
        },
"# Energy of the simulation (kJ/mol)
# Step Potential Kinetic Total
0 1.5000000000000027 949.9201593348566 951.4201593348566
"
        );
    }
}
