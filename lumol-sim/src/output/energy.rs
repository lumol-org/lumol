// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

use std::fs::File;
use std::io::{self, BufWriter};
use std::io::prelude::*;
use std::path::{Path, PathBuf};

use log::error;

use super::Output;
use lumol_core::System;
use lumol_core::units;

/// The `EnergyOutput` writes the energy of the system to a text file, organized
/// as: `steps PotentialEnergy KineticEnergy TotalEnergy`.
pub struct EnergyOutput {
    file: BufWriter<File>,
    path: PathBuf,
}

impl EnergyOutput {
    /// Create a new `EnergyOutput` writing to `filename`. The file is replaced
    /// if it already exists.
    pub fn new<P: AsRef<Path>>(filename: P) -> Result<EnergyOutput, io::Error> {
        Ok(EnergyOutput {
            file: BufWriter::new(File::create(filename.as_ref())?),
            path: filename.as_ref().to_owned(),
        })
    }
}

impl Output for EnergyOutput {
    fn setup(&mut self, _: &System) {
        writeln_or_log!(self, "# Energy of the simulation (kJ/mol)");
        writeln_or_log!(self, "# Step Potential Kinetic Total");
    }

    fn write(&mut self, system: &System) {
        let potential = units::to(system.potential_energy(), "kJ/mol").expect("bad unit");
        let kinetic = units::to(system.kinetic_energy(), "kJ/mol").expect("bad unit");
        let total = units::to(system.total_energy(), "kJ/mol").expect("bad unit");
        writeln_or_log!(self, "{} {} {} {}", system.step, potential, kinetic, total);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::tests::test_output;

    #[test]
    fn energy() {
        test_output(
            |path| Box::new(EnergyOutput::new(path).unwrap()),
            "# Energy of the simulation (kJ/mol)
            # Step Potential Kinetic Total
            42 1.5000000000000027 949.9201593348566 951.4201593348566
            ",
        );
    }
}
