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

/// The `PropertiesOutput` write various physical properties of the system to
/// a file. These properties are:
///
/// - volume of the unit cell;
/// - instant temperature;
/// - instant pressure;
pub struct PropertiesOutput {
    file: BufWriter<File>,
    path: PathBuf,
}

impl PropertiesOutput {
    /// Create a new `PropertiesOutput` writing to `filename`. The file is replaced
    /// if it already exists.
    pub fn new<P: AsRef<Path>>(filename: P) -> Result<PropertiesOutput, io::Error> {
        Ok(PropertiesOutput {
            file: BufWriter::new(File::create(filename.as_ref())?),
            path: filename.as_ref().to_owned(),
        })
    }
}

impl Output for PropertiesOutput {
    fn setup(&mut self, _: &System) {
        writeln_or_log!(self, "# Physical properties of the simulation");
        writeln_or_log!(self, "# Step Volume/A^3 Temperature/K Pressure/bar");
    }

    fn write(&mut self, system: &System) {
        let volume = units::to(system.volume(), "A^3").expect("bad unit");
        let temperature = units::to(system.temperature(), "K").expect("bad unit");
        let pressure = units::to(system.pressure(), "bar").expect("bad unit");
        writeln_or_log!(self, "{} {} {} {}", system.step, volume, temperature, pressure);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::tests::test_output;

    #[test]
    fn properties() {
        test_output(
            |path| Box::new(PropertiesOutput::new(path).unwrap()),
            "# Physical properties of the simulation
            # Step Volume/A^3 Temperature/K Pressure/bar
            42 1000 38083.04389172312 10299.991728079816
            ",
        );
    }
}
