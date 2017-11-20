// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

use std::fs::File;
use std::io;
use std::io::prelude::*;
use std::path::{Path, PathBuf};

use super::Output;

use sys::System;
use utils;

/// The `PropertiesOutput` write various physical properties of the system to
/// a file. These properties are:
///
/// - volume of the unit cell;
/// - instant temperature;
/// - instant pressure;
pub struct PropertiesOutput {
    file: File,
    path: PathBuf,
}

impl PropertiesOutput {
    /// Create a new `PropertiesOutput` writing to `filename`. The file is replaced
    /// if it already exists.
    pub fn new<P: AsRef<Path>>(filename: P) -> Result<PropertiesOutput, io::Error> {
        Ok(PropertiesOutput {
            file: try!(File::create(filename.as_ref())),
            path: filename.as_ref().to_owned(),
        })
    }
}

impl Output for PropertiesOutput {
    fn setup(&mut self, _: &System) {
        if let Err(err) = writeln!(&mut self.file, "# Physical properties of the simulation") {
            fatal_error!("Could not write to file '{}': {}", self.path.display(), err);
        }
        if let Err(err) = writeln!(&mut self.file, "# Step Volume/A^3 Temperature/K Pressure/bar") {
            fatal_error!("Could not write to file '{}': {}", self.path.display(), err);
        }
    }

    fn write(&mut self, system: &System) {
        let volume = utils::unit_to(system.volume(), "A^3");
        let temperature = utils::unit_to(system.temperature(), "K");
        let pressure = utils::unit_to(system.pressure(), "bar");
        if let Err(err) =
            writeln!(&mut self.file, "{} {} {} {}", system.step(), volume, temperature, pressure)
        {
            error!("Could not write to file '{}': {}", self.path.display(), err);
        }
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
