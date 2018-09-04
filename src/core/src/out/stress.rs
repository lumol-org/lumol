// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

use std::fs::File;
use std::io;
use std::io::prelude::*;
use std::path::{Path, PathBuf};

use super::Output;
use sys::System;
use utils;

/// The `StressOutput` writes the stress of the system to a text file, organized
/// as: `step stress.xx stress.yy stress.zz stress.xy stress.xz stress.yz`.
pub struct StressOutput {
    file: File,
    path: PathBuf,
}

impl StressOutput {
    /// Create a new `StressOutput` writing to `filename`. The file is replaced
    /// if it already exists.
    pub fn new<P: AsRef<Path>>(filename: P) -> Result<StressOutput, io::Error> {
        Ok(StressOutput {
            file: File::create(filename.as_ref())?,
            path: filename.as_ref().to_owned(),
        })
    }
}

impl Output for StressOutput {
    fn setup(&mut self, _: &System) {
        if let Err(err) = writeln!(&mut self.file, "# Stress tensor of the simulation (bar)") {
            fatal_error!("Could not write to file '{}': {}", self.path.display(), err);
        }
        if let Err(err) = writeln!(
            &mut self.file,
            "# step stress.xx stress.yy stress.zz stress.xy stress.xz stress.yz"
        ) {
            fatal_error!("Could not write to file '{}': {}", self.path.display(), err);
        }
    }

    fn write(&mut self, system: &System) {
        let stress = system.stress();
        let xx = utils::unit_to(stress[0][0], "bar");
        let yy = utils::unit_to(stress[1][1], "bar");
        let zz = utils::unit_to(stress[2][2], "bar");
        let xy = utils::unit_to(stress[0][1], "bar");
        let xz = utils::unit_to(stress[0][2], "bar");
        let yz = utils::unit_to(stress[1][2], "bar");
        writeln_or_log!(self, "{} {} {} {} {} {} {}", system.step(), xx, yy, zz, xy, xz, yz);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::tests::test_output;

    #[test]
    fn energy() {
        test_output(
            |path| Box::new(StressOutput::new(path).unwrap()),
            "# Stress tensor of the simulation (bar)
            # step stress.xx stress.yy stress.zz stress.xy stress.xz stress.yz
            42 30899.975184239443 0 0 0 0 0
            ",
        );
    }
}
