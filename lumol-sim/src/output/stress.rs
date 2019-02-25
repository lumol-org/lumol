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

/// The `StressOutput` writes the stress of the system to a text file, organized
/// as: `step stress.xx stress.yy stress.zz stress.xy stress.xz stress.yz`.
pub struct StressOutput {
    file: BufWriter<File>,
    path: PathBuf,
}

impl StressOutput {
    /// Create a new `StressOutput` writing to `filename`. The file is replaced
    /// if it already exists.
    pub fn new<P: AsRef<Path>>(filename: P) -> Result<StressOutput, io::Error> {
        Ok(StressOutput {
            file: BufWriter::new(File::create(filename.as_ref())?),
            path: filename.as_ref().to_owned(),
        })
    }
}

impl Output for StressOutput {
    fn setup(&mut self, _: &System) {
        if let Err(err) = writeln!(&mut self.file, "# Stress tensor of the simulation (bar)") {
            panic!("Could not write to file '{}': {}", self.path.display(), err);
        }
        if let Err(err) = writeln!(
            &mut self.file,
            "# step stress.xx stress.yy stress.zz stress.xy stress.xz stress.yz"
        ) {
            panic!("Could not write to file '{}': {}", self.path.display(), err);
        }
    }

    fn write(&mut self, system: &System) {
        let conversion = units::to(1.0, "bar").expect("bad unit");
        let stress = system.stress();
        let xx = stress[0][0] * conversion;
        let yy = stress[1][1] * conversion;
        let zz = stress[2][2] * conversion;
        let xy = stress[0][1] * conversion;
        let xz = stress[0][2] * conversion;
        let yz = stress[1][2] * conversion;
        writeln_or_log!(self, "{} {} {} {} {} {} {}", system.step, xx, yy, zz, xy, xz, yz);
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
