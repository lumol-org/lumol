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

/// The `ForcesOutput` writes the forces acting on the atoms using XYZ format
pub struct ForcesOutput {
    file: BufWriter<File>,
    path: PathBuf,
}

impl ForcesOutput {
    /// Create a new `ForcesOutput` writing to `filename`. The file is replaced
    /// if it already exists.
    pub fn new<P: AsRef<Path>>(filename: P) -> Result<ForcesOutput, io::Error> {
        Ok(ForcesOutput {
            file: BufWriter::new(File::create(filename.as_ref())?),
            path: filename.as_ref().to_owned(),
        })
    }
}

impl Output for ForcesOutput {
    fn setup(&mut self, _: &System) {}

    fn write(&mut self, system: &System) {
        let forces = system.forces();
        let names = system.particles().name;
        let conversion = units::to(1.0, "kJ/mol/A").expect("bad unit");

        writeln_or_log!(self, "{}", forces.len());
        writeln_or_log!(self, "forces in kJ/mol/A at step {}", system.step);
        for (i, force) in forces.iter().enumerate() {
            let x = conversion * force[0];
            let y = conversion * force[1];
            let z = conversion * force[2];
            writeln_or_log!(self, "{} {} {} {}", names[i], x, y, z);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::tests::test_output;

    #[test]
    fn energy() {
        test_output(
            |path| Box::new(ForcesOutput::new(path).unwrap()),
            "2
            forces in kJ/mol/A at step 42
            F 30.000000000000025 0 0
            F -30.000000000000025 0 0
            ",
        );
    }
}
