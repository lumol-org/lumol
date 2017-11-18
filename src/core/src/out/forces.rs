// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

use std::fs::File;
use std::io;
use std::io::prelude::*;
use std::path::{Path, PathBuf};

use super::Output;
use sys::System;
use utils;

/// The `ForcesOutput` writes the forces acting on the atoms using XYZ format
pub struct ForcesOutput {
    file: File,
    path: PathBuf,
}

impl ForcesOutput {
    /// Create a new `ForcesOutput` writing to `filename`. The file is replaced
    /// if it already exists.
    pub fn new<P: AsRef<Path>>(filename: P) -> Result<ForcesOutput, io::Error> {
        Ok(ForcesOutput {
            file: try!(File::create(filename.as_ref())),
            path: filename.as_ref().to_owned(),
        })
    }
}

macro_rules! try {
    ($e: expr, $path: expr) => (
        if let Err(err) = $e {
            error!("Could not write to file '{}': {}", $path.display(), err);
            return;
        }
    );
}

impl Output for ForcesOutput {
    fn setup(&mut self, _: &System) {}

    fn write(&mut self, system: &System) {
        let forces = system.forces();
        let names = system.particles().name;
        let conversion = utils::unit_to(1.0, "kJ/mol/A");

        try!(writeln!(&mut self.file, "{}", forces.len()), self.path);
        try!(
            writeln!(&mut self.file, "forces in kJ/mol/A at step {}", system.step()),
            self.path
        );
        for (i, force) in forces.iter().enumerate() {
            let x = conversion * force[0];
            let y = conversion * force[1];
            let z = conversion * force[2];
            try!(writeln!(&mut self.file, "{} {} {} {}", names[i], x, y, z), self.path);
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
            forces in kJ/mol/A at step 0
            F 30.000000000000025 0 0
            F -30.000000000000025 0 0
            ",
        );
    }
}
