// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

use std::io::prelude::*;
use std::io;
use std::fs::File;
use std::path::{Path, PathBuf};

use super::Output;
use sys::System;

/// The `ForcesOutput` writes the forces acting on the atoms using XYZ format
pub struct ForcesOutput {
    file: File,
    path: PathBuf
}

impl ForcesOutput {
    /// Create a new `ForcesOutput` writing to `filename`. The file is replaced
    /// if it already exists.
    pub fn new<P: AsRef<Path>>(filename: P) -> Result<ForcesOutput, io::Error> {
        Ok(ForcesOutput{
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

        try!(writeln!(&mut self.file, "{}", forces.len()), self.path);
        try!(writeln!(&mut self.file, "forces at step {}", system.step()), self.path);
        for (i, force) in forces.iter().enumerate() {
            try!(writeln!(&mut self.file, "{} {} {} {}", names[i], force[0], force[1], force[2]), self.path);
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
            Box::new(ForcesOutput::new(path).unwrap())
        },
"2
forces at step 0
F 0.0030000000021006322 0 0
F -0.0030000000021006322 0 0
"
        );
    }
}
