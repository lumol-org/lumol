// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

use std::io::prelude::*;
use std::io;
use std::fs::File;
use std::path::{Path, PathBuf};

use sys::System;
use super::Output;

/// The `CellOutput` writes all the components of a cell to a file . The columns
/// in the file contain the following values: `step A B C α β γ`.
pub struct CellOutput {
    file: File,
    path: PathBuf
}

impl CellOutput {
    /// Create a new `CellOutput` writing to `filename`. The file is replaced if
    /// it already exists.
    pub fn new<P: AsRef<Path>>(filename: P) -> Result<CellOutput, io::Error> {
        Ok(CellOutput{
            file: try!(File::create(filename.as_ref())),
            path: filename.as_ref().to_owned(),
        })
    }
}

impl Output for CellOutput {
    fn setup(&mut self, _: &System) {
        if let Err(err) = writeln!(&mut self.file, "# Unit cell of the simulation") {
            // Do panic in early time
            fatal_error!("Could not write to file '{}': {}", self.path.display(), err);
        }
        if let Err(err) = writeln!(&mut self.file, "# Step A/Å B/Å C/Å α/deg β/deg γ/deg") {
            fatal_error!("Could not write to file '{}': {}", self.path.display(), err);
        }
    }

    fn write(&mut self, system: &System) {
        let cell = &system.cell;
        if let Err(err) = writeln!(
            &mut self.file, "{} {} {} {} {} {} {}",
            system.step(), cell.a(), cell.b(), cell.c(), cell.alpha(), cell.beta(), cell.gamma()
        ) {
            // Do not panic during the simulation
            error!("Could not write to file '{}': {}", self.path.display(), err);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::tests::test_output;

    #[test]
    fn cell() {
        test_output(|path| {
            Box::new(CellOutput::new(path).unwrap())
        },
"# Unit cell of the simulation
# Step A/Å B/Å C/Å α/deg β/deg γ/deg
0 10 10 10 90 90 90
"
        );
    }
}
