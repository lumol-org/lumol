// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

use std::fs::File;
use std::io::{self, BufWriter};
use std::io::prelude::*;
use std::path::{Path, PathBuf};

use log::error;

use super::Output;
use lumol_core::System;

/// The `CellOutput` writes all the components of a cell to a file . The columns
/// in the file contain the following values: `step A B C α β γ`.
pub struct CellOutput {
    file: BufWriter<File>,
    path: PathBuf,
}

impl CellOutput {
    /// Create a new `CellOutput` writing to `filename`. The file is replaced if
    /// it already exists.
    pub fn new<P: AsRef<Path>>(filename: P) -> Result<CellOutput, io::Error> {
        Ok(CellOutput {
            file: BufWriter::new(File::create(filename.as_ref())?),
            path: filename.as_ref().to_owned(),
        })
    }
}

impl Output for CellOutput {
    #[allow(clippy::non_ascii_literal)]
    fn setup(&mut self, _: &System) {
        writeln_or_log!(self, "# Unit cell of the simulation");
        writeln_or_log!(self, "# Step A/Å B/Å C/Å α/deg β/deg γ/deg");
    }

    fn write(&mut self, system: &System) {
        writeln_or_log!(self, "{} {} {} {} {} {} {}",
            system.step,
            system.cell.a(),
            system.cell.b(),
            system.cell.c(),
            system.cell.alpha(),
            system.cell.beta(),
            system.cell.gamma()
        );
    }
}

#[cfg(test)]
#[allow(clippy::non_ascii_literal)]
mod tests {
    use super::*;
    use super::super::tests::test_output;

    #[test]
    fn cell() {
        test_output(
            |path| Box::new(CellOutput::new(path).unwrap()),
            "# Unit cell of the simulation
            # Step A/Å B/Å C/Å α/deg β/deg γ/deg
            42 10 10 10 90 90 90
            ",
        );
    }
}
