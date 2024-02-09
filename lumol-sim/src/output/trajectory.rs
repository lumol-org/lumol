// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

use std::path::Path;

use super::Output;

use lumol_core::{OpenMode, Trajectory, TrajectoryBuilder, TrajectoryError};
use lumol_core::System;

/// The `TrajectoryOutput` allows to write the trajectory of the system to a
/// file, using any format supported by the [Chemfiles][chemfiles] library.
///
/// [chemfiles]: http://chemfiles.github.io
pub struct TrajectoryOutput {
    file: Trajectory,
}

impl TrajectoryOutput {
    /// Create a new `TrajectoryOutput` writing to `filename`. The file is
    /// replaced if it already exists. The file format is guessed from the
    /// extension. Please refer to the list of [supported formats][formats]
    /// for more information.
    ///
    /// [formats]: http://chemfiles.org/chemfiles/latest/formats.html
    pub fn new<P>(path: P) -> Result<TrajectoryOutput, TrajectoryError>
    where
        P: AsRef<Path>,
    {
        let builder = TrajectoryBuilder::new().mode(OpenMode::Write);
        Ok(TrajectoryOutput {
            file: builder.open(path)?
        })
    }

    /// Create a new `TrajectoryOutput` writing to `filename` using the given
    /// `format`. The file is replaced if it already exists.
    ///
    /// Please refer to the list of [supported formats][formats] for more
    /// information.
    ///
    /// [formats]: http://chemfiles.org/chemfiles/latest/formats.html
    pub fn with_format<P>(path: P, format: &str) -> Result<TrajectoryOutput, TrajectoryError>
    where
        P: AsRef<Path>,
    {
        let builder = TrajectoryBuilder::new().mode(OpenMode::Write).format(format);
        Ok(TrajectoryOutput {
            file: builder.open(path)?
        })
    }
}

impl Output for TrajectoryOutput {
    fn write(&mut self, system: &System) {
        match self.file.write(system) {
            Ok(()) => (),
            Err(err) => {
                panic!("Error in while writing trajectory: {}", err);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::tests::test_output;

    #[test]
    fn cell() {
        test_output(
            |path| Box::new(TrajectoryOutput::with_format(path, "XYZ").unwrap()),
            "2
            Properties=species:S:1:pos:R:3 Lattice=\"10 0 0 0 10 0 0 0 10\"
            F 0 0 0
            F 1.3 0 0
            ",
        );
    }
}
