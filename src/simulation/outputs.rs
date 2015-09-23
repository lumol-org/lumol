/*
 * Cymbalum, Molecular Simulation in Rust
 * Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
*/

//! Saving properties of the simulation to some stream (file or stdout mainly)

extern crate chemharp;

use std::io::prelude::*;
use std::fs::File;

use ::units;
use ::universe::Universe;

use ::universe::chemharp::universe_to_frame;

/// The `Output` trait define the interface for all the quantities outputed by
/// the simulation during the run. An Output can be a text or a binary data
/// file, an image, a text log, …
pub trait Output {
    /// Function called once at the beggining of the simulation, which allow
    /// for some setup of the output if needed.
    fn setup(&mut self, _: &Universe) {}

    /// Write the output from the universe.
    fn write(&mut self, universe: &Universe);

    /// Function called once at the end of the simulation.
    fn finish(&mut self, _: &Universe) {}
}

/******************************************************************************/
/// The `TrajectoryOutput` allow to write the trajectory of the system to a
/// file, using the XYZ format.
pub struct TrajectoryOutput {
    file: chemharp::Trajectory,
}

impl TrajectoryOutput {
    /// Create a new `TrajectoryOutput`, using the file ate `filename` for
    /// as the trajectory file. The file is replaced if it already exists.
    pub fn new<'a, S>(filename: S) -> Result<TrajectoryOutput, chemharp::Error> where S: Into<&'a str> {
        Ok(TrajectoryOutput{
            file: try!(chemharp::Trajectory::create(filename.into()))
        })
    }
}

impl Output for TrajectoryOutput {
    fn write(&mut self, universe: &Universe) {
        let frame = match universe_to_frame(universe) {
            Ok(val) => val,
            Err(err) => {
                error!("Error in Chemharp runtime while converting data: {}", err.message());
                panic!();
            }
        };

        match self.file.write(&frame) {
            Ok(()) => {},
            Err(err) => {
                error!("Error in Chemharp runtime while writing file: {}", err.message());
                panic!();
            }
        }
    }
}

/******************************************************************************/
/// The `EnergyOutput` write the energy of the system to a text file, organized
/// as: `PotentialEnergy     KineticEnergy     TotalEnergy`.
pub struct EnergyOutput {
    file: File,
}

impl EnergyOutput {
    /// Create a new `TrajectoryOutput`, using the file ate `filename` for
    /// output.
    pub fn new<'a, S>(filename: S) -> EnergyOutput where S: Into<&'a str> {
        EnergyOutput{
            file: File::create(filename.into()).unwrap()
        }
    }
}

impl Output for EnergyOutput {
    fn setup(&mut self, _: &Universe) {
        writeln!(&mut self.file, "# Energy of the simulation (kJ/mol)").unwrap();
        writeln!(&mut self.file, "# Potential     Kinetic     Total").unwrap();
    }

    fn write(&mut self, universe: &Universe) {
        let potential = units::to(universe.potential_energy(), "kJ/mol").unwrap();
        let kinetic = units::to(universe.kinetic_energy(), "kJ/mol").unwrap();
        let total = units::to(universe.total_energy(), "kJ/mol").unwrap();
        writeln!(&mut self.file, "{}   {}   {}", potential, kinetic, total).unwrap();
    }
}

#[cfg(test)]
mod tests {
    use std::io::prelude::*;
    use std::fs;

    use super::*;
    use ::universe::*;
    use ::types::*;
    use ::potentials::*;
    use ::units;

    fn testing_universe() -> Universe {
        let mut universe = Universe::from_cell(UnitCell::cubic(10.0));;

        universe.add_particle(Particle::new("F"));
        universe[0].set_position(Vector3D::new(0.0, 0.0, 0.0));

        universe.add_particle(Particle::new("F"));
        universe[1].set_position(Vector3D::new(1.3, 0.0, 0.0));

        universe.add_pair_interaction("F", "F",
            Harmonic{k: units::from(300.0, "kJ/mol/A^2").unwrap(), x0: units::from(1.2, "A").unwrap()});
        return universe;
    }

    fn check_file_content(filename: &str, content: &str) {
        let mut file = fs::File::open(filename).unwrap();
        let mut buffer = String::new();
        file.read_to_string(&mut buffer).unwrap();

        assert_eq!(buffer, content);
    }

    #[test]
    fn trajectory() {
        let filename = "testing-trajectory-output.xyz";
        let universe = testing_universe();
        {
            let mut out = TrajectoryOutput::new(filename).unwrap();
            out.setup(&universe);
            out.write(&universe);
            out.finish(&universe);
        }

        check_file_content(filename, "2\nWritten by Chemharp\nF 0 0 0\nF 1.3 0 0\n");
        fs::remove_file(filename).unwrap();
    }

    #[test]
    fn energy() {
        let filename = "testing-energy-output.dat";
        let universe = testing_universe();
        {
            let mut out = EnergyOutput::new(filename);
            out.setup(&universe);
            out.write(&universe);
            out.finish(&universe);
        }

        check_file_content(filename, "# Energy of the simulation (kJ/mol)\n# Potential     Kinetic     Total\n1.5000000000000027   0   1.5000000000000027\n");
        fs::remove_file(filename).unwrap();
    }
}
