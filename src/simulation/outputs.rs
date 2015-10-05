/* Cymbalum, Molecular Simulation in Rust - Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
 */

//! Saving properties of the simulation to some stream (file or stdout mainly)

extern crate chemharp;

use std::io::prelude::*;
use std::io;
use std::fs::File;
use std::path::Path;

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
/// file, using any format supported by the
/// [Chemharp](http://chemharp.readthedocs.org/en/latest/formats.html) library.
pub struct TrajectoryOutput {
    file: chemharp::Trajectory,
}

impl TrajectoryOutput {
    /// Create a new `TrajectoryOutput` writing to `filename`. The file is
    /// replaced if it already exists.
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
/// The `CellOutput` write all the components of a cell to a file . The columns
/// in the file contains the following values: `A B C α β γ`.
pub struct CellOutput {
    file: File,
    path: String
}

impl CellOutput {
    /// Create a new `CellOutput` writing to `filename`. The file is replaced if
    /// it already exists.
    pub fn new<P: AsRef<Path>>(filename: P) -> Result<CellOutput, io::Error> {
        let path = String::from(filename.as_ref().to_str().unwrap());
        Ok(CellOutput{
            file: try!(File::create(filename)),
            path: path,
        })
    }
}

impl Output for CellOutput {
    fn setup(&mut self, _: &Universe) {
        if let Err(e) = writeln!(&mut self.file, "# Unit cell of the simulation") {
            // Do panic in early time
            panic!("Could not write to file '{}': {:?}", self.path, e);
        }
        if let Err(e) = writeln!(&mut self.file, "# Step A/Å B/Å C/Å α/deg β/deg γ/deg") {
            panic!("Could not write to file '{}': {:?}", self.path, e);
        }
    }

    fn write(&mut self, universe: &Universe) {
        let cell = universe.cell();
        if let Err(e) = writeln!(
            &mut self.file, "{} {} {} {} {} {} {}",
            universe.step(), cell.a(), cell.b(), cell.c(), cell.alpha(), cell.beta(), cell.gamma()
        ) {
            // Do not panic during the simulation
            error!("Could not write to file '{}': {:?}", self.path, e);
        }
    }
}

/******************************************************************************/
/// The `EnergyOutput` write the energy of the system to a text file, organized
/// as: `PotentialEnergy     KineticEnergy     TotalEnergy`.
pub struct EnergyOutput {
    file: File,
    path: String
}

impl EnergyOutput {
    /// Create a new `EnergyOutput` writing to `filename`. The file is replaced
    /// if it already exists.
    pub fn new<P: AsRef<Path>>(filename: P) -> Result<EnergyOutput, io::Error> {
        let path = String::from(filename.as_ref().to_str().unwrap());
        Ok(EnergyOutput{
            file: try!(File::create(filename)),
            path: path,
        })
    }
}

impl Output for EnergyOutput {
    fn setup(&mut self, _: &Universe) {
        if let Err(e) = writeln!(&mut self.file, "# Energy of the simulation (kJ/mol)") {
            panic!("Could not write to file '{}': {:?}", self.path, e);
        }
        if let Err(e) = writeln!(&mut self.file, "# Step Potential Kinetic Total") {
            panic!("Could not write to file '{}': {:?}", self.path, e);
        }
    }

    fn write(&mut self, universe: &Universe) {
        let potential = units::to(universe.potential_energy(), "kJ/mol").unwrap();
        let kinetic = units::to(universe.kinetic_energy(), "kJ/mol").unwrap();
        let total = units::to(universe.total_energy(), "kJ/mol").unwrap();
        if let Err(e) = writeln!(&mut self.file, "{} {} {} {}", universe.step(), potential, kinetic, total) {
            error!("Could not write to file '{}': {:?}", self.path, e);
        }
    }
}

/******************************************************************************/
/// The `PropertiesOutput` write various physical properties of the universe to
/// a file. These properties are:
///
/// - volume of the unit cell;
/// - instant temperature;
/// - instant pressure;
pub struct PropertiesOutput {
    file: File,
    path: String,
}

impl PropertiesOutput {
    /// Create a new `PropertiesOutput` writing to `filename`. The file is replaced
    /// if it already exists.
    pub fn new<P: AsRef<Path>>(filename: P) -> Result<PropertiesOutput, io::Error> {
        let path = String::from(filename.as_ref().to_str().unwrap());
        Ok(PropertiesOutput{
            file: try!(File::create(filename)),
            path: path,
        })
    }
}

impl Output for PropertiesOutput {
    fn setup(&mut self, _: &Universe) {
        if let Err(e) = writeln!(&mut self.file, "# Physical properties of the simulation") {
            panic!("Could not write to file '{}': {:?}", self.path, e);
        }
        if let Err(e) = writeln!(&mut self.file, "# Step Volume/A^3 Temperature/K Pressure/bar") {
            panic!("Could not write to file '{}': {:?}", self.path, e);
        }
    }

    fn write(&mut self, universe: &Universe) {
        let V = units::to(universe.volume(), "A^3").unwrap();
        let T = units::to(universe.temperature(), "K").unwrap();
        let P = units::to(universe.pressure(), "bar").unwrap();
        if let Err(e) = writeln!(&mut self.file, "{} {} {} {}", universe.step(), V, T, P) {
            error!("Could not write to file '{}': {:?}", self.path, e);
        }
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
            Box::new(Harmonic{k: units::from(300.0, "kJ/mol/A^2").unwrap(), x0: units::from(1.2, "A").unwrap()}));
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
            let mut out = EnergyOutput::new(filename).unwrap();
            out.setup(&universe);
            out.write(&universe);
            out.finish(&universe);
        }

        check_file_content(filename, "# Energy of the simulation (kJ/mol)\n# Step Potential Kinetic Total\n0 1.5000000000000027 0 1.5000000000000027\n");
        fs::remove_file(filename).unwrap();
    }

    #[test]
    fn cell() {
        let filename = "testing-cell-output.dat";
        let universe = testing_universe();
        {
            let mut out = CellOutput::new(filename).unwrap();
            out.setup(&universe);
            out.write(&universe);
            out.finish(&universe);
        }

        check_file_content(filename, "# Unit cell of the simulation\n# Step A/Å B/Å C/Å α/deg β/deg γ/deg\n0 10 10 10 90 90 90\n");
        fs::remove_file(filename).unwrap();
    }

    #[test]
    fn properties() {
        let filename = "testing-properties-output.dat";
        let universe = testing_universe();
        {
            let mut out = PropertiesOutput::new(filename).unwrap();
            out.setup(&universe);
            out.write(&universe);
            out.finish(&universe);
        }

        check_file_content(filename, "# Physical properties of the simulation\n# Step Volume/A^3 Temperature/K Pressure/bar\n0 1000 0 -431.7400836223091\n");
        fs::remove_file(filename).unwrap();
    }
}
