// Cymbalum, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

//! Saving properties of the simulation to some stream (file or stdout mainly)
use std::io::prelude::*;
use std::io;
use std::fs::File;
use std::path::{Path, PathBuf};

use utils;
use system::System;
use input::{Trajectory, TrajectoryError};

/// The `Output` trait define the interface for all the quantities outputed by
/// the simulation during the run. An Output can be a text or a binary data
/// file, an image, a text log, …
pub trait Output {
    /// Function called once at the beggining of the simulation, which allow
    /// for some setup of the output if needed.
    fn setup(&mut self, _: &System) {}

    /// Write the output from the system.
    fn write(&mut self, system: &System);

    /// Function called once at the end of the simulation.
    fn finish(&mut self, _: &System) {}
}

/******************************************************************************/
/// The `TrajectoryOutput` allow to write the trajectory of the system to a
/// file, using any format supported by the
/// [Chemfiles](http://chemfiles.readthedocs.org/en/latest/formats.html) library.
pub struct TrajectoryOutput {
    file: Trajectory,
}

impl TrajectoryOutput {
    /// Create a new `TrajectoryOutput` writing to `filename`. The file is
    /// replaced if it already exists.
    pub fn new<P>(path: P) -> Result<TrajectoryOutput, TrajectoryError> where P: AsRef<Path> {
        Ok(TrajectoryOutput{
            file: try!(Trajectory::create(path))
        })
    }
}

impl Output for TrajectoryOutput {
    fn write(&mut self, system: &System) {
        match self.file.write(system) {
            Ok(()) => (),
            Err(err) => {
                fatal_error!("Error in while writing trajectory: {}", err);
            }
        }
    }
}


/******************************************************************************/
/// The `CellOutput` write all the components of a cell to a file . The columns
/// in the file contains the following values: `A B C α β γ`.
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
            fatal_error!("Could not write to file '{}': {:?}", self.path.display(), err);
        }
        if let Err(err) = writeln!(&mut self.file, "# Step A/Å B/Å C/Å α/deg β/deg γ/deg") {
            fatal_error!("Could not write to file '{}': {:?}", self.path.display(), err);
        }
    }

    fn write(&mut self, system: &System) {
        let cell = system.cell();
        if let Err(err) = writeln!(
            &mut self.file, "{} {} {} {} {} {} {}",
            system.step(), cell.a(), cell.b(), cell.c(), cell.alpha(), cell.beta(), cell.gamma()
        ) {
            // Do not panic during the simulation
            error!("Could not write to file '{}': {:?}", self.path.display(), err);
        }
    }
}

/******************************************************************************/
/// The `EnergyOutput` write the energy of the system to a text file, organized
/// as: `PotentialEnergy     KineticEnergy     TotalEnergy`.
pub struct EnergyOutput {
    file: File,
    path: PathBuf
}

impl EnergyOutput {
    /// Create a new `EnergyOutput` writing to `filename`. The file is replaced
    /// if it already exists.
    pub fn new<P: AsRef<Path>>(filename: P) -> Result<EnergyOutput, io::Error> {
        Ok(EnergyOutput{
            file: try!(File::create(filename.as_ref())),
            path: filename.as_ref().to_owned(),
        })
    }
}

impl Output for EnergyOutput {
    fn setup(&mut self, _: &System) {
        if let Err(err) = writeln!(&mut self.file, "# Energy of the simulation (kJ/mol)") {
            fatal_error!("Could not write to file '{}': {:?}", self.path.display(), err);
        }
        if let Err(err) = writeln!(&mut self.file, "# Step Potential Kinetic Total") {
            fatal_error!("Could not write to file '{}': {:?}", self.path.display(), err);
        }
    }

    fn write(&mut self, system: &System) {
        let potential = utils::unit_to(system.potential_energy(), "kJ/mol");
        let kinetic = utils::unit_to(system.kinetic_energy(), "kJ/mol");
        let total = utils::unit_to(system.total_energy(), "kJ/mol");
        if let Err(err) = writeln!(&mut self.file, "{} {} {} {}", system.step(), potential, kinetic, total) {
            error!("Could not write to file '{}': {:?}", self.path.display(), err);
        }
    }
}

/******************************************************************************/
/// The `PropertiesOutput` write various physical properties of the system to
/// a file. These properties are:
///
/// - volume of the unit cell;
/// - instant temperature;
/// - instant pressure;
pub struct PropertiesOutput {
    file: File,
    path: PathBuf,
}

impl PropertiesOutput {
    /// Create a new `PropertiesOutput` writing to `filename`. The file is replaced
    /// if it already exists.
    pub fn new<P: AsRef<Path>>(filename: P) -> Result<PropertiesOutput, io::Error> {
        Ok(PropertiesOutput{
            file: try!(File::create(filename.as_ref())),
            path: filename.as_ref().to_owned(),
        })
    }
}

impl Output for PropertiesOutput {
    fn setup(&mut self, _: &System) {
        if let Err(err) = writeln!(&mut self.file, "# Physical properties of the simulation") {
            fatal_error!("Could not write to file '{}': {:?}", self.path.display(), err);
        }
        if let Err(err) = writeln!(&mut self.file, "# Step Volume/A^3 Temperature/K Pressure/bar") {
            fatal_error!("Could not write to file '{}': {:?}", self.path.display(), err);
        }
    }

    fn write(&mut self, system: &System) {
        let volume = utils::unit_to(system.volume(), "A^3");
        let temperature = utils::unit_to(system.temperature(), "K");
        let pressure = utils::unit_to(system.pressure(), "bar");
        if let Err(err) = writeln!(&mut self.file, "{} {} {} {}", system.step(), volume, temperature, pressure) {
            error!("Could not write to file '{}': {:?}", self.path.display(), err);
        }
    }
}

#[cfg(test)]
mod tests {
    extern crate tempfile;
    use self::tempfile::{NamedTempFile, NamedTempFileOptions};

    use std::io::prelude::*;
    use std::fs::File;

    use super::*;
    use system::*;
    use types::*;
    use potentials::*;
    use utils::unit_from;

    fn testing_system() -> System {
        let mut system = System::from_cell(UnitCell::cubic(10.0));;

        system.add_particle(Particle::new("F"));
        system[0].position = Vector3D::zero();

        system.add_particle(Particle::new("F"));
        system[1].position = Vector3D::new(1.3, 0.0, 0.0);

        system.add_pair_interaction("F", "F",
            Box::new(Harmonic{k: unit_from(300.0, "kJ/mol/A^2"), x0: unit_from(1.2, "A")}));
        return system;
    }

    fn check_file_content(mut file: File, content: &str) {
        let mut buffer = String::new();
        let _ = file.read_to_string(&mut buffer).unwrap();

        assert_eq!(buffer, content);
    }

    #[test]
    fn trajectory() {
        let tempfile = NamedTempFileOptions::new().suffix(".xyz").create().unwrap();
        let system = testing_system();
        {
            let mut out = TrajectoryOutput::new(tempfile.path()).unwrap();
            out.setup(&system);
            out.write(&system);
            out.finish(&system);
        }

        let content = "2
Written by the chemfiles library
F 0 0 0
F 1.3 0 0
";
        let file = tempfile.reopen().unwrap();
        check_file_content(file, content);
    }

    #[test]
    fn energy() {
        let tempfile = NamedTempFile::new().unwrap();
        let system = testing_system();
        {
            let mut out = EnergyOutput::new(tempfile.path()).unwrap();
            out.setup(&system);
            out.write(&system);
            out.finish(&system);
        }

        let content = "\
# Energy of the simulation (kJ/mol)
# Step Potential Kinetic Total
0 1.5000000000000027 0 1.5000000000000027
";

        let file = tempfile.reopen().unwrap();
        check_file_content(file, content);
    }

    #[test]
    fn cell() {
        let tempfile = NamedTempFile::new().unwrap();
        let system = testing_system();
        {
            let mut out = CellOutput::new(tempfile.path()).unwrap();
            out.setup(&system);
            out.write(&system);
            out.finish(&system);
        }

        let content = "\
# Unit cell of the simulation
# Step A/Å B/Å C/Å α/deg β/deg γ/deg
0 10 10 10 90 90 90
";

        let file = tempfile.reopen().unwrap();
        check_file_content(file, content);
    }

    #[test]
    fn properties() {
        let tempfile = NamedTempFile::new().unwrap();
        let system = testing_system();
        {
            let mut out = PropertiesOutput::new(tempfile.path()).unwrap();
            out.setup(&system);
            out.write(&system);
            out.finish(&system);
        }

        let content = "\
# Physical properties of the simulation
# Step Volume/A^3 Temperature/K Pressure/bar
0 1000 0 -215.87004181115455
";

        let file = tempfile.reopen().unwrap();
        check_file_content(file, content);
    }
}
