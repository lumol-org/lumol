// Cymbalum, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

//! Saving properties of the simulation to some stream (file or stdout mainly)
use std::io::prelude::*;
use std::io;
use std::fs::File;
use std::path::Path;

use units;
use system::System;
use io::{Trajectory, TrajectoryResult};

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
    pub fn new<'a, S>(filename: S) -> TrajectoryResult<TrajectoryOutput> where S: Into<&'a str> {
        Ok(TrajectoryOutput{
            file: try!(Trajectory::create(filename.into()))
        })
    }
}

impl Output for TrajectoryOutput {
    fn write(&mut self, system: &System) {
        match self.file.write(system) {
            Ok(()) => (),
            Err(err) => {
                error!("Error in while writing trajectory: {}", err.message());
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
    fn setup(&mut self, _: &System) {
        if let Err(e) = writeln!(&mut self.file, "# Unit cell of the simulation") {
            // Do panic in early time
            panic!("Could not write to file '{}': {:?}", self.path, e);
        }
        if let Err(e) = writeln!(&mut self.file, "# Step A/Å B/Å C/Å α/deg β/deg γ/deg") {
            panic!("Could not write to file '{}': {:?}", self.path, e);
        }
    }

    fn write(&mut self, system: &System) {
        let cell = system.cell();
        if let Err(e) = writeln!(
            &mut self.file, "{} {} {} {} {} {} {}",
            system.step(), cell.a(), cell.b(), cell.c(), cell.alpha(), cell.beta(), cell.gamma()
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
    fn setup(&mut self, _: &System) {
        if let Err(e) = writeln!(&mut self.file, "# Energy of the simulation (kJ/mol)") {
            panic!("Could not write to file '{}': {:?}", self.path, e);
        }
        if let Err(e) = writeln!(&mut self.file, "# Step Potential Kinetic Total") {
            panic!("Could not write to file '{}': {:?}", self.path, e);
        }
    }

    fn write(&mut self, system: &System) {
        let potential = units::to(system.potential_energy(), "kJ/mol").unwrap();
        let kinetic = units::to(system.kinetic_energy(), "kJ/mol").unwrap();
        let total = units::to(system.total_energy(), "kJ/mol").unwrap();
        if let Err(e) = writeln!(&mut self.file, "{} {} {} {}", system.step(), potential, kinetic, total) {
            error!("Could not write to file '{}': {:?}", self.path, e);
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
    fn setup(&mut self, _: &System) {
        if let Err(e) = writeln!(&mut self.file, "# Physical properties of the simulation") {
            panic!("Could not write to file '{}': {:?}", self.path, e);
        }
        if let Err(e) = writeln!(&mut self.file, "# Step Volume/A^3 Temperature/K Pressure/bar") {
            panic!("Could not write to file '{}': {:?}", self.path, e);
        }
    }

    fn write(&mut self, system: &System) {
        let volume = units::to(system.volume(), "A^3").unwrap();
        let temperature = units::to(system.temperature(), "K").unwrap();
        let pressure = units::to(system.pressure(), "bar").unwrap();
        if let Err(e) = writeln!(&mut self.file, "{} {} {} {}", system.step(), volume, temperature, pressure) {
            error!("Could not write to file '{}': {:?}", self.path, e);
        }
    }
}

#[cfg(test)]
mod tests {
    use std::io::prelude::*;
    use std::fs;

    use super::*;
    use system::*;
    use types::*;
    use potentials::*;
    use units;

    fn testing_system() -> System {
        let mut system = System::from_cell(UnitCell::cubic(10.0));;

        system.add_particle(Particle::new("F"));
        system[0].position = Vector3D::new(0.0, 0.0, 0.0);

        system.add_particle(Particle::new("F"));
        system[1].position = Vector3D::new(1.3, 0.0, 0.0);

        system.add_pair_interaction("F", "F",
            Box::new(Harmonic{k: units::from(300.0, "kJ/mol/A^2").unwrap(), x0: units::from(1.2, "A").unwrap()}));
        return system;
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
        let system = testing_system();
        {
            let mut out = TrajectoryOutput::new(filename).unwrap();
            out.setup(&system);
            out.write(&system);
            out.finish(&system);
        }

        check_file_content(filename, "2\nWritten by the chemfiles library\nF 0 0 0\nF 1.3 0 0\n");
        fs::remove_file(filename).unwrap();
    }

    #[test]
    fn energy() {
        let filename = "testing-energy-output.dat";
        let system = testing_system();
        {
            let mut out = EnergyOutput::new(filename).unwrap();
            out.setup(&system);
            out.write(&system);
            out.finish(&system);
        }

        check_file_content(filename, "# Energy of the simulation (kJ/mol)\n# Step Potential Kinetic Total\n0 1.5000000000000027 0 1.5000000000000027\n");
        fs::remove_file(filename).unwrap();
    }

    #[test]
    fn cell() {
        let filename = "testing-cell-output.dat";
        let system = testing_system();
        {
            let mut out = CellOutput::new(filename).unwrap();
            out.setup(&system);
            out.write(&system);
            out.finish(&system);
        }

        check_file_content(filename, "# Unit cell of the simulation\n# Step A/Å B/Å C/Å α/deg β/deg γ/deg\n0 10 10 10 90 90 90\n");
        fs::remove_file(filename).unwrap();
    }

    #[test]
    fn properties() {
        let filename = "testing-properties-output.dat";
        let system = testing_system();
        {
            let mut out = PropertiesOutput::new(filename).unwrap();
            out.setup(&system);
            out.write(&system);
            out.finish(&system);
        }

        check_file_content(filename, "# Physical properties of the simulation\n# Step Volume/A^3 Temperature/K Pressure/bar\n0 1000 0 -431.7400836223091\n");
        fs::remove_file(filename).unwrap();
    }
}
