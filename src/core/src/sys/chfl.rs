// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

//! [Chemfiles][Chemfiles] conversion for Lumol.
//!
//! [Chemfiles]: https://chemfiles.org/
use chemfiles;
use sys::{CellShape, Molecule, Particle, ParticleRef, ParticleVec, System, UnitCell};
use types::Vector3D;

use std::error;
use std::fmt;
use std::path::Path;
use std::sync::{Once, ONCE_INIT};

/// Possible error when reading and writing to trajectories
#[derive(Debug)]
pub struct TrajectoryError(chemfiles::Error);

impl From<chemfiles::Error> for TrajectoryError {
    fn from(err: chemfiles::Error) -> TrajectoryError {
        TrajectoryError(err)
    }
}

impl fmt::Display for TrajectoryError {
    fn fmt(&self, fmt: &mut fmt::Formatter) -> Result<(), fmt::Error> {
        self.0.fmt(fmt)
    }
}

impl error::Error for TrajectoryError {
    fn description(&self) -> &str {
        self.0.description()
    }
}

/// Convert Chemfiles types to Lumol types
pub trait ToLumol {
    /// Output type
    type Output;
    /// Conversion function
    fn to_lumol(self) -> TrajectoryResult<Self::Output>;
}

impl ToLumol for chemfiles::Atom {
    type Output = Particle;
    fn to_lumol(self) -> TrajectoryResult<Particle> {
        let name = try!(self.name());
        let mut part = Particle::new(name);
        let mass = try!(self.mass());
        part.mass = mass as f64;
        Ok(part)
    }
}

impl ToLumol for chemfiles::UnitCell {
    type Output = UnitCell;
    fn to_lumol(self) -> TrajectoryResult<UnitCell> {
        let cell_type = try!(self.shape());
        let cell = match cell_type {
            chemfiles::CellShape::Infinite => UnitCell::infinite(),
            chemfiles::CellShape::Orthorhombic => {
                let lengths = try!(self.lengths());
                UnitCell::ortho(lengths[0], lengths[1], lengths[2])
            }
            chemfiles::CellShape::Triclinic => {
                let lengths = try!(self.lengths());
                let angles = try!(self.angles());
                UnitCell::triclinic(
                    lengths[0], lengths[1], lengths[2],
                    angles[0], angles[1], angles[2]
                )
            }
        };
        Ok(cell)
    }
}

impl ToLumol for chemfiles::Frame {
    type Output = System;
    fn to_lumol(self) -> TrajectoryResult<System> {
        let cell = try!(self.cell());
        let cell = try!(cell.to_lumol());
        let mut system = System::with_cell(cell);
        let topology = try!(self.topology());
        let natoms = try!(self.size()) as usize;

        let positions = try!(self.positions());
        for i in 0..natoms {
            let atom = try!(topology.atom(i as u64));
            let particle = try!(atom.to_lumol());

            system.add_particle(particle);
            let position = Vector3D::new(positions[i][0], positions[i][1], positions[i][2]);
            system.particles_mut().position[i] = position;
        }

        let mut bonds = try!(topology.bonds());
        while let Some(bond) = bonds.pop() {
            let permutations = system.add_bond(bond[0] as usize, bond[1] as usize);
            apply_particle_permutation(&mut bonds, &permutations);
        }
        Ok(system)
    }
}

fn apply_particle_permutation(bonds: &mut Vec<[u64; 2]>, permutations: &[(usize, usize)]) {
    for bond in bonds {
        // Search for a permutation applying to the first atom of the bond. We
        // need to stop just after the first permutations is found, because we
        // can have a permutation looking like this: [1 -> 2, 2 -> 3, 3 -> 4].
        // If we do not stop after the first match, then all indexes in 1-3
        // range will become 4.
        for permutation in permutations {
            if bond[0] == permutation.0 as u64 {
                bond[0] = permutation.1 as u64;
                break;
            }
        }

        // Now we look for permutations applying to the second atom of the bond
        for permutation in permutations {
            if bond[1] == permutation.0 as u64 {
                bond[1] = permutation.1 as u64;
                break;
            }
        }
    }
}

/// Convert Lumol types to Chemfiles types
pub trait ToChemfiles {
    /// Output type
    type Output;
    /// Conversion function
    fn to_chemfiles(&self) -> TrajectoryResult<Self::Output>;
}

impl<'a> ToChemfiles for ParticleRef<'a> {
    type Output = chemfiles::Atom;
    fn to_chemfiles(&self) -> TrajectoryResult<chemfiles::Atom> {
        let mut atom = try!(chemfiles::Atom::new(&**self.name));
        try!(atom.set_mass(*self.mass));
        return Ok(atom);
    }
}

impl ToChemfiles for UnitCell {
    type Output = chemfiles::UnitCell;
    fn to_chemfiles(&self) -> TrajectoryResult<chemfiles::UnitCell> {
        let res = match self.shape() {
            CellShape::Infinite => try!(chemfiles::UnitCell::infinite()),
            CellShape::Orthorhombic => {
                let lengths = [self.a(), self.b(), self.c()];
                try!(chemfiles::UnitCell::new(lengths))
            }
            CellShape::Triclinic => {
                let lengths = [self.a(), self.b(), self.c()];
                let angles = [self.alpha(), self.beta(), self.gamma()];
                try!(chemfiles::UnitCell::triclinic(lengths, angles))
            }
        };
        return Ok(res);
    }
}

impl ToChemfiles for System {
    type Output = chemfiles::Frame;
    fn to_chemfiles(&self) -> TrajectoryResult<chemfiles::Frame> {
        let mut frame = try!(chemfiles::Frame::new());
        try!(frame.resize(self.size() as u64));
        try!(frame.set_step(self.step() as u64));

        {
            let chfl_positions = try!(frame.positions_mut());
            let positions = self.particles().position;
            for (chfl_position, position) in izip!(chfl_positions, positions) {
                chfl_position[0] = position[0];
                chfl_position[1] = position[1];
                chfl_position[2] = position[2];
            }
        }

        {
            try!(frame.add_velocities());
            let chfl_velocities = try!(frame.velocities_mut());
            let velocities = self.particles().velocity;
            for (chfl_velocity, velocity) in izip!(chfl_velocities, velocities) {
                chfl_velocity[0] = velocity[0];
                chfl_velocity[1] = velocity[1];
                chfl_velocity[2] = velocity[2];
            }
        }

        let mut topology = try!(chemfiles::Topology::new());
        for particle in self.particles().iter() {
            let atom = try!(particle.to_chemfiles());
            try!(topology.add_atom(&atom));
        }

        for molecule in self.molecules() {
            for bond in molecule.bonds() {
                try!(topology.add_bond(bond.i() as u64, bond.j() as u64));
            }
        }

        try!(frame.set_topology(&topology));
        let cell = try!(self.cell.to_chemfiles());
        try!(frame.set_cell(&cell));
        Ok(frame)
    }
}

/// A Trajectory is a file containing one or more successive simulation steps.
///
/// One should use the [`TrajectoryBuilder`](struct.TrajectoryBuilder.html) to
/// create a new trajectory.
///
/// # Examples
///
/// ```no_run
/// # use lumol_core::sys::TrajectoryBuilder;
/// let mut trajectory = TrajectoryBuilder::new()
///                                       .open("file.xyz")
///                                       .unwrap();
///
/// let system = trajectory.read().unwrap();
/// ```
pub struct Trajectory(chemfiles::Trajectory);

/// Result type for all Trajectory operations
type TrajectoryResult<T> = Result<T, TrajectoryError>;

/// Possible modes when opening a [`Trajectory`](struct.Trajectory.html).
pub enum OpenMode {
    /// Open the file as read-only
    Read,
    /// Open the file as write-only, and overwrite any existing file
    Write,
    /// Open the file as read-write, keeping existing files
    Append,
}

/// A [`Trajectory`](struct.Trajectory.html) builder, to set some options
/// before opening a trajectory.
///
/// # Examples
///
/// ```no_run
/// # use lumol_core::sys::{TrajectoryBuilder, OpenMode};
/// let trajectory = TrajectoryBuilder::new()
///                                    .mode(OpenMode::Read)
///                                    .open("file.xyz")
///                                    .unwrap();
/// ```
pub struct TrajectoryBuilder<'a> {
    mode: OpenMode,
    format: &'a str,
}

impl<'a> TrajectoryBuilder<'a> {
    /// Create a new builder in read mode and with automatic format detection.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use lumol_core::sys::{TrajectoryBuilder, OpenMode};
    /// let trajectory = TrajectoryBuilder::new()
    ///                                    .open("file.xyz")
    ///                                    .unwrap();
    /// ```
    pub fn new() -> TrajectoryBuilder<'a> {
        TrajectoryBuilder {
            mode: OpenMode::Read,
            format: "",
        }
    }

    /// Use a specific `format` when opening a file. See the [chemfiles]
    /// documentation for a format list.
    ///
    /// [chemfiles]: http://chemfiles.org/chemfiles/latest/formats.html#list-of-supported-formats
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use lumol_core::sys::TrajectoryBuilder;
    /// let trajectory = TrajectoryBuilder::new()
    ///                                    .format("PDB")
    ///                                    .open("file.mol")
    ///                                    .unwrap();
    /// ```
    pub fn format(self, format: &'a str) -> TrajectoryBuilder<'a> {
        TrajectoryBuilder {
            format: format,
            mode: self.mode,
        }
    }

    /// Use a specific `mode` when opening a file.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use lumol_core::sys::{TrajectoryBuilder, OpenMode};
    /// let trajectory = TrajectoryBuilder::new()
    ///                                    .mode(OpenMode::Write)
    ///                                    .open("file.nc")
    ///                                    .unwrap();
    /// ```
    pub fn mode(self, mode: OpenMode) -> TrajectoryBuilder<'a> {
        TrajectoryBuilder {
            mode: mode,
            format: self.format,
        }
    }

    /// Open the trajectory at the given `path`.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use lumol_core::sys::TrajectoryBuilder;
    /// let trajectory = TrajectoryBuilder::new()
    ///                                    .open("file.nc")
    ///                                    .unwrap();
    /// ```
    pub fn open<P: AsRef<Path>>(self, path: P) -> TrajectoryResult<Trajectory> {
        redirect_chemfiles_warnings();
        let mode = match self.mode {
            OpenMode::Read => 'r',
            OpenMode::Write => 'w',
            OpenMode::Append => 'a',
        };
        let trajectory = try!(chemfiles::Trajectory::open_with_format(path, mode, self.format));
        return Ok(Trajectory(trajectory));
    }
}

impl Trajectory {
    /// Read the next step of the trajectory
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use lumol_core::sys::TrajectoryBuilder;
    /// let mut trajectory = TrajectoryBuilder::new()
    ///                                       .open("file.nc")
    ///                                       .unwrap();
    ///
    /// let system = trajectory.read().unwrap();
    /// ```
    pub fn read(&mut self) -> TrajectoryResult<System> {
        let mut frame = try!(chemfiles::Frame::new());
        try!(self.0.read(&mut frame));
        return frame.to_lumol();
    }

    /// Read the next step of the trajectory, and guess the bonds of the
    /// resulting [`System`][struct.System.html].
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use lumol_core::sys::TrajectoryBuilder;
    /// let mut trajectory = TrajectoryBuilder::new()
    ///                                       .open("file.nc")
    ///                                       .unwrap();
    ///
    /// let system = trajectory.read_guess_bonds().unwrap();
    /// ```
    pub fn read_guess_bonds(&mut self) -> TrajectoryResult<System> {
        let mut frame = try!(chemfiles::Frame::new());
        try!(self.0.read(&mut frame));
        try!(frame.guess_topology());
        return frame.to_lumol();
    }

    /// Write the system to the trajectory.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use lumol_core::sys::{System, TrajectoryBuilder, OpenMode};
    /// # let system = System::new();
    /// let mut trajectory = TrajectoryBuilder::new()
    ///                                       .mode(OpenMode::Write)
    ///                                       .open("file.xyz")
    ///                                       .unwrap();
    ///
    /// trajectory.write(&system).unwrap();
    /// ```
    pub fn write(&mut self, system: &System) -> TrajectoryResult<()> {
        let frame = try!(system.to_chemfiles());
        try!(self.0.write(&frame));
        Ok(())
    }

    /// Set the unit cell associated with a trajectory. This cell will be used
    /// when reading and writing the files, replacing any unit cell in the
    /// frames or files.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use lumol_core::sys::{TrajectoryBuilder, UnitCell};
    /// let mut trajectory = TrajectoryBuilder::new()
    ///                                       .open("file.xyz")
    ///                                       .unwrap();
    ///
    /// trajectory.set_cell(&UnitCell::cubic(10.0)).unwrap();
    /// let system = trajectory.read().unwrap();
    ///
    /// assert_eq!(system.cell, UnitCell::cubic(10.0));
    /// ```
    pub fn set_cell(&mut self, cell: &UnitCell) -> TrajectoryResult<()> {
        let cell = try!(cell.to_chemfiles());
        try!(self.0.set_cell(&cell));
        Ok(())
    }

    /// Set the topology associated with this trajectory by reading the first
    /// frame of the file at the given `path` and extracting the topology of
    /// this frame. This topology will be used to replace any existing topology
    /// when reading or writing with this trajectory.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use lumol_core::sys::TrajectoryBuilder;
    /// let mut trajectory = TrajectoryBuilder::new()
    ///                                       .open("file.xyz")
    ///                                       .unwrap();
    ///
    /// trajectory.set_topology_file("topology.pdb").unwrap();
    /// // The system will contain the topology from topology.pdb
    /// let system = trajectory.read().unwrap();
    /// ```
    pub fn set_topology_file(&mut self, path: &str) -> TrajectoryResult<()> {
        try!(self.0.set_topology_file(path));
        Ok(())
    }
}

/// Read a the first molecule from the file at `path`. If no bond information
/// exists in the file, bonds are guessed.
///
/// # Examples
///
/// ```no_run
/// # use lumol_core::sys::read_molecule;
/// let (molecule, particles) = read_molecule("file.xyz").unwrap();
/// assert_eq!(molecule.size(), particles.len());
/// ```
pub fn read_molecule<P: AsRef<Path>>(path: P) -> TrajectoryResult<(Molecule, ParticleVec)> {
    let mut trajectory = try!(chemfiles::Trajectory::open(&path, 'r'));
    let mut frame = try!(chemfiles::Frame::new());
    try!(trajectory.read(&mut frame));

    // Only guess the topology when we have no bond information
    let topology = try!(frame.topology());
    if try!(topology.bonds_count()) == 0 {
        try!(frame.guess_topology());
    }
    let system = try!(frame.to_lumol());

    assert!(!system.is_empty(), "No molecule in the file at {}", path.as_ref().display());
    let molecule = system.molecule(0).clone();
    return Ok((molecule, system.particles().to_vec()));
}

static REDIRECT_CHEMFILES_WARNING: Once = ONCE_INIT;

fn redirect_chemfiles_warnings() {
    fn warning_callback(message: &str) {
        warn!("[chemfiles] {}", message);
    }

    REDIRECT_CHEMFILES_WARNING.call_once(|| {
        chemfiles::set_warning_callback(warning_callback).expect("could not redirect chemfiles warning");
    });
}


#[cfg(test)]
mod tests {
    extern crate tempfile;

    use super::*;
    use std::io::prelude::*;
    use sys::{Angle, Bond};
    use sys::molecule_type;

    static WATER: &'static str = "3

O 0.0 0.0 0.0
H 1.0 0.0 0.0
H 0.0 1.0 0.0
";

    static PROPANE: &'static str = "11

H     -1.306761     0.917434     0.885782
C     -1.277333     0.278101     0.000000
H     -2.172669    -0.348524    -0.000126
H     -1.306638     0.917616    -0.885655
C      0.000000    -0.596258     0.000000
H      0.000000    -1.247870    -0.879715
H      0.000000    -1.247870     0.879715
C      1.277333     0.278101     0.000000
H      1.306675     0.917562     0.885693
H      1.306724     0.917489    -0.885744
H      2.172669     -0.348524    0.000051
";

    #[test]
    fn read_water() {
        let mut file = tempfile::Builder::new().suffix(".xyz").tempfile().unwrap();
        write!(file, "{}", WATER).unwrap();

        let (molecule, atoms) = read_molecule(file.path()).unwrap();

        assert_eq!(molecule.size(), 3);
        assert_eq!(molecule.size(), atoms.len());

        assert_eq!(molecule.bonds().len(), 2);
        assert!(molecule.bonds().contains(&Bond::new(0, 1)));
        assert!(molecule.bonds().contains(&Bond::new(0, 2)));

        assert_eq!(molecule.angles().len(), 1);
        assert!(molecule.angles().contains(&Angle::new(1, 0, 2)));

        assert!(molecule.dihedrals().is_empty());

        assert_eq!(atoms.name[0], "O");
        assert_eq!(atoms.name[1], "H");
        assert_eq!(atoms.name[2], "H");

        // This is only a simple regression test on the moltype function. Feel
        // free to change the value if the molecule type algorithm change.
        assert_eq!(molecule_type(&molecule, atoms.as_slice()), 2727145596042757306);
    }

    #[test]
    fn read_propane() {
        let mut file = tempfile::Builder::new().suffix(".xyz").tempfile().unwrap();
        write!(file, "{}", PROPANE).unwrap();

        let (molecule, atoms) = read_molecule(file.path()).unwrap();

        assert_eq!(molecule.size(), 11);
        assert_eq!(molecule.size(), atoms.len());

        assert_eq!(molecule.bonds().len(), 10);
        assert_eq!(molecule.angles().len(), 18);
        assert_eq!(molecule.dihedrals().len(), 18);

        // This is only a simple regression test on the moltype function. Feel
        // free to change the value if the molecule type algorithm change.
        assert_eq!(molecule_type(&molecule, atoms.as_slice()), 2028056351119909064);
    }
}
