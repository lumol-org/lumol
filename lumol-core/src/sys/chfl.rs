// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

//! [Chemfiles](https://chemfiles.org/) conversion for Lumol.
use std::path::Path;
use std::sync::Once;

use soa_derive::soa_zip;
use log::warn;

use crate::sys::Permutation;
use crate::{Molecule, Particle, ParticleRef, System, UnitCell, CellShape};
use crate::Vector3D;

impl<'a> From<&'a chemfiles::Atom> for Particle {
    fn from(atom: &'a chemfiles::Atom) -> Particle {
        let name = atom.atomic_type();
        let mut particle = Particle::new(name);
        particle.mass = atom.mass();
        return particle;
    }
}

impl<'a> From<&'a chemfiles::UnitCell> for UnitCell {
    fn from(cell: &'a chemfiles::UnitCell) -> UnitCell {
        match cell.shape() {
            chemfiles::CellShape::Infinite => UnitCell::infinite(),
            chemfiles::CellShape::Orthorhombic => {
                let lengths = cell.lengths();
                UnitCell::ortho(lengths[0], lengths[1], lengths[2])
            }
            chemfiles::CellShape::Triclinic => {
                let lengths = cell.lengths();
                let angles = cell.angles();
                UnitCell::triclinic(
                    lengths[0], lengths[1], lengths[2],
                    angles[0], angles[1], angles[2]
                )
            }
        }
    }
}

impl From<chemfiles::Frame> for System {
    fn from(frame: chemfiles::Frame) -> System {
        let cell = UnitCell::from(&*frame.cell());
        let mut system = System::with_cell(cell);

        for (i, position) in frame.positions().iter().enumerate() {
            let mut particle = Particle::from(&*frame.atom(i));
            particle.position = Vector3D::from(*position);

            system.add_molecule(Molecule::new(particle));
        }


        if let Some(velocities) = frame.velocities() {
            for (i, velocity) in velocities.iter().enumerate() {
                system.particles_mut().velocity[i] = Vector3D::from(*velocity);
            }
        }

        let mut bonds = frame.topology().bonds();
        while let Some(bond) = bonds.pop() {
            let permutations = system.add_bond(bond[0], bond[1]);
            apply_particle_permutation(&mut bonds, &permutations);
        }
        return system;
    }
}

fn apply_particle_permutation(bonds: &mut Vec<[usize; 2]>, permutations: &[Permutation]) {
    for bond in bonds {
        // Search for a permutation applying to the first atom of the bond. We
        // need to stop just after the first permutations is found, because we
        // can have a permutation looking like this: [1 -> 2, 2 -> 3, 3 -> 4].
        // If we do not stop after the first match, then all indexes in 1-3
        // range will become 4.
        for permutation in permutations {
            if bond[0] == permutation.old {
                bond[0] = permutation.new;
                break;
            }
        }

        // Now we look for permutations applying to the second atom of the bond
        for permutation in permutations {
            if bond[1] == permutation.old {
                bond[1] = permutation.new;
                break;
            }
        }
    }
}

impl<'a> From<ParticleRef<'a>> for chemfiles::Atom {
    fn from(particle: ParticleRef<'a>) -> chemfiles::Atom {
        let mut atom = chemfiles::Atom::new(&**particle.name);
        atom.set_mass(*particle.mass);
        return atom;
    }
}

impl<'a> From<&'a UnitCell> for chemfiles::UnitCell {
    fn from(cell: &'a UnitCell) -> chemfiles::UnitCell {
        match cell.shape() {
            CellShape::Infinite => chemfiles::UnitCell::infinite(),
            CellShape::Orthorhombic => {
                let lengths = [cell.a(), cell.b(), cell.c()];
                chemfiles::UnitCell::new(lengths)
            }
            CellShape::Triclinic => {
                let lengths = [cell.a(), cell.b(), cell.c()];
                let angles = [cell.alpha(), cell.beta(), cell.gamma()];
                chemfiles::UnitCell::triclinic(lengths, angles)
            }
        }
    }
}

impl<'a> From<&'a System> for chemfiles::Frame {
    fn from(system: &'a System) -> chemfiles::Frame {
        let mut frame = chemfiles::Frame::new();
        frame.resize(system.size());
        frame.set_step(system.step as usize);

        for (position, chfl_position) in soa_zip!(system.particles(), [position], frame.positions_mut()) {
            *chfl_position = **position;
        }

        frame.add_velocities();
        for (velocity, chfl_velocity) in soa_zip!(system.particles(), [position], frame.positions_mut()) {
            *chfl_velocity = **velocity;
        }

        let mut topology = chemfiles::Topology::new();
        for particle in system.particles() {
            topology.add_atom(&particle.into());
        }

        for molecule in system.molecules() {
            for bond in molecule.bonds() {
                topology.add_bond(bond.i(), bond.j());
            }
        }

        frame.set_topology(&topology).expect("wrong topology size");
        let cell = (&system.cell).into();
        frame.set_cell(&cell);
        return frame;
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
///     .open("file.xyz")
///     .unwrap();
///
/// let system = trajectory.read().unwrap();
/// ```
pub struct Trajectory(chemfiles::Trajectory);

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
///     .mode(OpenMode::Read)
///     .open("file.xyz")
///     .unwrap();
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
    /// # use lumol_core::sys::TrajectoryBuilder;
    /// let trajectory = TrajectoryBuilder::new()
    ///     .open("file.xyz")
    ///     .unwrap();
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
    ///     .format("PDB")
    ///     .open("file.mol")
    ///     .unwrap();
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
    ///     .mode(OpenMode::Write)
    ///     .open("file.nc")
    ///     .unwrap();
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
    ///     .open("file.nc")
    ///     .unwrap();
    /// ```
    pub fn open<P: AsRef<Path>>(self, path: P) -> Result<Trajectory, chemfiles::Error> {
        redirect_chemfiles_warnings();
        let mode = match self.mode {
            OpenMode::Read => 'r',
            OpenMode::Write => 'w',
            OpenMode::Append => 'a',
        };
        let trajectory = chemfiles::Trajectory::open_with_format(path, mode, self.format)?;
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
    ///     .open("file.nc")
    ///     .unwrap();
    ///
    /// let system = trajectory.read().unwrap();
    /// ```
    pub fn read(&mut self) -> Result<System, chemfiles::Error> {
        let mut frame = chemfiles::Frame::new();
        self.0.read(&mut frame)?;
        return Ok(frame.into());
    }

    /// Read the next step of the trajectory, and guess the bonds of the
    /// resulting [`System`][struct.System.html].
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use lumol_core::sys::TrajectoryBuilder;
    /// let mut trajectory = TrajectoryBuilder::new()
    ///     .open("file.nc")
    ///     .unwrap();
    ///
    /// let system = trajectory.read_guess_bonds().unwrap();
    /// ```
    pub fn read_guess_bonds(&mut self) -> Result<System, chemfiles::Error> {
        let mut frame = chemfiles::Frame::new();
        self.0.read(&mut frame)?;
        frame.guess_bonds()?;
        return Ok(frame.into());
    }

    /// Write the system to the trajectory.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use lumol_core::sys::{System, TrajectoryBuilder, OpenMode};
    /// # let system = System::new();
    /// let mut trajectory = TrajectoryBuilder::new()
    ///     .mode(OpenMode::Write)
    ///     .open("file.xyz")
    ///     .unwrap();
    ///
    /// trajectory.write(&system).unwrap();
    /// ```
    pub fn write(&mut self, system: &System) -> Result<(), chemfiles::Error> {
        self.0.write(&system.into())
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
    ///     .open("file.xyz")
    ///     .unwrap();
    ///
    /// trajectory.set_cell(&UnitCell::cubic(10.0));
    /// let system = trajectory.read().unwrap();
    ///
    /// assert_eq!(system.cell, UnitCell::cubic(10.0));
    /// ```
    pub fn set_cell(&mut self, cell: &UnitCell) {
        self.0.set_cell(&cell.into());
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
    ///     .open("file.xyz")
    ///     .unwrap();
    ///
    /// trajectory.set_topology_file("topology.pdb").unwrap();
    /// // The system will contain the topology from topology.pdb
    /// let system = trajectory.read().unwrap();
    /// ```
    pub fn set_topology_file(&mut self, path: &str) -> Result<(), chemfiles::Error> {
        self.0.set_topology_file(path)?;
        Ok(())
    }
}

/// Read a the first molecule from the file at `path`. If no bond information
/// exists in the file, bonds are guessed.
pub fn read_molecule<P: AsRef<Path>>(path: P) -> Result<Molecule, chemfiles::Error> {
    let mut trajectory = chemfiles::Trajectory::open(&path, 'r')?;
    let mut frame = chemfiles::Frame::new();
    trajectory.read(&mut frame)?;

    // Only guess the topology when we have no bond information
    if frame.topology().bonds_count() == 0 {
        frame.guess_bonds()?;
    }

    let system: System = frame.into();
    assert!(!system.is_empty(), "No molecule in the file at {}", path.as_ref().display());

    return Ok(system.molecule(0).to_owned());
}

static REDIRECT_CHEMFILES_WARNING: Once = Once::new();

fn redirect_chemfiles_warnings() {
    fn warning_callback(message: &str) {
        warn!("[chemfiles] {}", message);
    }

    REDIRECT_CHEMFILES_WARNING.call_once(|| {
        chemfiles::set_warning_callback(warning_callback);
    });
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Angle, Bond, MoleculeHash};
    use std::io::prelude::*;

    static WATER: &str = "3

O 0.0 0.0 0.0
H 1.0 0.0 0.0
H 0.0 1.0 0.0
";

    static PROPANE: &str = "11

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

    static PDB_WATER: &str = "
CRYST1   28.000   28.000   28.000  90.00  90.00  90.00 P 1           1
HETATM    1  O   RES X   1       0.000   0.000   0.000  1.00  0.00           O
HETATM    2  Hw  RES X   2       0.635   0.898   0.000  1.00  0.00           H
HETATM    3  Hw  RES X   3      -0.635   0.898   0.000  1.00  0.00           H
HETATM    4  Ow  RES X   4       0.000   0.000   9.333  1.00  0.00           O
HETATM    5  H   RES X   5       0.635   0.898   9.333  1.00  0.00           H
HETATM    6  H   RES X   6      -0.635   0.898   9.333  1.00  0.00           H
CONECT    1    2    3
CONECT    2    1
CONECT    3    1
END
";

    #[test]
    #[allow(clippy::unreadable_literal)]
    fn read_water() {
        let mut file = tempfile::Builder::new().suffix(".xyz").tempfile().unwrap();
        write!(file, "{WATER}").unwrap();

        let molecule = read_molecule(file.path()).unwrap();

        assert_eq!(molecule.size(), 3);

        assert_eq!(molecule.bonds().len(), 2);
        assert!(molecule.bonds().contains(&Bond::new(0, 1)));
        assert!(molecule.bonds().contains(&Bond::new(0, 2)));

        assert_eq!(molecule.angles().len(), 1);
        assert!(molecule.angles().contains(&Angle::new(1, 0, 2)));

        assert!(molecule.dihedrals().is_empty());

        assert_eq!(molecule.particles().name[0], "O");
        assert_eq!(molecule.particles().name[1], "H");
        assert_eq!(molecule.particles().name[2], "H");

        // This is only a simple regression test on the moltype function. Feel
        // free to change the value if the molecule type algorithm change.
        assert_eq!(molecule.hash(), MoleculeHash::new(3988311241583852942));
    }

    #[test]
    fn read_pdb_water() {
        let mut file = tempfile::Builder::new().suffix(".pdb").tempfile().unwrap();
        write!(file, "{PDB_WATER}").unwrap();

        let system = TrajectoryBuilder::new()
            .open(&file).unwrap()
            .read().unwrap();

        assert_eq!(system.size(), 6);
        assert_eq!(system.molecules().count(), 4);
        assert_eq!(system.cell, UnitCell::cubic(28.0));

        let molecule = system.molecule(0);
        assert_eq!(molecule.bonds().len(), 2);
        assert!(molecule.bonds().contains(&Bond::new(0, 1)));
        assert!(molecule.bonds().contains(&Bond::new(0, 2)));

        assert_eq!(system.particles().name[0], "O");
        assert_eq!(system.particles().name[1], "H");
        assert_eq!(system.particles().name[2], "H");
        assert_eq!(system.particles().name[3], "O");
        assert_eq!(system.particles().name[4], "H");
        assert_eq!(system.particles().name[5], "H");
    }

    #[test]
    #[allow(clippy::unreadable_literal)]
    fn read_propane() {
        let mut file = tempfile::Builder::new().suffix(".xyz").tempfile().unwrap();
        write!(file, "{PROPANE}").unwrap();

        let molecule = read_molecule(file.path()).unwrap();

        assert_eq!(molecule.size(), 11);
        assert_eq!(molecule.bonds().len(), 10);
        assert_eq!(molecule.angles().len(), 18);
        assert_eq!(molecule.dihedrals().len(), 18);

        // This is only a simple regression test on the moltype function. Feel
        // free to change the value if the molecule type algorithm change.
        assert_eq!(molecule.hash(), MoleculeHash::new(10634064187773497961));
    }
}
