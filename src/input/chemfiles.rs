// Cymbalum, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

//! [Chemfiles](https://github.com/chemfiles/chemfiles/) adaptators for Cymbalum.
use system::{System, Particle, Molecule, UnitCell, CellShape};
use types::Vector3D;
use chemfiles;
use input::error::TrajectoryError;

use std::path::Path;

/// Convert chemfiles types to Cymbalum types
pub trait ToCymbalum {
    /// Output type
    type Output;
    /// Conversion function
    fn to_cymbalum(self) -> TrajectoryResult<Self::Output>;
}

impl ToCymbalum for chemfiles::Atom {
    type Output = Particle;
    fn to_cymbalum(self) -> TrajectoryResult<Particle> {
        let name = try!(self.name());
        let mut part = Particle::new(name);
        let mass = try!(self.mass());
        part.mass = mass as f64;
        Ok(part)
    }
}

impl ToCymbalum for chemfiles::UnitCell {
    type Output = UnitCell;
    fn to_cymbalum(self) -> TrajectoryResult<UnitCell> {
        let cell_type = try!(self.cell_type());
        let cell = match cell_type {
            chemfiles::CellType::Infinite => UnitCell::new(),
            chemfiles::CellType::Orthorhombic => {
                let (a, b, c) = try!(self.lengths());
                UnitCell::ortho(a, b, c)
            },
            chemfiles::CellType::Triclinic => {
                let (a, b, c) = try!(self.lengths());
                let (alpha, beta, gamma) = try!(self.angles());
                UnitCell::triclinic(a, b, c, alpha, beta, gamma)
            }
        };
        Ok(cell)
    }
}

impl ToCymbalum for chemfiles::Frame {
    type Output = System;
    fn to_cymbalum(self) -> TrajectoryResult<System> {
        let cell = try!(self.cell());
        let cell = try!(cell.to_cymbalum());
        let mut system = System::from_cell(cell);
        let topology = try!(self.topology());
        let natoms = try!(self.natoms());

        let positions = try!(self.positions());
        for i in 0..natoms {
            let atom = try!(topology.atom(i));
            let particle = try!(atom.to_cymbalum());

            system.add_particle(particle);
            let position = Vector3D::new(
                positions[i][0] as f64,
                positions[i][1] as f64,
                positions[i][2] as f64
            );
            system[i].position = position;
        }

        let mut bonds = try!(topology.bonds());
        while let Some(bond) = bonds.pop() {
            let permutations = system.add_bond(bond[0] as usize, bond[1] as usize);
            apply_particle_permutation(&mut bonds, permutations);
        }
        Ok(system)
    }
}

fn apply_particle_permutation(bonds: &mut Vec<[usize; 2]>, perms: Vec<(usize, usize)>) {
    for bond in bonds {
        // Search for a permutation applying to the first atom of the bond. We
        // need to stop just after the first permutations is found, because we
        // can have a permutation looking like this: [1 -> 2, 2 -> 3, 3 -> 4].
        // If we do not stop after the first match, then all indexes in 1-3
        // range will become 4.
        for perm in &perms {
            if bond[0] == perm.0 {
                bond[0] = perm.1;
                break;
            }
        }

        // Now we look for permutations applying to the second atom of the bond
        for perm in &perms {
            if bond[1] == perm.0 {
                bond[1] = perm.1;
                break;
            }
        }
    }
}

/******************************************************************************/

/// Convert Cymbalum types to chemfiles types
pub trait ToChemfiles {
    /// Output type
    type Output;
    /// Conversion function
    fn to_chemfiles(&self) -> TrajectoryResult<Self::Output>;
}

impl ToChemfiles for Particle {
    type Output = chemfiles::Atom;
    fn to_chemfiles(&self) -> TrajectoryResult<chemfiles::Atom> {
        let mut res = try!(chemfiles::Atom::new(self.name()));
        try!(res.set_mass(self.mass as f32));
        return Ok(res);
    }
}

impl ToChemfiles for UnitCell {
    type Output = chemfiles::UnitCell;
    fn to_chemfiles(&self) -> TrajectoryResult<chemfiles::UnitCell> {
        let res = match self.shape() {
            CellShape::Infinite => {
                try!(chemfiles::UnitCell::infinite())
            }
            CellShape::Orthorombic => {
                let (a, b, c) = (self.a(), self.b(), self.c());
                try!(chemfiles::UnitCell::new(a, b, c))
            },
            CellShape::Triclinic => {
                let (a, b, c) = (self.a(), self.b(), self.c());
                let (alpha, beta, gamma) = (self.alpha(), self.beta(), self.gamma());
                try!(chemfiles::UnitCell::triclinic(a, b, c, alpha, beta, gamma))
            },
        };
        return Ok(res);
    }
}

impl ToChemfiles for System {
    type Output = chemfiles::Frame;
    fn to_chemfiles(&self) -> TrajectoryResult<chemfiles::Frame> {

        let natoms = self.size();
        let mut frame = try!(chemfiles::Frame::new(natoms));

        try!(frame.set_step(self.step() as usize));

        {
            let positions = try!(frame.positions_mut());
            for (i, p) in self.iter().enumerate() {
                let pos = p.position;
                positions[i][0] = pos[0] as f32;
                positions[i][1] = pos[1] as f32;
                positions[i][2] = pos[2] as f32;
            }
        }

        {
            try!(frame.add_velocities());
            let velocities = try!(frame.velocities_mut());
            for (i, p) in self.iter().enumerate() {
                let vel = p.velocity;
                velocities[i][0] = vel[0] as f32;
                velocities[i][1] = vel[1] as f32;
                velocities[i][2] = vel[2] as f32;
            }
        }

        let mut topology = try!(chemfiles::Topology::new());
        for particle in self {
            let atom = try!(particle.to_chemfiles());
            try!(topology.push(&atom));
        }


        for molecule in self.molecules() {
            for bond in molecule.bonds() {
                try!(topology.add_bond(bond.i(), bond.j()));
            }
        }

        try!(frame.set_topology(&topology));
        let cell = try!(self.cell().to_chemfiles());
        try!(frame.set_cell(&cell));
        Ok(frame)
    }
}

/******************************************************************************/

/// A Trajectory is a file containing one or more successives simulation steps
pub struct Trajectory(chemfiles::Trajectory);

/// Result type for all Trajectory operations
pub type TrajectoryResult<T> = Result<T, TrajectoryError>;

impl Trajectory {
    /// Open an existing file at `path` for reading.
    pub fn open<P: AsRef<Path>>(path: P) -> TrajectoryResult<Trajectory> {
        let traj = try!(chemfiles::Trajectory::open(path));
        return Ok(Trajectory(traj));
    }

    /// Create a new file at `path` for writing, and overwrite any existing file.
    pub fn create<P: AsRef<Path>>(path: P) -> TrajectoryResult<Trajectory> {
        let traj = try!(chemfiles::Trajectory::create(path));
        return Ok(Trajectory(traj));
    }

    /// Read the next step of the trajectory
    pub fn read(&mut self) -> TrajectoryResult<System> {
        let mut frame = try!(chemfiles::Frame::new(0));
        try!(self.0.read(&mut frame));
        return frame.to_cymbalum();
    }

    /// Read the next step of the trajectory, and guess the bonds of the
    /// resulting System.
    pub fn read_guess_bonds(&mut self) -> TrajectoryResult<System> {
        let mut frame = try!(chemfiles::Frame::new(0));
        try!(self.0.read(&mut frame));
        try!(frame.guess_topology());
        return frame.to_cymbalum();
    }

    /// Write the system to the trajectory.
    pub fn write(&mut self, system: &System) -> TrajectoryResult<()> {
        let frame = try!(system.to_chemfiles());
        try!(self.0.write(&frame));
        Ok(())
    }

    /// Get access to the chemfiles trajectory, and the associated features
    // TODO: use partial privacy for this function
    pub fn as_chemfiles(&mut self) -> &mut chemfiles::Trajectory {
        &mut self.0
    }
}

/// Read a the first molecule from the file at `path`. If no bond information
/// exists in the file, bonds are guessed.
pub fn read_molecule<P: AsRef<Path>>(path: P) -> TrajectoryResult<(Molecule, Vec<Particle>)> {
    let path = path.as_ref();

    let mut trajectory = try!(chemfiles::Trajectory::open(path));
    let mut frame = try!(chemfiles::Frame::new(0));
    try!(trajectory.read(&mut frame));

    // Only guess the topology when we have no bond information
    let topology = try!(frame.topology());
    if try!(topology.bonds_count()) == 0 {
        try!(frame.guess_topology());
    }
    let system = try!(frame.to_cymbalum());

    assert!(system.size() != 0, "No molecule in the file at {}", path.display());
    let molecule = system.molecule(0).clone();
    let mut particles = Vec::new();
    for i in &molecule {
        particles.push(system[i].clone());
    }
    return Ok((molecule, particles));
}

/// Guess the bonds in a system
pub fn guess_bonds(system: System) -> TrajectoryResult<System> {
    let mut frame = try!(system.to_chemfiles());
    try!(frame.guess_topology());
    return frame.to_cymbalum();
}


#[cfg(test)]
mod tests {
    extern crate tempfile;
    use self::tempfile::NamedTempFileOptions;

    use super::*;
    use std::io::prelude::*;
    use system::moltype;

    static MOLECULE: &'static str = "3

O 0.0 0.0 0.0
H 1.0 0.0 0.0
H 0.0 1.0 0.0
";

    #[test]
    fn read() {
        let mut file = NamedTempFileOptions::new().suffix(".xyz").create().unwrap();
        write!(file, "{}", MOLECULE).unwrap();

        let (molecule, atoms) = read_molecule(file.path()).unwrap();

        assert_eq!(molecule.size(), 3);
        assert_eq!(molecule.size(), atoms.len());
        assert_eq!(molecule.bonds().len(), 2);

        assert_eq!(atoms[0].name(), "O");
        assert_eq!(atoms[1].name(), "H");
        assert_eq!(atoms[2].name(), "H");

        // This is only a simple regression test on the moltype function. Feel
        // free to change the value if the molecule type algorithm change.
        assert_eq!(moltype(&molecule, &atoms), 4144180246440175497);
    }
}
