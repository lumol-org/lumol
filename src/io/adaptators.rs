// Cymbalum, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

//! [Chemfiles](https://github.com/chemfiles/chemfiles/) adaptators for Cymbalum.
use system::{System, Particle, UnitCell, CellType};
use types::Vector3D;
use chemfiles;

use super::TrajectoryResult;

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
            chemfiles::CellType::Orthorombic => {
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
            if let Some(perms) = system.add_bond(bond[0] as usize, bond[1] as usize) {
                apply_particle_permutation(&mut bonds, perms);
            }
        }
        Ok(system)
    }
}

fn apply_particle_permutation(bonds: &mut Vec<[usize; 2]>, perms: Vec<(usize, usize)>) {
    'bonds: for bond in bonds {
        for perm in &perms {
            if bond[0] == perm.0 {
                bond[0] = perm.1;
                continue 'bonds;
            } else if bond[1] == perm.0 {
                bond[1] = perm.1;
                continue 'bonds;
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
        let res = match self.celltype() {
            CellType::Infinite => {
                try!(chemfiles::UnitCell::infinite())
            }
            CellType::Orthorombic => {
                let (a, b, c) = (self.a(), self.b(), self.c());
                try!(chemfiles::UnitCell::new(a, b, c))
            },
            CellType::Triclinic => {
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
        // Guessing angles and dihedrals
        try!(frame.guess_topology(false));

        let cell = try!(self.cell().to_chemfiles());
        try!(frame.set_cell(&cell));
        Ok(frame)
    }
}
