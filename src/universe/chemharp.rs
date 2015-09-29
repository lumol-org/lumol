/*
 * Cymbalum, Molecular Simulation in Rust
 * Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
*/

//! [Chemharp](https://github.com/Luthaf/Chemharp/) adaptators for Cymbalum.

extern crate chemharp;
use ::universe::{Particle, Universe};
use ::universe::cells::{UnitCell, CellType};
use ::types::Vector3D;

use self::chemharp::Error;

/// Convert chemharp types to Cymbalum types
trait ToCymbalum {
    /// Output type
    type Output;
    /// Conversion function
    fn to_cymbalum(self) -> Result<Self::Output, Error>;
}

impl ToCymbalum for chemharp::Atom {
    type Output = Particle;
    fn to_cymbalum(self) -> Result<Particle, Error> {
        let name = try!(self.name());
        let mut part = Particle::new(name);
        let mass = try!(self.mass());
        part.set_mass(mass as f64);
        Ok(part)
    }
}

impl ToCymbalum for chemharp::UnitCell {
    type Output = UnitCell;
    fn to_cymbalum(self) -> Result<UnitCell, Error> {
        let cell_type = try!(self.cell_type());
        let cell = match cell_type {
            chemharp::CellType::Infinite => UnitCell::new(),
            chemharp::CellType::Orthorombic => {
                let (a, b, c) = try!(self.lengths());
                UnitCell::ortho(a, b, c)
            },
            chemharp::CellType::Triclinic => {
                let (a, b, c) = try!(self.lengths());
                let (alpha, beta, gamma) = try!(self.angles());
                UnitCell::triclinic(a, b, c, alpha, beta, gamma)
            }
        };
        Ok(cell)
    }
}

impl ToCymbalum for chemharp::Frame {
    type Output = Universe;
    fn to_cymbalum(self) -> Result<Universe, Error> {
        let cell = try!(self.cell());
        let cell = try!(cell.to_cymbalum());
        let mut universe = Universe::from_cell(cell);
        let topology = try!(self.topology());
        let positions = try!(self.positions());
        let natoms = try!(self.natoms());
        let step = try!(self.step());

        universe.set_step(step);

        for i in 0..natoms {
            let atom = try!(topology.atom(i as u64));
            let particle = try!(atom.to_cymbalum());

            universe.add_particle(particle);
            let position = Vector3D::new(
                positions[i][0] as f64,
                positions[i][1] as f64,
                positions[i][2] as f64
            );
            universe[i].set_position(position);
        }

        for bond in try!(topology.bonds()).iter() {
            universe.add_bond(bond[0] as usize, bond[1] as usize);
        }
        Ok(universe)
    }
}

/// Convert a chemharp `Frame` to an `Universe`
pub fn frame_to_universe(frame: chemharp::Frame) -> Result<Universe, Error> {
    frame.to_cymbalum()
}


/******************************************************************************/

/// Convert Cymbalum types to chemharp types
trait ToChemharp {
    /// Output type
    type Output;
    /// Conversion function
    fn to_chemharp(&self) -> Result<Self::Output, Error>;
}

impl ToChemharp for Particle {
    type Output = chemharp::Atom;
    fn to_chemharp(&self) -> Result<chemharp::Atom, Error> {
        let mut res = try!(chemharp::Atom::new(self.name()));
        try!(res.set_mass(self.mass() as f32));
        return Ok(res);
    }
}

impl ToChemharp for UnitCell {
    type Output = chemharp::UnitCell;
    fn to_chemharp(&self) -> Result<chemharp::UnitCell, Error> {
        let res = match self.celltype() {
            CellType::Infinite => {
                unimplemented!()
            }
            CellType::Orthorombic => {
                let (a, b, c) = (self.a(), self.b(), self.c());
                try!(chemharp::UnitCell::new(a, b, c))
            },
            CellType::Triclinic => {
                let (a, b, c) = (self.a(), self.b(), self.c());
                let (alpha, beta, gamma) = (self.alpha(), self.beta(), self.gamma());
                try!(chemharp::UnitCell::triclinic(a, b, c, alpha, beta, gamma))
            },
        };
        return Ok(res);
    }
}

impl ToChemharp for Universe {
    type Output = chemharp::Frame;
    fn to_chemharp(&self) -> Result<chemharp::Frame, Error> {

        let natoms = self.size();
        let mut frame = try!(chemharp::Frame::new(natoms as u64));

        try!(frame.set_step(self.step()));

        let mut topology = try!(chemharp::Topology::new());
        let mut positions = vec![[0.0f32; 3]; natoms];
        let mut velocities = vec![[0.0f32; 3]; natoms];

        for (i, p) in self.iter().enumerate() {
            let pos = p.position();
            positions[i][0] = pos.x as f32;
            positions[i][1] = pos.y as f32;
            positions[i][2] = pos.z as f32;

            let vel = p.velocity();
            velocities[i][0] = vel.x as f32;
            velocities[i][1] = vel.y as f32;
            velocities[i][2] = vel.z as f32;

            let atom = try!(p.to_chemharp());
            try!(topology.push(&atom));
        }
        try!(frame.set_positions(positions));
        try!(frame.set_velocities(velocities));

        for bond in self.topology().bonds().iter() {
            try!(topology.add_bond(bond.i() as u64, bond.j() as u64));
        }
        try!(frame.set_topology(&topology));
        // Guessing angles and dihedrals
        try!(frame.guess_topology(false));

        let cell = try!(self.cell().to_chemharp());
        try!(frame.set_cell(&cell));
        Ok(frame)
    }
}

/// Convert an `Universe` to a chemharp `Frame`
pub fn universe_to_frame(universe: &Universe) -> Result<chemharp::Frame, Error> {
    universe.to_chemharp()
}
