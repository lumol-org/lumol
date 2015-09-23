/*
 * Cymbalum, Molecular Simulation in Rust
 * Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
*/

//! [Chemharp](https://github.com/Luthaf/Chemharp/) library binding and
//! adaptators for Cymbalum.

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

        {
            let topo = universe.topology_mut();
            for bond in try!(topology.bonds()).iter() {
                topo.add_bond(bond[0] as usize, bond[1] as usize);
            }
        }
        Ok(universe)
    }
}

pub fn frame_to_universe(frame: chemharp::Frame) -> Universe {
    match frame.to_cymbalum() {
        Ok(val) => val,
        Err(err) => {
            error!("Error in Chemharp runtime: {}", err.message());
            panic!();
        }
    }
}
