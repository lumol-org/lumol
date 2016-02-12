/* Cymbalum, Molecular Simulation in Rust - Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
 */

extern crate chemfiles;
use self::chemfiles::{Trajectory, Frame};

use std::path::Path;

use universe::chemfiles::Error;
use universe::chemfiles::frame_to_universe;
use universe::{Molecule, Particle};

/// Read the first molecule in the first frame of the file at `path`, and return
/// this molecule and the corresponding list of particles. If `guess_bonds` is
/// `true`, then the bonds are automatically guessed. This functionality uses
/// [chemfiles](https://crates.io/crates/chemfiles), read the corresponding
/// documentation for more informations.
pub fn molecule_from_file<P: AsRef<Path>>(path: P, guess_bonds: bool) -> Result<(Molecule, Vec<Particle>), Error> {
    let mut trajectory = try!(Trajectory::open(path));
    let mut frame = try!(Frame::new(0));
    try!(trajectory.read(&mut frame));
    try!(frame.guess_topology(guess_bonds));
    let universe = try!(frame_to_universe(frame));

    let molecule = universe.molecule(0).clone();
    let mut particles = Vec::new();
    for i in &molecule {
        particles.push(universe[i].clone());
    }

    return Ok((molecule, particles));
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;
    use universe::moltype;

    #[test]
    fn read() {
        let file = Path::new(file!()).parent().unwrap().join("data").join("H20.xyz");
        let (molecule, atoms) = molecule_from_file(file, true).unwrap();

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
