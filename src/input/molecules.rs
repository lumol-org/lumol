// Cymbalum, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

extern crate chemfiles;
use self::chemfiles::{Trajectory, Frame};

use std::path::Path;

use system::chemfiles::Error;
use system::chemfiles::frame_to_system;
use system::{Molecule, Particle};

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
    let system = try!(frame_to_system(frame));

    let molecule = system.molecule(0).clone();
    let mut particles = Vec::new();
    for i in &molecule {
        particles.push(system[i].clone());
    }

    return Ok((molecule, particles));
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;
    use system::moltype;

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
