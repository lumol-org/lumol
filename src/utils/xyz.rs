/* Cymbalum, Molecular Simulation in Rust - Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
 */
//! Read static string using the XYZ file format, and create the corresponding
//! system.

use system::{System, Particle};
use system::chemfiles::{frame_to_system, system_to_frame};
use types::Vector3D;

/// Read the `content` string, assuming XYZ format, and create the corresponding
/// system. This function is intended for testing purposes only, and will
/// panic if the string is not well-formated.
///
/// If the comment line contains `bonds`, chemfiles will be used to guess the
/// bonds in the system.
pub fn system_from_xyz(content: &str) -> System {
    let mut system = System::new();

    let lines = content.lines().collect::<Vec<_>>();
    let natoms = lines[0].trim().parse::<usize>().expect("Could not parse integer");
    for i in 0..natoms {
        let splitted = lines[i + 2].split_whitespace().collect::<Vec<_>>();
        let name = splitted[0];
        let x = splitted[1].parse::<f64>().expect("Could not parse float");
        let y = splitted[2].parse::<f64>().expect("Could not parse float");
        let z = splitted[3].parse::<f64>().expect("Could not parse float");
        let mut particle = Particle::new(name);
        particle.position = Vector3D::new(x, y, z);
        system.add_particle(particle);
    }

    if lines[1].contains("bonds") {
        let frame = system_to_frame(&system).unwrap();
        frame.guess_topology(true).expect("Could not guess the topology");
        system = frame_to_system(frame).unwrap();
    }

    return system;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bonds() {
        let file = "3
        bonds
        O 0 0 -1.5
        C 0 0 0
        O 0 0 1.5";

        let system = system_from_xyz(file);
        assert_eq!(system.size(), 3);
        assert_eq!(system[0].name(), "O");
        assert_eq!(system[1].name(), "C");
        assert_eq!(system[2].name(), "O");

        assert_eq!(system[0].position.x, 0.0);
        assert_eq!(system[0].position.y, 0.0);
        assert_eq!(system[0].position.z, -1.5);

        assert_eq!(system[1].position.x, 0.0);
        assert_eq!(system[1].position.y, 0.0);
        assert_eq!(system[1].position.z, 0.0);

        assert_eq!(system[2].position.x, 0.0);
        assert_eq!(system[2].position.y, 0.0);
        assert_eq!(system[2].position.z, 1.5);

        assert_eq!(system.molecules().len(), 1);
        assert_eq!(system.molecule(0).bonds().len(), 2);
    }

    #[test]
    fn no_bonds() {
        let file = "4

        He 0 0 0
        He 1 0 0
        He 0 1 0
        He 0 0 1";

        let system = system_from_xyz(file);
        assert_eq!(system.size(), 4);
        assert_eq!(system.molecules().len(), 4);
    }
}
