// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

//! Read static string using the XYZ file format, and create the corresponding
//! system.

use system::{System, Particle};
use types::Vector3D;
use input::guess_bonds;

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
        return guess_bonds(system).expect("Could not guess the bonds");
    }

    return system;
}

#[cfg(test)]
mod tests {
    use super::*;
    use types::{Vector3D, Zero};

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

        assert_eq!(system[0].position, Vector3D::new(0.0, 0.0, -1.5));

        assert_eq!(system[1].position, Vector3D::zero());

        assert_eq!(system[2].position, Vector3D::new(0.0, 0.0, 1.5));

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
