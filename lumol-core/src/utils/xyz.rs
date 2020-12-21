// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

//! Read static string using the XYZ file format, and create the corresponding
//! system.

use crate::{Molecule, Particle, System, UnitCell};

/// Read the `content` string, assuming XYZ format, and create the corresponding
/// system. This function is intended for testing purposes only, and will
/// panic if the string is not well-formatted.
///
/// If the comment line contains `cell: <a>`, the system will have a cubic unit
/// cell of size a.
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
        let mut particle = Particle::with_position(name, [x, y, z].into());
        if splitted.len() == 7 {
            let vx = splitted[4].parse::<f64>().expect("Could not parse float");
            let vy = splitted[5].parse::<f64>().expect("Could not parse float");
            let vz = splitted[6].parse::<f64>().expect("Could not parse float");
            particle.velocity = [vx, vy, vz].into();
        }
        system.add_molecule(Molecule::new(particle));
    }

    if lines[1].contains("cell:") {
        let cell = lines[1].split("cell:").nth(1).expect("Missing cell size");
        let cell = cell.split_whitespace().next().expect("Missing cell size");
        let cell = UnitCell::cubic(cell.parse().expect("Could not parse float"));
        system.cell = cell;
    }

    return system;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Vector3D;

    #[test]
    fn cell() {
        let system = system_from_xyz(
            "4
            cell: 67
            He 0 0 0
            He 1 0 0
            He 0 1 0
            He 0 0 1",
        );
        assert_eq!(system.size(), 4);
        assert_eq!(system.molecules().count(), 4);
        assert_eq!(system.cell, UnitCell::cubic(67.0));

        assert_eq!(system.particles().position[0], Vector3D::new(0.0, 0.0, 0.0));
        assert_eq!(system.particles().position[1], Vector3D::new(1.0, 0.0, 0.0));
        assert_eq!(system.particles().position[2], Vector3D::new(0.0, 1.0, 0.0));
        assert_eq!(system.particles().position[3], Vector3D::new(0.0, 0.0, 1.0));
    }

    #[test]
    fn velocities() {
        let system = system_from_xyz(
            "4
            cell: 67
            He 0 0 0 0 0 0
            He 1 0 0 1 2 3
            He 0 1 0 0 1 0
            He 0 0 1 2 2 3
            ",
        );
        assert_eq!(system.size(), 4);
        assert_eq!(system.molecules().count(), 4);
        assert_eq!(system.cell, UnitCell::cubic(67.0));

        assert_eq!(system.particles().position[0], Vector3D::new(0.0, 0.0, 0.0));
        assert_eq!(system.particles().position[1], Vector3D::new(1.0, 0.0, 0.0));
        assert_eq!(system.particles().position[2], Vector3D::new(0.0, 1.0, 0.0));
        assert_eq!(system.particles().position[3], Vector3D::new(0.0, 0.0, 1.0));

        assert_eq!(system.particles().velocity[0], Vector3D::new(0.0, 0.0, 0.0));
        assert_eq!(system.particles().velocity[1], Vector3D::new(1.0, 2.0, 3.0));
        assert_eq!(system.particles().velocity[2], Vector3D::new(0.0, 1.0, 0.0));
        assert_eq!(system.particles().velocity[3], Vector3D::new(2.0, 2.0, 3.0));
    }
}
