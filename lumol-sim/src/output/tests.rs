// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

//! Utilities to test the output algorithms

#![cfg(test)]

use tempfile::NamedTempFile;

use std::fs::File;
use std::io::prelude::*;
use std::path::Path;

use super::Output;
use lumol_core::energy::{Harmonic, PairInteraction};
use lumol_core::{System, Molecule, Particle, UnitCell};
use lumol_core::units;

pub fn test_output<F>(function: F, expected: &str)
where
    F: Fn(&Path) -> Box<dyn Output>,
{
    let tempfile = NamedTempFile::new().unwrap();
    let system = testing_system();
    {
        let mut output = function(tempfile.path());
        output.setup(&system);
        output.write(&system);
        output.finish(&system);
    }

    let file = tempfile.reopen().unwrap();
    check_file_content(file, expected);
}

pub fn testing_system() -> System {
    let mut system = System::with_cell(UnitCell::cubic(10.0));
    system.add_molecule(Molecule::new(Particle::with_position("F", [0.0, 0.0, 0.0].into())));
    system.add_molecule(Molecule::new(Particle::with_position("F", [1.3, 0.0, 0.0].into())));

    system.particles_mut().velocity[0] = [0.1, 0.0, 0.0].into();
    system.particles_mut().velocity[1] = [0.0, 0.0, 0.0].into();

    let harmonic = Box::new(Harmonic {
        k: units::from(300.0, "kJ/mol/A^2").unwrap(),
        x0: units::from(1.2, "A").unwrap(),
    });
    system.set_pair_potential(("F", "F"), PairInteraction::new(harmonic, 5.0));
    system.step = 42;
    return system;
}

fn check_file_content(mut file: File, content: &str) {
    let mut buffer = String::new();
    let _ = file.read_to_string(&mut buffer).unwrap();

    for (l1, l2) in buffer.lines().zip(content.lines()) {
        assert_eq!(l1, l2.trim_start());
    }
}
