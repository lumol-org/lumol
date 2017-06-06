// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

//! Utilities to test the output algorithms

#![cfg(test)]
extern crate tempfile;
use self::tempfile::NamedTempFile;

use std::io::prelude::*;
use std::fs::File;
use std::path::Path;

use super::Output;
use sys::System;
use energy::{PairInteraction, Harmonic};
use utils::{unit_from, system_from_xyz};

pub fn test_output<F>(function: F, expected: &str) where F: Fn(&Path) -> Box<Output> {
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

fn testing_system() -> System {
    let mut system = system_from_xyz("2
    cell: 10
    F 0.0 0.0 0.0 0.1 0.0 0.0
    F 1.3 0.0 0.0 0.0 0.0 0.0
    ");

    let harmonic = Box::new(Harmonic{
        k: unit_from(300.0, "kJ/mol/A^2"),
        x0: unit_from(1.2, "A")
    });
    system.add_pair_potential("F", "F",
        PairInteraction::new(harmonic, 5.0)
    );
    return system;
}

fn check_file_content(mut file: File, content: &str) {
    let mut buffer = String::new();
    let _ = file.read_to_string(&mut buffer).unwrap();

    assert_eq!(buffer, content);
}
