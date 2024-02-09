// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

//! Checking that the documentation tutorials run
use lumol::input::Input;

use std::path::Path;
use std::sync::Once;
static START: Once = Once::new();

struct Cleaner {
    files: Vec<&'static str>
}

impl Cleaner {
    fn new(files: Vec<&'static str>) -> Cleaner {
        Cleaner {
            files
        }
    }
}

impl Drop for Cleaner {
    fn drop(&mut self) {
        for file in &self.files {
            let _ = std::fs::remove_file(file);
        }
    }
}

#[test]
fn argon() {
    START.call_once(::env_logger::init);
    let path = Path::new(file!()).parent()
                                 .unwrap()
                                 .join("..")
                                 .join("doc")
                                 .join("src")
                                 .join("data")
                                 .join("argon.toml");
    let mut config = Input::new(path).unwrap().read().unwrap();

    let _ = Cleaner::new(vec!["energy.dat", "trajectory.xyz"]);

    config.simulation.run(&mut config.system, 1);
}

#[test]
fn nacl() {
    START.call_once(::env_logger::init);
    let path = Path::new(file!()).parent()
                                 .unwrap()
                                 .join("..")
                                 .join("doc")
                                 .join("src")
                                 .join("data")
                                 .join("nacl.toml");
    let mut config = Input::new(path).unwrap().read().unwrap();

    let _ = Cleaner::new(vec!["trajectory.xyz"]);

    config.simulation.run(&mut config.system, 1);
}

#[test]
fn water() {
    START.call_once(::env_logger::init);
    let path = Path::new(file!()).parent()
                                 .unwrap()
                                 .join("..")
                                 .join("doc")
                                 .join("src")
                                 .join("data")
                                 .join("water.toml");
    let mut config = Input::new(path).unwrap().read().unwrap();

    let _ = Cleaner::new(vec!["trajectory.xyz"]);

    config.simulation.run(&mut config.system, 1);
}
