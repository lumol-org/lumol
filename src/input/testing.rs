// Cymbalum, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

//! Utilities for testing the input module

use std::path::{Path, PathBuf};
use std::io::ErrorKind;
use std::fs;

pub fn bad_inputs(input: &str, motif: &str) -> Vec<PathBuf> {
    let data_root = Path::new(file!()).parent().unwrap()
                                      .join(input)
                                      .join("data")
                                      .join("bad")
                                      .join(motif);
    return fs::read_dir(data_root).unwrap()
                                  .filter_map(|entry| entry.ok())
                                  .map(|entry| entry.path())
                                  .collect();
}

pub fn cleanup(path: &Path) {
    match path.to_str().unwrap() {
        "src/input/simulations/data/md.toml" => {
            for file in &["energy.dat", "filename.xyz"] {
                if let Err(err) = fs::remove_file(file) {
                    match err.kind() {
                        ErrorKind::NotFound => {},
                        _ => panic!("IO error: {}", err)
                    }
                }
            }
        }
        _ => unreachable!()
    }
}
