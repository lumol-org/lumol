// Cymbalum, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

//! Utilities for testing the input module

use std::path::{Path, PathBuf};
use std::fs;

pub fn bad_interactions(motif: &str) -> Vec<PathBuf> {
    let data_root = Path::new(file!()).parent().unwrap()
                                      .join("data")
                                      .join("bad")
                                      .join(motif);
    return fs::read_dir(data_root).unwrap()
                                  .filter_map(|entry| entry.ok())
                                  .map(|entry| entry.path())
                                  .collect();
}
