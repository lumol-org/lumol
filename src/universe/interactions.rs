/*
 * Cymbalum, Molecular Simulation in Rust
 * Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
*/
use std::collections::HashMap;

use ::potentials::{PairPotential, AnglePotential, DihedralPotential};

/// The Interaction type hold all data about the potentials in the system,
/// indexed by particle type.
pub struct Interactions {
    /// Pair potentials
    pub pairs: HashMap<(u16, u16), Vec<Box<PairPotential>>>,
    /// Angle potentials
    pub angles: HashMap<(u16, u16, u16), Vec<Box<AnglePotential>>>,
    /// Dihedral angles potentials
    pub dihedrals: HashMap<(u16, u16, u16, u16), Vec<Box<DihedralPotential>>>,
}

impl Interactions {
    pub fn new() -> Interactions {
        Interactions{
            pairs: HashMap::new(),
            angles: HashMap::new(),
            dihedrals: HashMap::new(),
        }
    }
}
