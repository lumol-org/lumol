/*
 * Cymbalum, Molecular Simulation in Rust
 * Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
*/

use ::types::Vector3D;
use ::types::Matrix3;

pub trait PairPotential {
    fn energy(&self, r: f64) -> f64;
    fn force(&self, r: f64) -> f64;

    fn virial(&self, r: &Vector3D) -> Matrix3 {
        let fact = self.force(r.norm());
        let rn = r.normalize();
        let force = fact * rn;
        force.tensorial(r)
    }
}
