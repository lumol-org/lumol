/* Cymbalum, Molecular Simulation in Rust - Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
 */
extern crate rand;
use self::rand::distributions::{Normal, Range, Sample};
use self::rand::Rng;

use std::usize;
use std::f64;

use super::MCMove;
use super::select_molecule;

use types::Vector3D;
use universe::Universe;

/// Monte-Carlo move for rotating a rigid molecule
pub struct Rotate {
    /// Type of molecule to rotate. `None` means all molecules.
    moltype: Option<u64>,
    /// Index of the molecule to rotate
    molid: usize,
    /// Previous positions for all the particles in the molecule
    prevpos: Vec<Vector3D>,
    /// Normal distribution, for generation of the axis
    axis_rng: Normal,
    /// Range distribution, for generation of the angle
    angle_rng: Range<f64>,
    /// Potential energy before the move
    e_before: f64,
}

impl Rotate {
    /// Create a new `Rotate` move, with maximum angular displacement of `theta`,
    /// rotating all the molecules in the system.
    pub fn new(theta: f64) -> Rotate {
        Rotate::create(theta, None)
    }

    /// Create a new `Rotate` move, with maximum angular displacement of `theta`,
    /// rotating only molecules with `moltype` type.
    pub fn with_moltype(theta: f64, moltype: u64) -> Rotate {
        Rotate::create(theta, Some(moltype))
    }

    // Factorizing the constructors
    fn create(theta: f64, moltype: Option<u64>) -> Rotate {
        assert!(theta > 0.0, "theta must be positive in Rotate move");
        Rotate {
            moltype: moltype,
            molid: usize::MAX,
            prevpos: Vec::new(),
            axis_rng: Normal::new(0.0, 1.0),
            angle_rng: Range::new(-theta, theta),
            e_before: 0.0,
        }
    }
}

impl Default for Rotate {
    fn default() -> Rotate {
        Rotate::new(0.2)
    }
}

impl MCMove for Rotate {
    fn describe(&self) -> &str {
        "molecular rotation"
    }

    fn prepare(&mut self, universe: &mut Universe, rng: &mut Box<Rng>) -> bool {
        if let Some(id) = select_molecule(universe, self.moltype, rng) {
            self.molid = id;
        } else {
            warn!("Can not rotate molecule: no molecule of this type in the universe.");
            return false;
        }

        self.e_before = universe.potential_energy();

        // Getting values from a 3D normal distribution gives an uniform
        // distribution on the R3 sphere.
        let axis = Vector3D::new(
            self.axis_rng.sample(rng),
            self.axis_rng.sample(rng),
            self.axis_rng.sample(rng)
        ).normalized();
        let theta = self.angle_rng.sample(rng);

        self.prevpos = vec![Vector3D::new(0.0, 0.0, 0.0); universe.molecule(self.molid).size()];
        for (i, idx) in universe.molecule(self.molid).iter().enumerate() {
            self.prevpos[i] = universe[idx].position;
        }

        rotate_molecule_around_axis(universe, self.molid, axis, theta);
        return true;
    }

    fn cost(&self, universe: &Universe, beta: f64) -> f64 {
        let e_after = universe.potential_energy();
        return (e_after - self.e_before)/beta;
    }

    fn apply(&mut self, _: &mut Universe) {
        // Nothing to do
    }

    fn restore(&mut self, universe: &mut Universe) {
        for (i, idx) in universe.molecule(self.molid).iter().enumerate() {
            universe[idx].position = self.prevpos[i];
        }
    }
}

/// Rotate the molecule at index `molid` in universe around the `axis` axis by
/// `angle`.
fn rotate_molecule_around_axis(universe: &mut Universe, molid: usize, axis: Vector3D, angle: f64) {
    let molecule = universe.molecule(molid).clone();

    // Get center of mass (com) of the molecule
    let total_mass = molecule.iter().fold(0.0, |total_mass, i| total_mass + universe[i].mass);
    let com = molecule.iter().fold(Vector3D::new(0.0, 0.0, 0.0), |com, i| com + universe[i].position / total_mass);

    for i in &molecule {
        let oldpos = universe[i].position.clone() - com;
        let newpos = com + f64::cos(angle) * oldpos
                         + f64::sin(angle) * axis ^ oldpos
                         + (axis * oldpos) * (1.0 - f64::cos(angle)) * axis;
        universe[i].position = newpos;
    }
}
