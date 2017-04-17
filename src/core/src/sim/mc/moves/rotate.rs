// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license
use rand::distributions::{Normal, Range, Sample};
use rand::Rng;

use std::usize;
use std::f64;

use super::MCMove;
use super::select_molecule;

use types::{Matrix3, Vector3D};
use sys::{System, EnergyCache};

/// Monte-Carlo move for rotating a rigid molecule
pub struct Rotate {
    /// Type of molecule to rotate. `None` means all molecules.
    moltype: Option<u64>,
    /// Index of the molecule to rotate
    molid: usize,
    /// New positions of the atom in the rotated molecule
    newpos: Vec<Vector3D>,
    /// Normal distribution, for generation of the axis
    axis_rng: Normal,
    /// Maximum values for the range of the range distribution of the angle
    theta: f64,
    /// Range distribution, for generation of the angle
    range: Range<f64>,
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
            newpos: Vec::new(),
            axis_rng: Normal::new(0.0, 1.0),
            theta: theta,
            range: Range::new(-theta, theta),
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

    fn setup(&mut self, _: &System) { }

    fn prepare(&mut self, system: &mut System, rng: &mut Box<Rng>) -> bool {
        if let Some(id) = select_molecule(system, self.moltype, rng) {
            self.molid = id;
        } else {
            warn!("Can not rotate molecule: no molecule of this type in the system.");
            return false;
        }

        // Getting values from a 3D normal distribution gives an uniform
        // distribution on the unit sphere.
        let axis = Vector3D::new(
            self.axis_rng.sample(rng),
            self.axis_rng.sample(rng),
            self.axis_rng.sample(rng)
        ).normalized();
        let theta = self.range.sample(rng);

        // get indices of particles of selected molecule
        let molecule = system.molecule(self.molid);
        // store positions of selected molecule
        self.newpos.clear();
        for pi in molecule.iter() {
            self.newpos.push(system[pi].position);
        }
        // get center-of-mass of molecule
        let com = system.molecule_com(self.molid);
        rotate_around_axis(&mut self.newpos, com, axis, theta);
        true
    }

    fn cost(&self, system: &System, beta: f64, cache: &mut EnergyCache) -> f64 {
        let idxes = system.molecule(self.molid).iter().collect::<Vec<_>>();
        let cost = cache.move_particles_cost(system, idxes, &self.newpos);
        return cost*beta;
    }

    fn apply(&mut self, system: &mut System) {
        for (i, pi) in system.molecule(self.molid).iter().enumerate() {
            system[pi].position = self.newpos[i];
        }
    }

    fn restore(&mut self, _: &mut System) {
        // Nothing to do
    }

    fn update_amplitude(&mut self, scaling_factor: Option<f64>) {
        if let Some(s) = scaling_factor {
            if (s * self.theta).abs().to_degrees() <= 180.0 {
                self.theta *= s;
                self.range = Range::new(-self.theta, self.theta);
            } else {
                warn_once!(
                    "Tried to increase the maximum amplitude for rotations to more than 180°."
                );
            }
        }
    }
}

/// Rotate the particles at `positions` with the center-of-mass position
/// `com` around the `axis` axis by `angle`. The `positions` array is
/// overwritten with the new positions.
fn rotate_around_axis(positions: &mut [Vector3D], com: Vector3D, axis: Vector3D, angle: f64) {
    let rotation = Matrix3::rotation(&axis, angle);
    for position in positions {
        let oldpos = *position - com;
        *position = com + rotation * oldpos;
    }
}
