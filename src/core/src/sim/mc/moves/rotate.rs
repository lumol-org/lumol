// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license
use rand::distributions::{Normal, Range, Sample};
use rand::Rng;

use std::usize;
use std::f64;

use super::MCMove;
use super::select_molecule;

use types::{Matrix3, Vector3D, Zero};
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
        self.newpos.clear();

        let molecule = system.molecule(self.molid);

        let mut masses = vec![0.0; molecule.size()];
        for (i, pi) in molecule.iter().enumerate() {
            masses[i] = system[pi].mass;
            self.newpos.push(system[pi].position);
        }

        rotate_around_axis(&mut self.newpos, &masses, axis, theta);
        return true;
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
            self.theta *= s;
            self.range = Range::new(-self.theta, self.theta);
        }
    }
}

/// Rotate the particles at `positions` with the masses in `masses` around the
/// `axis` axis by `angle`. The `positions` array is overwritten with the new
/// positions.
fn rotate_around_axis(positions: &mut [Vector3D], masses: &[f64], axis: Vector3D, angle: f64) {
    debug_assert!(positions.len() == masses.len());
    // Get center of mass (com) of the molecule
    let total_mass = masses.iter().fold(0.0, |total, m| total + m);
    let com = positions.iter().zip(masses).fold(Vector3D::zero(),
        |com, (&position, &mass)| com + position * mass / total_mass
    );

    let rotation = rotation_matrix(&axis, angle);
    for position in positions {
        let oldpos = *position - com;
        *position = com + rotation * oldpos;
    }
}

fn rotation_matrix(axis: &Vector3D, angle: f64) -> Matrix3 {
    let sin = f64::sin(angle);
    let cos = f64::cos(angle);

    let x_sin = axis[0] * sin;
    let y_sin = axis[1] * sin;
    let z_sin = axis[2] * sin;
    let one_minus_cos = 1.0 - cos;
    let xy = axis[0] * axis[1] * one_minus_cos;
    let xz = axis[0] * axis[2] * one_minus_cos;
    let yz = axis[1] * axis[2] * one_minus_cos;

    // Build the rotation matrix
    let rotation = Matrix3::new(
        (axis[0] * axis[0]) * one_minus_cos + cos,
        xy + z_sin,
        xz - y_sin,
        xy - z_sin,
        (axis[1] * axis[1]) * one_minus_cos + cos,
        yz + x_sin,
        xz + y_sin,
        yz - x_sin,
        (axis[2] * axis[2]) * one_minus_cos + cos
    );

    return rotation;
}
