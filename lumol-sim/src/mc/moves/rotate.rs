// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license
use rand::RngCore;
use rand_distr::{Distribution, Uniform, UnitSphere};

use std::collections::BTreeSet;
use std::f64;
use std::usize;

use log::warn;
use log_once::warn_once;

use soa_derive::soa_zip;

use super::{MCDegreeOfFreedom, MCMove};
use super::select_molecule;

use lumol_core::{EnergyCache, System, MoleculeHash, Matrix3, Vector3D};

/// Monte Carlo move for rotating a rigid molecule
pub struct Rotate {
    /// Hash of molecule to rotate. `None` means all molecules.
    hash: Option<MoleculeHash>,
    /// Index of the molecule to rotate
    molid: usize,
    /// New positions of the atom in the rotated molecule
    newpos: Vec<Vector3D>,
    /// Maximum values for the range of the range distribution of the angle
    theta: f64,
    /// Range distribution, for generation of the angle
    range: Uniform<f64>,
}

impl Rotate {
    /// Create a new `Rotate` move, with maximum angular displacement of
    /// `theta`. This move will apply to the molecules with the given `hash`,
    /// or all molecules if `hash` is `None`.
    pub fn new<H: Into<Option<MoleculeHash>>>(theta: f64, hash: H) -> Rotate {
        assert!(theta > 0.0, "theta must be positive in Rotate move");
        Rotate {
            hash: hash.into(),
            molid: usize::max_value(),
            newpos: Vec::new(),
            theta: theta,
            range: Uniform::new(-theta, theta),
        }
    }
}

impl MCMove for Rotate {
    fn describe(&self) -> &str {
        "molecular rotation"
    }

    fn degrees_of_freedom(&self) -> MCDegreeOfFreedom {
        match self.hash {
            Some(hash) => {
                let mut all = BTreeSet::new();
                let _ = all.insert(hash);
                MCDegreeOfFreedom::Molecules(all)
            }
            None => MCDegreeOfFreedom::AllMolecules,
        }
    }

    fn setup(&mut self, _: &System) {}

    fn prepare(&mut self, system: &mut System, rng: &mut dyn RngCore) -> bool {
        if let Some(id) = select_molecule(system, self.hash, rng) {
            self.molid = id;
        } else {
            warn!("Can not rotate molecule: no molecule of this type in the system.");
            return false;
        }

        let axis = Vector3D::from(UnitSphere.sample(rng));
        let theta = self.range.sample(rng);

        // store positions of selected molecule
        self.newpos = system.molecule(self.molid).particles().position.to_vec();
        // get center-of-mass of molecule
        let com = system.molecule(self.molid).center_of_mass();
        rotate_around_axis(&mut self.newpos, com, axis, theta);
        true
    }

    fn cost(&self, system: &System, beta: f64, cache: &mut EnergyCache) -> f64 {
        return beta * cache.move_molecule_cost(system, self.molid, &self.newpos);
    }

    fn apply(&mut self, system: &mut System) {
        let mut molecule = system.molecule_mut(self.molid);
        for (position, newpos) in soa_zip!(molecule.particles_mut(), [mut position], &self.newpos) {
            *position = *newpos;
        }
    }

    fn restore(&mut self, _: &mut System) {
        // Nothing to do
    }

    fn update_amplitude(&mut self, scaling_factor: Option<f64>) {
        if let Some(s) = scaling_factor {
            if (s * self.theta).abs().to_degrees() <= 180.0 {
                self.theta *= s;
                self.range = Uniform::new(-self.theta, self.theta);
            } else {
                warn_once!(
                    "Tried to increase the maximum amplitude for rotations to more than 180 degrees."
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
