/*
 * Cymbalum, Molecular Simulation in Rust
 * Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
*/

use ::types::Vector3D;
use ::universe::Universe;

/// The compute trait allow to compute properties of an universe, whithout
/// modifying this universe. The Output type is the type of the computed
/// property.
pub trait Compute {
    type Output;
    /// Compute the property
    fn compute(&self, universe: &Universe) -> Self::Output;
}

/// Compute all the forces acting on the system, and return a vector of
/// force acting on each particles
pub struct Forces;
impl Compute for Forces {
    type Output = Vec<Vector3D>;
    fn compute(&self, universe: &Universe) -> Vec<Vector3D> {
        let mut res: Vec<Vector3D> = Vec::with_capacity(universe.size());
        for _ in 0..universe.size() {
            res.push(Vector3D::new(0.0, 0.0, 0.0));
        }

        for i in 0..universe.size() {
            for j in (i+1)..universe.size() {
                for potential in universe.pairs(i, j) {
                    let d = universe.wrap_vector(i, j);
                    let dn = d.normalized();
                    let f = potential.force(d.norm());
                    res[i] = res[i] + f * dn;
                    res[j] = res[j] - f * dn;
                }
            }
        }
        return res;
    }
}



#[cfg(test)]
mod test {
    use super::*;
    use ::types::Vector3D;
    use ::universe::{Universe, Particle};
    use ::potentials::Harmonic;

    #[test]
    fn forces() {
        let mut universe = Universe::new();
        universe.add_particle(Particle::new("F"));
        universe[0].set_position(Vector3D::new(0.0, 0.0, 0.0));
        universe.add_particle(Particle::new("F"));
        universe[1].set_position(Vector3D::new(1.3, 0.0, 0.0));

        universe.add_pair_interaction("F", "F", Harmonic{k: 300.0, r0: 1.2});
        let res = Forces.compute(&universe);

        let mut forces_tot = Vector3D::new(0.0, 0.0, 0.0);
        forces_tot.x += res[0].x + res[1].x;
        forces_tot.x += res[0].y + res[1].y;
        forces_tot.x += res[0].z + res[1].z;

        assert_eq!(forces_tot, Vector3D::new(0.0, 0.0, 0.0));

        assert_approx_eq!(res[0].x, 30.0, 1e-12);
        assert_approx_eq!(res[0].y, 0.0, 1e-12);
        assert_approx_eq!(res[0].y, 0.0, 1e-12);

        assert_approx_eq!(res[1].x, -30.0, 1e-12);
        assert_approx_eq!(res[1].y, 0.0, 1e-12);
        assert_approx_eq!(res[1].y, 0.0, 1e-12);
    }
}
