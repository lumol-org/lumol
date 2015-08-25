/*
 * Cymbalum, Molecular Simulation in Rust
 * Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
*/

use ::constants::K_BOLTZMANN;
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

/// Compute the potential energy of the system
pub struct PotentialEnergy;
impl Compute for PotentialEnergy {
    type Output = f64;
    fn compute(&self, universe: &Universe) -> f64 {
        let mut res = 0.0;
        for i in 0..universe.size() {
            for j in (i+1)..universe.size() {
                for potential in universe.pairs(i, j) {
                    let d = universe.wrap_vector(i, j);
                    res += potential.energy(d.norm());
                }
            }
        }
        return res;
    }
}

/// Compute the kinetic energy of the system
pub struct KineticEnergy;
impl Compute for KineticEnergy {
    type Output = f64;
    fn compute(&self, universe: &Universe) -> f64 {
        let mut res = 0.0;
        for particle in universe.iter() {
            res += 0.5 * particle.mass() * particle.velocity().norm2();
        }
        return res;
    }
}

/// Compute the total energy of the system
pub struct TotalEnergy;
impl Compute for TotalEnergy {
    type Output = f64;
    fn compute(&self, universe: &Universe) -> f64 {
        let kinetic = KineticEnergy.compute(universe);
        let potential = PotentialEnergy.compute(universe);
        return kinetic + potential;
    }
}

/// Compute the instananeous temperature of the system
pub struct Temperature;
impl Compute for Temperature {
    type Output = f64;
    fn compute(&self, universe: &Universe) -> f64 {
        let kinetic = KineticEnergy.compute(universe);
        let natoms = universe.size() as f64;
        return 1.0/K_BOLTZMANN * 2.0 * kinetic/(3.0 * natoms);
    }
}


#[cfg(test)]
mod test {
    use super::*;
    use ::types::Vector3D;
    use ::universe::{Universe, Particle, UnitCell};
    use ::universe::{InitVelocities, BoltzmanVelocities};
    use ::potentials::Harmonic;
    use ::units;

    const EPS : f64 = 1e-8;

    fn testing_universe() -> Universe {
        let mut universe = Universe::new();
        universe.add_particle(Particle::new("F"));
        universe[0].set_position(Vector3D::new(0.0, 0.0, 0.0));

        universe.add_particle(Particle::new("F"));
        universe[1].set_position(Vector3D::new(1.3, 0.0, 0.0));

        let mut velocities = BoltzmanVelocities::new(units::from(300.0, "K").unwrap());
        velocities.init(&mut universe);

        universe.add_pair_interaction("F", "F",
            Harmonic{k: units::from(300.0, "kJ/mol/A^2").unwrap(), r0: units::from(1.2, "A").unwrap()});
        return universe;
    }

    #[test]
    fn forces() {
        let universe = &testing_universe();
        let res = Forces.compute(universe);

        let mut forces_tot = Vector3D::new(0.0, 0.0, 0.0);
        forces_tot.x += res[0].x + res[1].x;
        forces_tot.x += res[0].y + res[1].y;
        forces_tot.x += res[0].z + res[1].z;

        assert_eq!(forces_tot, Vector3D::new(0.0, 0.0, 0.0));

        assert_approx_eq!(res[0].x, 3e-3, EPS);
        assert_approx_eq!(res[0].y, 0.0, EPS);
        assert_approx_eq!(res[0].y, 0.0, EPS);

        assert_approx_eq!(res[1].x, -3e-3, EPS);
        assert_approx_eq!(res[1].y, 0.0, EPS);
        assert_approx_eq!(res[1].y, 0.0, EPS);
    }

    #[test]
    fn energy() {
        let universe = &testing_universe();
        let kinetic = KineticEnergy.compute(universe);
        let potential = PotentialEnergy.compute(universe);
        let total = TotalEnergy.compute(universe);

        assert_eq!(kinetic + potential, total);
        assert_eq!(kinetic, 0.0007483016557453699);
        assert_approx_eq!(potential, 1.5e-4, EPS);

        assert_eq!(kinetic, universe.kinetic_energy());
        assert_eq!(potential, universe.potential_energy());
        assert_eq!(total, universe.total_energy());
    }

    #[test]
    fn temperature() {
        let universe = &testing_universe();
        let T = Temperature.compute(universe);
        assert_approx_eq!(T, 300.0, 1e-9);
        assert_eq!(T, universe.temperature());
    }
}
