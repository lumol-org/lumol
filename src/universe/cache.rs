/* Cymbalum, Molecular Simulation in Rust - Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
 */
//! Caching energy components for speeding up the energy computation in
//! Monte-Carlo simulations.
//!
//! In most of Monte-Carlo moves, only a very small subset of the system changes.
//! We can use that property to remove the need of recomputing most of the
//! energy components, by storing them and providing update callbacks.
use std::mem;

use super::Universe;
use types::{Vector3D, Array2};
use potentials::global::GlobalCache;

/// Callback for updating a cache. It also take an `&mut Universe` argument for
/// updating the cache inside the global potentials.
type UpdateCallback = Box<Fn(&mut EnergyCache, &mut Universe)>;

/// This is a cache for energy computation.
///
/// Cache integrity is left up to the user of this structure: any function with
/// `update_` prefix must be called as needed to ensure cache consistency.
pub struct EnergyCache {
    /// 2-D array containing the pairs interactions between particles `i` and
    /// `j` at index `i, j`
    pairs_cache: Array2<f64>,
    /// Energy of all the pairs in the system
    pairs: f64,
    /// Energy of all the bonds in the system
    bonds: f64,
    /// Energy of all the angles in the system
    angles: f64,
    /// Energy of all the dihedrals angles in the system
    dihedrals: f64,
    /// Energy of coulombic interactions
    coulomb: f64,
    /// Energy of global interactions
    global: f64,
    /// Callback to be called to update the cache if the universe is modified
    updater: Option<UpdateCallback>
}

impl EnergyCache {
    /// Create a new empty energy cache.
    pub fn new() -> EnergyCache {
        EnergyCache {
            pairs_cache: Array2::new(),
            pairs: 0.0,
            bonds: 0.0,
            angles: 0.0,
            dihedrals: 0.0,
            coulomb: 0.0,
            global: 0.0,
            updater: None,
        }
    }

    /// Clear all values in the cache by setting them to 0
    fn clear(&mut self) {
        self.pairs_cache.fill(0.0);
        self.pairs = 0.0;
        self.bonds = 0.0;
        self.angles = 0.0;
        self.dihedrals = 0.0;
        self.coulomb = 0.0;
        self.global = 0.0;
    }

    /// Initialize the cache to be used with `universe`. After a call to this
    /// function, the cache is only usable with the same universe. To change
    /// the associated universe, one must call this fucntion again.
    pub fn init(&mut self, universe: &Universe) {
        self.clear();
        self.pairs_cache.resize((universe.size(), universe.size()));
        let evaluator = universe.energy_evaluator();

        for i in 0..universe.size() {
            for j in (i + 1)..universe.size() {
                let r = universe.wraped_vector(i, j).norm();
                let energy = evaluator.pair(r, i, j);
                self.pairs_cache[(i, j)] = energy;
                self.pairs_cache[(j, i)] = energy;
                self.pairs += energy;
            }
        }

        self.bonds = evaluator.bonds();
        self.angles = evaluator.angles();
        self.dihedrals = evaluator.dihedrals();
        self.coulomb = evaluator.coulomb();
        self.global = evaluator.global();
    }

    /// Get the cached energy
    pub fn energy(&self) -> f64 {
        let mut energy = 0.0;
        energy += self.pairs;

        energy += self.bonds;
        energy += self.angles;
        energy += self.dihedrals;

        energy += self.coulomb;
        energy += self.global;
        return energy;
    }

    /// Update the cache after a call to a `EnergyCache::*_cost` function or
    /// `EnergyCache::unused`.
    pub fn update(&mut self, universe: &mut Universe) {
        let updater = mem::replace(&mut self.updater, None);
        if let Some(updater) = updater {
            updater(self, universe);
        } else {
            error!("Error in `EnergyCache::update`: \
                    This function MUST be called after a call to a `*_cost` function");
            panic!()
        }
    }

    /// This function should be called whenever the cache is not used, but one
    /// still want it to be updated. Future call to `EnergyCache::update` will
    /// recompute the full cache.
    pub fn unused(&mut self) {
        self.updater = Some(Box::new(|cache, universe| {
            cache.init(universe);
        }))
    }
}

impl EnergyCache {
    /// Get the cost of moving the set of particles with indexes in `idxes` to
    /// `newpos`. This function DO NOT update the cache, the
    /// `update_particles_moved` function MUST be called if the particles are
    /// effectively moved.
    pub fn move_particles_cost(&mut self, universe: &Universe, idxes: Vec<usize>, newpos: &[Vector3D]) -> f64 {
        let evaluator = universe.energy_evaluator();

        // First, go for pair interactions
        let mut new_pairs = Array2::<f64>::with_size((universe.size(), universe.size()));
        let mut pairs_delta = 0.0;
        // Interactions with the sub-system not being moved
        for (i, &part_i) in idxes.iter().enumerate() {
            for part_j in 0..universe.size() {
                // Exclude interactions inside the sub-system.
                if idxes.contains(&part_j) {continue}

                let r = universe.cell().distance(&universe[part_j].position, &newpos[i]);
                let energy = evaluator.pair(r, part_i, part_j);

                pairs_delta += energy;
                new_pairs[(part_i, part_j)] += energy;
                new_pairs[(part_j, part_i)] += energy;

                pairs_delta -= self.pairs_cache[(part_i, part_j)];
            }
        }

        // Interactions withing the sub-system being moved
        for (i, &part_i) in idxes.iter().enumerate() {
            for (j, &part_j) in idxes.iter().enumerate().skip(i + 1) {
                let r = universe.cell().distance(&newpos[i], &newpos[j]);
                let energy = evaluator.pair(r, part_i, part_j);

                pairs_delta += energy;
                new_pairs[(part_i, part_j)] += energy;
                new_pairs[(part_j, part_i)] += energy;

                pairs_delta -= self.pairs_cache[(part_i, part_j)];
            }
        }


        // Recompute everything here. One day we could store some terms and do
        // the same thing than for pair interactions.
        let mut bonds = 0.0;
        let mut angles = 0.0;
        let mut dihedrals = 0.0;
        for molecule in universe.molecules() {
            for bond in molecule.bonds() {
                let (i, j) = (bond.i(), bond.j());
                let ri = new_position(universe, i, &idxes, newpos);
                let rj = new_position(universe, j, &idxes, newpos);
                let r = universe.cell().distance(ri, rj);
                bonds += evaluator.bond(r, i, j);
            }

            for angle in molecule.angles() {
                let (i, j, k) = (angle.i(), angle.j(), angle.k());
                let ri = new_position(universe, i, &idxes, newpos);
                let rj = new_position(universe, j, &idxes, newpos);
                let rk = new_position(universe, k, &idxes, newpos);
                let theta = universe.cell().angle(ri, rj, rk);
                angles += evaluator.angle(theta, i, j, k);
            }

            for dihedral in molecule.dihedrals() {
                let (i, j, k, m) = (dihedral.i(), dihedral.j(), dihedral.k(), dihedral.m());
                let ri = new_position(universe, i, &idxes, newpos);
                let rj = new_position(universe, j, &idxes, newpos);
                let rk = new_position(universe, k, &idxes, newpos);
                let rm = new_position(universe, m, &idxes, newpos);
                let phi = universe.cell().dihedral(ri, rj, rk, rm);
                dihedrals += evaluator.dihedral(phi, i, j, k, m);
            }
        }

        let coulomb_delta = universe.coulomb_potential()
                                    .as_ref()
                                    .map_or(0.0, |coulomb|
            coulomb.move_particles_cost(universe, &idxes, newpos)
        );

        let mut global_delta = 0.0;
        for potential in universe.global_potentials() {
            global_delta += potential.move_particles_cost(universe, &idxes, newpos);
        }

        let cost = pairs_delta + (bonds - self.bonds)
                               + (angles - self.angles)
                               + (dihedrals - self.dihedrals)
                               + coulomb_delta + global_delta;

        self.updater = Some(Box::new(move |cache, universe| {
            cache.bonds = bonds;
            cache.angles = angles;
            cache.dihedrals = dihedrals;

            cache.pairs += pairs_delta;
            cache.coulomb += coulomb_delta;
            cache.global += global_delta;

            let (n, m) = new_pairs.size();
            debug_assert!(n == m);
            debug_assert!(n == cache.pairs_cache.size().0);
            debug_assert!(n == cache.pairs_cache.size().1);
            for i in 0..n {
                for j in 0..n {
                    if new_pairs[(i, j)] != 0.0 {
                        cache.pairs_cache[(i, j)] = new_pairs[(i, j)];
                    }
                }
            }
            // Update the cache for the global potentials
            universe.coulomb_potential_mut().as_mut().and_then(|coulomb| Some(coulomb.update()));
            for potential in universe.global_potentials_mut() {
                potential.update();
            }
        }));
        return cost;
    }
}

/// Return either the new position of a particle (from `newpos`) if its index
/// is in `idxes`, or its old position in the universe.
fn new_position<'a>(universe: &'a Universe, i: usize, idxes: &[usize], newpos: &'a[Vector3D]) -> &'a Vector3D {
    for (idx, &particle) in idxes.iter().enumerate() {
        if particle == i {
            return &newpos[idx];
        }
    }
    return &universe[i].position;
}

#[cfg(test)]
mod tests {
    use super::*;
    use universe::{Universe, UnitCell};
    use utils::universe_from_xyz;
    use input::read_interactions_string;
    use types::Vector3D;

    fn testing_universe() -> Universe {
        let mut universe = universe_from_xyz("4
        bonds
        O     0.000000     0.000000     0.000000
        O     0.000000     0.000000     1.480000
        H     0.895669     0.000000    -0.316667
        H    -0.895669     0.000000     1.796667");
        universe.set_cell(UnitCell::cubic(10.0));
        assert!(universe.molecules().len() == 1);

        read_interactions_string(&mut universe, "
        # Values are completely random, just having a bit of all the types
        pairs:
            - atoms: [H, H]
              type: LennardJones
              sigma: 3.0 A
              epsilon: 0.2 kJ/mol
            - atoms: [O, O]
              type: Null
            - atoms: [O, H]
              type: LennardJones
              sigma: 1.0 A
              epsilon: 0.33 kJ/mol
        bonds:
            - atoms: [O, H]
              type: Harmonic
              x0: 3.4 A
              k: 522 kJ/mol/A^2
        angles:
            - atoms: [O, O, H]
              type: Harmonic
              x0: 120 deg
              k: 150 kJ/mol/deg^2
        dihedrals:
            - atoms: [H, O, O, H]
              type: Harmonic
              x0: 180 deg
              k: 800 kJ/mol/deg^2
        coulomb:
            type: Wolf
            cutoff: 8 A
            charges:
                O: -0.5
                H: 0.5
        ").unwrap();

        return universe;
    }

    #[test]
    fn cache_energy() {
        let universe = testing_universe();
        let mut cache = EnergyCache::new();
        cache.init(&universe);

        assert_approx_eq!(cache.energy(), universe.potential_energy());
    }

    #[test]
    fn move_atoms() {
        let mut universe = testing_universe();
        let mut cache = EnergyCache::new();
        let old_e = universe.potential_energy();
        cache.init(&universe);

        let idxes = vec![0, 3];
        let newpos = &[Vector3D::new(0.0, 0.0, 0.5), Vector3D::new(-0.7, 0.2, 1.5)];

        let cost = cache.move_particles_cost(&universe, idxes, newpos);

        universe[0].position = newpos[0];
        universe[3].position = newpos[1];
        let new_e = universe.potential_energy();
        assert_approx_eq!(cost, new_e - old_e);

        cache.update(&mut universe);
        assert_approx_eq!(cache.energy(), new_e);

        // Check that the cache is really updated
        let old_e = new_e;
        let idxes = vec![2, 3];
        let newpos = &[Vector3D::new(0.9, 0.2, -0.4), Vector3D::new(-0.9, 0.0, 1.8)];
        let cost = cache.move_particles_cost(&universe, idxes, newpos);
        universe[2].position = newpos[0];
        universe[3].position = newpos[1];
        let new_e = universe.potential_energy();
        assert_approx_eq!(cost, new_e - old_e);
    }
}
