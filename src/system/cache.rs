// Cymbalum, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

//! Caching energy components for speeding up the energy computation in
//! Monte-Carlo simulations.
//!
//! In most of Monte-Carlo moves, only a very small subset of the system changes.
//! We can use that property to remove the need of recomputing most of the
//! energy components, by storing them and providing update callbacks.
use std::mem;

use super::System;
use types::{Vector3D, Array2};
use potentials::global::GlobalCache;

/// Callback for updating a cache. It also take an `&mut System` argument for
/// updating the cache inside the global potentials.
type UpdateCallback = Box<Fn(&mut EnergyCache, &mut System)>;

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
    /// Callback to be called to update the cache if the system is modified
    updater: Option<UpdateCallback>
}

impl EnergyCache {
    /// Create a new empty energy cache.
    pub fn new() -> EnergyCache {
        EnergyCache {
            pairs_cache: Array2::zeros((0, 0)),
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
        self.pairs_cache.assign(0.0);
        self.pairs = 0.0;
        self.bonds = 0.0;
        self.angles = 0.0;
        self.dihedrals = 0.0;
        self.coulomb = 0.0;
        self.global = 0.0;
    }

    /// Initialize the cache to be used with `system`. After a call to this
    /// function, the cache is only usable with the same system. To change
    /// the associated system, one must call this function again.
    pub fn init(&mut self, system: &System) {
        self.clear();
        self.pairs_cache.resize_if_different((system.size(), system.size()));

        let evaluator = system.energy_evaluator();

        for i in 0..system.size() {
            for j in (i + 1)..system.size() {
                let r = system.wraped_vector(i, j).norm();
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
    pub fn update(&mut self, system: &mut System) {
        let updater = mem::replace(&mut self.updater, None);
        if let Some(updater) = updater {
            updater(self, system);
        } else {
            fatal_error!("Error in `EnergyCache::update`: \
                         This function MUST be called after a call to a `*_cost` function");
        }
    }

    /// This function should be called whenever the cache is not used, but one
    /// still want it to be updated. Future call to `EnergyCache::update` will
    /// recompute the full cache.
    pub fn unused(&mut self) {
        self.updater = Some(Box::new(|cache, system| {
            cache.init(system);
        }))
    }
}

impl EnergyCache {
    /// Get the cost of moving the set of particles with indexes in `idxes` to
    /// `newpos`. This function DO NOT update the cache, the
    /// `update_particles_moved` function MUST be called if the particles are
    /// effectively moved.
    pub fn move_particles_cost(&mut self, system: &System, idxes: Vec<usize>, newpos: &[Vector3D]) -> f64 {
        let evaluator = system.energy_evaluator();

        // First, go for pair interactions
        let mut new_pairs = Array2::<f64>::zeros((system.size(), system.size()));
        let mut pairs_delta = 0.0;
        // Interactions with the sub-system not being moved
        for (i, &part_i) in idxes.iter().enumerate() {
            for part_j in 0..system.size() {
                // Exclude interactions inside the sub-system.
                if idxes.contains(&part_j) {continue}

                let r = system.cell().distance(&system[part_j].position, &newpos[i]);
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
                let r = system.cell().distance(&newpos[i], &newpos[j]);
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
        for molecule in system.molecules() {
            for bond in molecule.bonds() {
                let (i, j) = (bond.i(), bond.j());
                let ri = new_position(system, i, &idxes, newpos);
                let rj = new_position(system, j, &idxes, newpos);
                let r = system.cell().distance(ri, rj);
                bonds += evaluator.bond(r, i, j);
            }

            for angle in molecule.angles() {
                let (i, j, k) = (angle.i(), angle.j(), angle.k());
                let ri = new_position(system, i, &idxes, newpos);
                let rj = new_position(system, j, &idxes, newpos);
                let rk = new_position(system, k, &idxes, newpos);
                let theta = system.cell().angle(ri, rj, rk);
                angles += evaluator.angle(theta, i, j, k);
            }

            for dihedral in molecule.dihedrals() {
                let (i, j, k, m) = (dihedral.i(), dihedral.j(), dihedral.k(), dihedral.m());
                let ri = new_position(system, i, &idxes, newpos);
                let rj = new_position(system, j, &idxes, newpos);
                let rk = new_position(system, k, &idxes, newpos);
                let rm = new_position(system, m, &idxes, newpos);
                let phi = system.cell().dihedral(ri, rj, rk, rm);
                dihedrals += evaluator.dihedral(phi, i, j, k, m);
            }
        }

        let coulomb_delta = system.coulomb_potential()
                                    .as_ref()
                                    .map_or(0.0, |coulomb|
            coulomb.borrow_mut().move_particles_cost(system, &idxes, newpos)
        );

        let mut global_delta = 0.0;
        for potential in system.global_potentials() {
            global_delta += potential.borrow_mut().move_particles_cost(system, &idxes, newpos);
        }

        let cost = pairs_delta + (bonds - self.bonds)
                               + (angles - self.angles)
                               + (dihedrals - self.dihedrals)
                               + coulomb_delta + global_delta;

        self.updater = Some(Box::new(move |cache, system| {
            cache.bonds = bonds;
            cache.angles = angles;
            cache.dihedrals = dihedrals;

            cache.pairs += pairs_delta;
            cache.coulomb += coulomb_delta;
            cache.global += global_delta;

            let (n, m) = new_pairs.shape();
            debug_assert!(n == m);
            debug_assert!((n, m) == cache.pairs_cache.shape());
            for i in 0..n {
                for j in 0..n {
                    if new_pairs[(i, j)] != 0.0 {
                        cache.pairs_cache[(i, j)] = new_pairs[(i, j)];
                    }
                }
            }
            // Update the cache for the global potentials
            if let Some(coulomb) = system.coulomb_potential() {
                coulomb.borrow_mut().update();
            }
            for potential in system.global_potentials() {
                potential.borrow_mut().update();
            }
        }));
        return cost;
    }
}

/// Return either the new position of a particle (from `newpos`) if its index
/// is in `idxes`, or its old position in the system.
fn new_position<'a>(system: &'a System, i: usize, idxes: &[usize], newpos: &'a[Vector3D]) -> &'a Vector3D {
    for (idx, &particle) in idxes.iter().enumerate() {
        if particle == i {
            return &newpos[idx];
        }
    }
    return &system[i].position;
}

#[cfg(test)]
mod tests {
    use super::*;
    use system::{System, UnitCell};
    use utils::system_from_xyz;
    use input::read_interactions_string;
    use types::Vector3D;

    fn testing_system() -> System {
        let mut system = system_from_xyz("4
        bonds
        O     0.000000     0.000000     0.000000
        O     0.000000     0.000000     1.480000
        H     0.895669     0.000000    -0.316667
        H    -0.895669     0.000000     1.796667");
        system.set_cell(UnitCell::cubic(10.0));
        assert!(system.molecules().len() == 1);

        read_interactions_string(&mut system, "
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

        return system;
    }

    #[test]
    fn cache_energy() {
        let system = testing_system();
        let mut cache = EnergyCache::new();
        cache.init(&system);

        assert_approx_eq!(cache.energy(), system.potential_energy());
    }

    #[test]
    fn move_atoms() {
        let mut system = testing_system();
        let mut cache = EnergyCache::new();
        let old_e = system.potential_energy();
        cache.init(&system);

        let idxes = vec![0, 3];
        let newpos = &[Vector3D::new(0.0, 0.0, 0.5), Vector3D::new(-0.7, 0.2, 1.5)];

        let cost = cache.move_particles_cost(&system, idxes, newpos);

        system[0].position = newpos[0];
        system[3].position = newpos[1];
        let new_e = system.potential_energy();
        assert_approx_eq!(cost, new_e - old_e);

        cache.update(&mut system);
        assert_approx_eq!(cache.energy(), new_e);

        // Check that the cache is really updated
        let old_e = new_e;
        let idxes = vec![2, 3];
        let newpos = &[Vector3D::new(0.9, 0.2, -0.4), Vector3D::new(-0.9, 0.0, 1.8)];
        let cost = cache.move_particles_cost(&system, idxes, newpos);
        system[2].position = newpos[0];
        system[3].position = newpos[1];
        let new_e = system.potential_energy();
        assert_approx_eq!(cost, new_e - old_e);
    }
}
