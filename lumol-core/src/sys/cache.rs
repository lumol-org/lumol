// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

//! Caching energy components for speeding up the energy computation in
//! Monte Carlo simulations.
//!
//! In most of Monte Carlo moves, only a very small subset of the system changes.
//! We can use that property to remove the need of recomputing most of the
//! energy components, by storing them and providing update callbacks.

use crate::System;
use crate::{Array2, Vector3D};

/// Callback for updating a cache. It also take an `&mut System` argument for
/// updating the cache inside the global potentials.
type UpdateCallback = Box<dyn Fn(&mut EnergyCache, &mut System) + Send + Sync>;

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
    /// Contribution of long range corrections
    pairs_tail: f64,
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
    updater: Option<UpdateCallback>,
}

impl EnergyCache {
    /// Create a new empty energy cache.
    pub fn new() -> EnergyCache {
        EnergyCache {
            pairs_cache: Array2::zeros((0, 0)),
            pairs: 0.0,
            pairs_tail: 0.0,
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
        self.pairs_tail = 0.0;
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
                let r = system.nearest_image(i, j).norm();
                let path = system.bond_path(i, j);
                let energy = evaluator.pair(path, r, i, j);
                self.pairs_cache[(i, j)] = energy;
                self.pairs_cache[(j, i)] = energy;
                self.pairs += energy;
            }
        }

        self.pairs_tail = evaluator.pairs_tail();
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
        energy += self.pairs_tail;

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
        let updater = self.updater.take();
        if let Some(updater) = updater {
            updater(self, system);
        } else {
            panic!(
                "called EnergyCache::update without call a `*_cost` function first"
            );
        }
    }

    /// This function should be called whenever the cache is not used, but one
    /// still want it to be updated. Future call to `EnergyCache::update` will
    /// recompute the full cache.
    pub fn unused(&mut self) {
        self.updater = Some(Box::new(|cache, system| {
            cache.init(system);
        }));
    }
}

impl EnergyCache {
    /// Get the cost of moving a rigid molecule at `molecule_id` in the system
    /// to `new_positions`.
    ///
    /// This function ***DOES NOT*** update the cache, the `update` function
    /// MUST be called if the particles are effectively moved.
    pub fn move_molecule_cost(
        &mut self,
        system: &System,
        molecule_id: usize,
        new_positions: &[Vector3D],
    ) -> f64 {
        let evaluator = system.energy_evaluator();
        let positions = system.particles().position;
        let molecule = system.molecule(molecule_id);

        let mut new_pairs = Array2::<f64>::zeros((system.size(), system.size()));
        let mut pairs_delta = 0.0;

        // Iterate over all interactions between a particle in the moved
        // molecule and a particle in another molecule
        for (i, part_i) in molecule.indexes().enumerate() {
            for (_, other_molecule) in system.molecules().enumerate().filter(|(id, _)| molecule_id != *id) {
                for part_j in other_molecule.indexes() {
                    let r = system.cell.distance(&positions[part_j], &new_positions[i]);
                    let path = system.bond_path(part_i, part_j);
                    let energy = evaluator.pair(path, r, part_i, part_j);

                    pairs_delta += energy;
                    new_pairs[(part_i, part_j)] += energy;
                    new_pairs[(part_j, part_i)] += energy;

                    pairs_delta -= self.pairs_cache[(part_i, part_j)];
                }
            }
        }

        // Pairs tail correction do not change when moving a single molecule

        // Bonds / Angles / Dihedrals terms do not change

        let coulomb_delta = system.coulomb_potential()
            .map_or(0.0, |coulomb| coulomb.move_molecule_cost(system, molecule_id, new_positions));

        let mut global_delta = 0.0;
        for global in system.global_potentials() {
            global_delta += global.move_molecule_cost(system, molecule_id, new_positions);
        }

        let cost = pairs_delta + coulomb_delta + global_delta;

        self.updater = Some(Box::new(move |cache, system| {
            cache.pairs += pairs_delta;
            cache.coulomb += coulomb_delta;
            cache.global += global_delta;

            let (n, m) = new_pairs.dim();
            debug_assert_eq!(n, m);
            debug_assert_eq!((n, m), cache.pairs_cache.dim());

            let molecule = system.molecule(molecule_id);
            for i in molecule.indexes() {
                for j in 0..n {
                    if molecule.contains(j) {
                        continue;
                    }
                    cache.pairs_cache[(i, j)] = new_pairs[(i, j)];
                    cache.pairs_cache[(j, i)] = new_pairs[(i, j)];
                }
            }

            // Update the cache for the global potentials
            if let Some(coulomb) = system.coulomb_potential() {
                coulomb.update();
            }

            for global in system.global_potentials() {
                global.update();
            }
        }));
        return cost;
    }

    /// Return the cost for moving all **rigid** molecules of the system.
    ///
    /// This function is intended for use when all the molecules in the system
    /// are moved rigidly, for example when resizing the system in NPT
    /// Monte Carlo. It computes energy changes due to:
    ///
    /// - non bonded pairs interactions;
    /// - Coulomb interactions;
    /// - global interactions;
    ///
    /// It **DOES NOT** recompute bonds, angles and dihedral interactions. You
    /// must not use this function when the intramolecular configuration
    /// changed.
    ///
    /// This function ***DOES NOT*** update the cache, the `update` function
    /// MUST be called if the molecules are effectively moved.
    pub fn move_all_molecules_cost(&mut self, system: &System) -> f64 {
        let evaluator = system.energy_evaluator();

        let mut new_pairs = Array2::<f64>::zeros((system.size(), system.size()));
        let mut pairs_delta = 0.0;
        // Loop over all molecule pairs
        for (i, mol_i) in system.molecules().enumerate() {
            for mol_j in system.molecules().skip(i + 1) {
                // Loop over all particles in the molecules
                for part_i in mol_i.indexes() {
                    for part_j in mol_j.indexes() {
                        let r = system.distance(part_i, part_j);
                        let path = system.bond_path(part_i, part_j);
                        let energy = evaluator.pair(path, r, part_i, part_j);
                        pairs_delta += energy;
                        new_pairs[(part_i, part_j)] += energy;
                        new_pairs[(part_j, part_i)] += energy;
                        pairs_delta -= self.pairs_cache[(part_i, part_j)];
                    }
                }
            }
        }

        // temporarily, recompute all interactions
        let new_coulomb = evaluator.coulomb();
        let new_global = evaluator.global();

        // compute the new tail correction
        let pairs_tail = evaluator.pairs_tail();

        let cost = pairs_delta + (pairs_tail - self.pairs_tail) + (new_coulomb - self.coulomb)
            + (new_global - self.global);

        self.updater = Some(Box::new(move |cache, system| {
            cache.pairs += pairs_delta;
            cache.pairs_tail = pairs_tail;
            cache.coulomb = new_coulomb;
            cache.global = new_global;

            let (n, m) = new_pairs.dim();
            debug_assert_eq!(n, m);
            debug_assert_eq!((n, m), cache.pairs_cache.dim());
            for (i, mol_i) in system.molecules().enumerate() {
                for mol_j in system.molecules().skip(i + 1) {
                    for part_i in mol_i.indexes() {
                        for part_j in mol_j.indexes() {
                            cache.pairs_cache[(part_i, part_j)] = new_pairs[(part_i, part_j)];
                            cache.pairs_cache[(part_j, part_i)] = new_pairs[(part_i, part_j)];
                        }
                    }
                }
            }
        }));
        cost
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Harmonic, LennardJones, NullPotential, Wolf};
    use crate::PairInteraction;
    use crate::System;
    use crate::Vector3D;
    use crate::utils::system_from_xyz;
    use crate::units;

    use approx::{assert_ulps_eq, assert_relative_eq};

    fn testing_system() -> System {
        let mut system = system_from_xyz(
            "8
            cell: 10.0
            H     0.895669     0.000000    -0.316667
            O     0.000000     0.000000     0.000000
            O     0.000000     0.000000     1.480000
            H    -0.895669     0.000000     1.796667
            O     3.000000     0.000000     0.000000
            O     3.000000     0.000000     1.480000
            H     3.895669     0.000000    -0.316667
            H     2.104330     0.000000     1.796667",
        );
        assert!(system.add_bond(0, 1).is_empty());
        assert!(system.add_bond(0, 2).is_empty());
        assert!(system.add_bond(1, 3).is_empty());

        assert!(system.add_bond(4, 5).is_empty());
        assert!(system.add_bond(4, 6).is_empty());
        assert!(system.add_bond(5, 7).is_empty());
        assert!(system.molecules().count() == 2);

        system.set_pair_potential(
            ("H", "H"),
            PairInteraction::new(
                Box::new(LennardJones {
                    sigma: 3.0,
                    epsilon: units::from(0.5, "kJ/mol").unwrap(),
                }),
                3.0,
            ),
        );

        system.set_pair_potential(("O", "O"), PairInteraction::new(Box::new(NullPotential), 3.0));

        system.set_pair_potential(
            ("O", "H"),
            PairInteraction::new(
                Box::new(LennardJones {
                    sigma: 1.0,
                    epsilon: units::from(0.3, "kJ/mol").unwrap(),
                }),
                3.0,
            ),
        );

        system.set_bond_potential(
            ("O", "O"),
            Box::new(Harmonic {
                x0: 2.4,
                k: units::from(522.0, "kJ/mol/A^2").unwrap(),
            }),
        );

        system.set_bond_potential(
            ("O", "H"),
            Box::new(Harmonic {
                x0: 1.4,
                k: units::from(122.0, "kJ/mol/A^2").unwrap(),
            }),
        );

        system.set_angle_potential(
            ("O", "O", "H"),
            Box::new(Harmonic {
                x0: f64::to_radians(120.0),
                k: units::from(150.0, "kJ/mol/deg^2").unwrap(),
            }),
        );

        system.set_dihedral_potential(
            ("H", "O", "O", "H"),
            Box::new(Harmonic {
                x0: f64::to_radians(180.0),
                k: units::from(800.0, "kJ/mol/deg^2").unwrap(),
            }),
        );

        system.set_coulomb_potential(Box::new(Wolf::new(5.0)));

        for particle in system.particles_mut() {
            if particle.name == "O" {
                *particle.charge = -0.5;
            } else if particle.name == "H" {
                *particle.charge = 0.5;
            }
        }
        system
    }

    #[test]
    fn cache_energy() {
        let system = testing_system();
        let mut cache = EnergyCache::new();
        cache.init(&system);

        assert_ulps_eq!(cache.energy(), system.potential_energy());
    }

    #[test]
    #[allow(clippy::unreadable_literal)]
    fn move_molecule() {
        let mut system = testing_system();
        let mut cache = EnergyCache::new();
        let old_energy = system.potential_energy();
        cache.init(&system);
        assert_ulps_eq!(cache.energy(), old_energy);

        let new_positions = &[
            Vector3D::new(-0.987061, 0.59401, 0.427533),
            Vector3D::new(-1.0744137409578138, 1.2111820514074991, -0.2893833856814936),
            Vector3D::new(-1.4352068561309008, 2.5425486908430286, 0.24698514382209652),
            Vector3D::new(-1.5225595970887147, 3.159720742250528, -0.46993124185939705),
        ];
        let cost = cache.move_molecule_cost(&system, 0, new_positions);

        system.particles_mut().position[0] = new_positions[0];
        system.particles_mut().position[1] = new_positions[1];
        system.particles_mut().position[2] = new_positions[2];
        system.particles_mut().position[3] = new_positions[3];
        let new_energy = system.potential_energy();
        assert_relative_eq!(cost, new_energy - old_energy, max_relative = 1e-9);

        cache.update(&mut system);
        assert_ulps_eq!(cache.energy(), new_energy);

        // Check that the cache is really updated
        let old_energy = new_energy;
        let new_positions = &[
            Vector3D::new(-0.49138099999999996, 1.08969, 0.923213),
            Vector3D::new(-1.0139188773839494, 1.7555242433058806, 1.3546291257885033),
            Vector3D::new(-2.4014453359903767, 1.2663193861239206, 1.5154051637316108),
            Vector3D::new(-2.923983213374326, 1.9321536294298012, 1.946821289520114),
        ];
        let cost = cache.move_molecule_cost(&system, 0, new_positions);
        system.particles_mut().position[0] = new_positions[0];
        system.particles_mut().position[1] = new_positions[1];
        system.particles_mut().position[2] = new_positions[2];
        system.particles_mut().position[3] = new_positions[3];
        let new_energy = system.potential_energy();
        assert_relative_eq!(cost, new_energy - old_energy, max_relative = 1e-9);
    }

    #[test]
    fn move_all_molecules() {
        let system = testing_system();
        let mut cache = EnergyCache::new();
        let old_energy = system.potential_energy();
        cache.init(&system);

        let delta = Vector3D::new(1.0, 0.5, -0.5);

        let mut new_system = system.clone();
        // translate the center of mass
        for position in new_system.molecule_mut(0).particles_mut().position {
            *position += delta;
        }
        let cost = cache.move_all_molecules_cost(&new_system);
        let new_energy = new_system.potential_energy();
        assert_ulps_eq!(cost, new_energy - old_energy, epsilon = 1e-12);
        cache.update(&mut new_system);
        assert_ulps_eq!(cache.energy(), new_energy, epsilon = 1e-12);

        // Check that the cache is really updated
        // move the other molecule
        let old_energy = new_energy;
        let delta = Vector3D::new(-0.9, 0.0, 1.8);
        let mut new_system = system;
        // translate the center of mass
        for mut molecule in new_system.molecules_mut() {
            for position in molecule.particles_mut().position {
                *position += delta;
            }
        }
        let cost = cache.move_all_molecules_cost(&new_system);
        let new_energy = new_system.potential_energy();
        assert_ulps_eq!(cost, new_energy - old_energy, epsilon = 1e-12);
    }
}
