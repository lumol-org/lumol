// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

#[macro_use]
extern crate bencher;
extern crate lumol;
extern crate lumol_input;
extern crate rand;

#[macro_use]
mod utils;

mod ewald {
    use bencher::Bencher;
    use rand::Rng;

    use lumol::energy::{CoulombicPotential, Ewald, GlobalPotential, PairRestriction, SharedEwald};
    use lumol::sys::EnergyCache;
    use lumol::types::{Vector3D, Zero};

    use utils;

    fn get_ewald() -> SharedEwald {
        let mut ewald = SharedEwald::new(Ewald::new(8.0, 7, None));
        ewald.set_restriction(PairRestriction::InterMolecular);
        ewald
    }


    pub fn energy(bencher: &mut Bencher) {
        let system = utils::get_system("water");
        let ewald = get_ewald();
        ewald.energy(&system);

        bencher.iter(|| {
            let _ = ewald.energy(&system);
        })
    }

    pub fn forces(bencher: &mut Bencher) {
        let system = utils::get_system("water");
        let ewald = get_ewald();
        let mut forces = vec![Vector3D::zero(); system.size()];
        ewald.forces(&system, &mut forces);

        bencher.iter(|| {
            ewald.forces(&system, &mut forces);
        })
    }

    pub fn virial(bencher: &mut Bencher) {
        let system = utils::get_system("water");
        let ewald = get_ewald();
        ewald.virial(&system);

        bencher.iter(|| {
            let _ = ewald.virial(&system);
        })
    }

    pub fn cache_move_particles(bencher: &mut Bencher) {
        let mut system = utils::get_system("water");
        system.set_coulomb_potential(Box::new(get_ewald()));

        let mut cache = EnergyCache::new();
        cache.init(&system);

        let mut rng = utils::get_rng([
            206, 1, 245, 36, 62, 147, 30, 213, 177, 131, 94, 148, 239, 154, 161, 1
        ]);

        let molid = rng.gen_range(0, system.molecules_count());
        let molecule = system.molecule(molid);
        let mut delta = vec![];
        for position in molecule.particles().position {
            delta.push(position + Vector3D::new(rng.gen(), rng.gen(), rng.gen()));
        }

        cache.move_particles_cost(&system, molecule.indexes().collect(), &delta);

        bencher.iter(|| cache.move_particles_cost(&system, molecule.indexes().collect(), &delta))
    }

    pub fn cache_move_all_rigid_molecules(bencher: &mut Bencher) {
        let mut system = utils::get_system("water");
        system.set_coulomb_potential(Box::new(get_ewald()));

        let mut cache = EnergyCache::new();
        cache.init(&system);

        let mut rng = utils::get_rng([
            79, 129, 118, 38, 44, 204, 227, 6, 233, 6, 7, 216, 192, 77, 33, 85
        ]);

        for mut molecule in system.molecules_mut() {
            let delta = Vector3D::new(rng.gen(), rng.gen(), rng.gen());
            for position in molecule.particles_mut().position {
                *position += delta;
            }
        }

        cache.move_all_rigid_molecules_cost(&system);

        bencher.iter(|| cache.move_all_rigid_molecules_cost(&system))
    }
}

mod wolf {
    use bencher::Bencher;
    use rand::Rng;

    use lumol::energy::{CoulombicPotential, GlobalPotential, PairRestriction, Wolf};
    use lumol::sys::EnergyCache;
    use lumol::types::{Vector3D, Zero};

    use utils;

    fn get_wolf() -> Wolf {
        let mut wolf = Wolf::new(9.0);
        wolf.set_restriction(PairRestriction::InterMolecular);
        wolf
    }

    pub fn energy(bencher: &mut Bencher) {
        let system = utils::get_system("water");
        let wolf = get_wolf();
        wolf.energy(&system);

        bencher.iter(|| {
            let _ = wolf.energy(&system);
        })
    }

    pub fn forces(bencher: &mut Bencher) {
        let system = utils::get_system("water");
        let wolf = get_wolf();
        let mut forces = vec![Vector3D::zero(); system.size()];
        wolf.forces(&system, &mut forces);

        bencher.iter(|| {
            wolf.forces(&system, &mut forces);
        })
    }

    pub fn virial(bencher: &mut Bencher) {
        let system = utils::get_system("water");
        let wolf = get_wolf();
        wolf.virial(&system);

        bencher.iter(|| {
            let _ = wolf.virial(&system);
        })
    }

    pub fn cache_move_particles(bencher: &mut Bencher) {
        let mut system = utils::get_system("water");
        system.set_coulomb_potential(Box::new(get_wolf()));

        let mut cache = EnergyCache::new();
        cache.init(&system);

        let mut rng = utils::get_rng([
            215, 235, 194, 22, 205, 151, 210, 241, 188, 67, 241, 2, 204, 62, 11, 201
        ]);

        let molid = rng.gen_range(0, system.molecules_count());
        let molecule = system.molecule(molid);
        let mut delta = vec![];
        for position in molecule.particles().position {
            delta.push(position + Vector3D::new(rng.gen(), rng.gen(), rng.gen()));
        }

        cache.move_particles_cost(&system, molecule.indexes().collect(), &delta);

        bencher.iter(|| cache.move_particles_cost(&system, molecule.indexes().collect(), &delta))
    }

    pub fn cache_move_all_rigid_molecules(bencher: &mut Bencher) {
        let mut system = utils::get_system("water");
        system.set_coulomb_potential(Box::new(get_wolf()));

        let mut cache = EnergyCache::new();
        cache.init(&system);

        let mut rng = utils::get_rng([
            89, 208, 141, 72, 208, 131, 249, 179, 77, 243, 111, 32, 176, 194, 79, 44
        ]);

        for mut molecule in system.molecules_mut() {
            let delta = Vector3D::new(rng.gen(), rng.gen(), rng.gen());
            for position in molecule.particles_mut().position {
                *position += delta;
            }
        }

        cache.move_all_rigid_molecules_cost(&system);

        bencher.iter(|| cache.move_all_rigid_molecules_cost(&system))
    }
}

benchmark_group!(ewald, ewald::energy, ewald::forces, ewald::virial);
benchmark_group!(wolf, wolf::energy, wolf::forces, wolf::virial);
benchmark_group!(monte_carlo_cache,
    ewald::cache_move_particles,
    ewald::cache_move_all_rigid_molecules,
    wolf::cache_move_particles,
    wolf::cache_move_all_rigid_molecules
);

benchmark_main!(ewald, wolf, monte_carlo_cache);
