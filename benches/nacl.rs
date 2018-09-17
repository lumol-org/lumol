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

    use lumol::energy::{GlobalPotential, Ewald, SharedEwald};
    use lumol::sys::EnergyCache;
    use lumol::types::{Vector3D, Zero};

    use utils;

    pub fn energy(bencher: &mut Bencher) {
        let system = utils::get_system("nacl");
        let ewald = SharedEwald::new(Ewald::new(9.5, 7, None));
        ewald.energy(&system);

        bencher.iter(|| {
            let _ = ewald.energy(&system);
        })
    }

    pub fn forces(bencher: &mut Bencher) {
        let system = utils::get_system("nacl");
        let ewald = SharedEwald::new(Ewald::new(9.5, 7, None));
        let mut forces = vec![Vector3D::zero(); system.size()];
        ewald.forces(&system, &mut forces);

        bencher.iter(|| {
            ewald.forces(&system, &mut forces);
        })
    }

    pub fn virial(bencher: &mut Bencher) {
        let system = utils::get_system("nacl");
        let ewald = SharedEwald::new(Ewald::new(9.5, 7, None));
        ewald.virial(&system);

        bencher.iter(|| {
            let _ = ewald.virial(&system);
        })
    }

    pub fn cache_move_particle(bencher: &mut Bencher) {
        let mut system = utils::get_system("nacl");
        system.set_coulomb_potential(Box::new(SharedEwald::new(Ewald::new(9.5, 7, None))));

        let mut cache = EnergyCache::new();
        cache.init(&system);

        let mut rng = utils::get_rng([
            228, 140, 229, 238, 195, 151, 56, 106, 68, 11, 24, 143, 231, 175, 55, 52
        ]);

        let molid = rng.gen_range(0, system.size());
        let mut new_position = system.particles().position[molid];
        new_position += Vector3D::new(rng.gen(), rng.gen(), rng.gen());

        cache.move_molecule_cost(&system, molid, &[new_position]);
        bencher.iter(|| cache.move_molecule_cost(&system, molid, &[new_position]))
    }
}

mod wolf {
    use bencher::Bencher;
    use rand::Rng;

    use lumol::energy::{GlobalPotential, Wolf};
    use lumol::sys::EnergyCache;
    use lumol::types::{Vector3D, Zero};

    use utils;

    pub fn energy(bencher: &mut Bencher) {
        let system = utils::get_system("nacl");
        let wolf = Wolf::new(12.0);
        wolf.energy(&system);

        bencher.iter(|| {
            let _ = wolf.energy(&system);
        })
    }

    pub fn forces(bencher: &mut Bencher) {
        let system = utils::get_system("nacl");
        let wolf = Wolf::new(12.0);
        let mut forces = vec![Vector3D::zero(); system.size()];
        wolf.forces(&system, &mut forces);

        bencher.iter(|| {
            wolf.forces(&system, &mut forces);
        })
    }

    pub fn virial(bencher: &mut Bencher) {
        let system = utils::get_system("nacl");
        let wolf = Wolf::new(12.0);
        wolf.virial(&system);

        bencher.iter(|| {
            let _ = wolf.virial(&system);
        })
    }

    pub fn cache_move_particle(bencher: &mut Bencher) {
        let mut system = utils::get_system("nacl");
        system.set_coulomb_potential(Box::new(Wolf::new(12.0)));

        let mut cache = EnergyCache::new();
        cache.init(&system);

        let mut rng = utils::get_rng([
            12, 197, 68, 124, 239, 99, 89, 228, 140, 170, 228, 215, 97, 218, 201, 24
        ]);

        let molid = rng.gen_range(0, system.size());
        let mut new_position = system.particles().position[molid];
        new_position += Vector3D::new(rng.gen(), rng.gen(), rng.gen());

        cache.move_molecule_cost(&system, molid, &[new_position]);
        bencher.iter(|| cache.move_molecule_cost(&system, molid, &[new_position]))
    }
}

benchmark_group!(ewald, ewald::energy, ewald::forces, ewald::virial);
benchmark_group!(wolf, wolf::energy, wolf::forces, wolf::virial);
benchmark_group!(monte_carlo_cache, ewald::cache_move_particle, wolf::cache_move_particle);

benchmark_main!(ewald, wolf, monte_carlo_cache);
