// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

#[macro_use]
extern crate bencher;
extern crate rand;
extern crate lumol;
extern crate lumol_input;

use bencher::Bencher;
use rand::Rng;

use lumol::energy::{Ewald, SharedEwald, Wolf, GlobalPotential};
use lumol::sys::EnergyCache;
use lumol::types::Vector3D;

#[macro_use]
mod utils;

fn energy_ewald(bencher: &mut Bencher) {
    let system = utils::get_system("nacl");
    let ewald = SharedEwald::new(Ewald::new(9.5, 7));

    bencher.iter(||{
        let _ = ewald.energy(&system);
    })
}

fn forces_ewald(bencher: &mut Bencher) {
    let system = utils::get_system("nacl");
    let ewald = SharedEwald::new(Ewald::new(9.5, 7));

    bencher.iter(||{
        let _ = ewald.forces(&system);
    })
}

fn virial_ewald(bencher: &mut Bencher) {
    let system = utils::get_system("nacl");
    let ewald = SharedEwald::new(Ewald::new(9.5, 7));

    bencher.iter(||{
        let _ = ewald.virial(&system);
    })
}

fn energy_wolf(bencher: &mut Bencher) {
    let system = utils::get_system("nacl");
    let wolf = Wolf::new(12.0);

    bencher.iter(||{
        let _ = wolf.energy(&system);
    })
}

fn forces_wolf(bencher: &mut Bencher) {
    let system = utils::get_system("nacl");
    let wolf = Wolf::new(12.0);

    bencher.iter(||{
        let _ = wolf.forces(&system);
    })
}

fn virial_wolf(bencher: &mut Bencher) {
    let system = utils::get_system("nacl");
    let wolf = Wolf::new(12.0);

    bencher.iter(||{
        let _ = wolf.virial(&system);
    })
}

fn cache_move_particle_ewald(bencher: &mut Bencher) {
    let mut system = utils::get_system("nacl");
    system.set_coulomb_potential(
        Box::new(SharedEwald::new(Ewald::new(9.5, 7)))
    );

    let mut cache = EnergyCache::new();
    cache.init(&system);

    let mut rng = utils::get_rng(41201154);

    let i: usize = rng.gen_range(0, system.size());
    let mut delta = system.particle(i).position;
    delta += Vector3D::new(rng.gen(), rng.gen(), rng.gen());

    bencher.iter(||{
        cache.move_particles_cost(&system, vec![i], &[delta])
    })
}

fn cache_move_particle_wolf(bencher: &mut Bencher) {
    let mut system = utils::get_system("nacl");
    system.set_coulomb_potential(Box::new(Wolf::new(12.0)));

    let mut cache = EnergyCache::new();
    cache.init(&system);

    let mut rng = utils::get_rng(474114);

    let i: usize = rng.gen_range(0, system.size());
    let mut delta = system.particle(i).position;
    delta += Vector3D::new(rng.gen(), rng.gen(), rng.gen());

    bencher.iter(||{
        cache.move_particles_cost(&system, vec![i], &[delta])
    })
}

benchmark_group!(ewald, energy_ewald, forces_ewald, virial_ewald);
benchmark_group!(wolf, energy_wolf, forces_wolf, virial_wolf);
benchmark_group!(monte_carlo_cache, cache_move_particle_ewald, cache_move_particle_wolf);

benchmark_main!(ewald, wolf, monte_carlo_cache);
