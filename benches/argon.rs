// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 Lumol's contributors — BSD license

#[macro_use]
extern crate bencher;
extern crate rand;
extern crate lumol;
extern crate lumol_input;

use bencher::Bencher;
use rand::Rng;

use lumol::sys::EnergyCache;
use lumol::types::Vector3D;

mod utils;


fn energy(bencher: &mut Bencher) {
    let system = utils::get_system("argon");
    bencher.iter(||{
        let _ = system.potential_energy();
    })
}

fn forces(bencher: &mut Bencher) {
    let system = utils::get_system("argon");
    bencher.iter(||{
        let _ = system.forces();
    })
}

fn virial(bencher: &mut Bencher) {
    let system = utils::get_system("argon");
    bencher.iter(||{
        let _ = system.virial();
    })
}

fn cache_move_particle(bencher: &mut Bencher) {
    let system = utils::get_system("argon");
    let mut cache = EnergyCache::new();
    cache.init(&system);

    let mut rng = rand::weak_rng();

    let particle: usize = rng.gen_range(0, system.size());
    let mut delta = system[particle].position;
    delta += Vector3D::new(rng.gen(), rng.gen(), rng.gen());

    bencher.iter(||{
        cache.move_particles_cost(&system, vec![particle], &[delta])
    })
}

benchmark_group!(energy_computation, energy, forces, virial);
benchmark_group!(monte_carlo_cache, cache_move_particle);

benchmark_main!(energy_computation, monte_carlo_cache);
