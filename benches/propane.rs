// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

#[macro_use]
extern crate bencher;
extern crate lumol;
extern crate lumol_input;
extern crate rand;

use bencher::Bencher;
use rand::Rng;

use lumol::sys::EnergyCache;
use lumol::types::Vector3D;

#[macro_use]
mod utils;

fn energy(bencher: &mut Bencher) {
    let system = utils::get_system("propane");
    bencher.iter(|| {
        let _ = system.potential_energy();
    })
}

fn forces(bencher: &mut Bencher) {
    let system = utils::get_system("propane");
    bencher.iter(|| {
        let _ = system.forces();
    })
}

fn virial(bencher: &mut Bencher) {
    let system = utils::get_system("propane");
    bencher.iter(|| {
        let _ = system.virial();
    })
}

fn cache_move_particles(bencher: &mut Bencher) {
    let system = utils::get_system("propane");
    let mut cache = EnergyCache::new();
    cache.init(&system);

    let mut rng = utils::get_rng([
        145, 59, 58, 50, 238, 182, 97, 28, 107, 149, 227, 40, 90, 109, 196, 129
    ]);

    let molid = rng.gen_range(0, system.molecules_count());
    let molecule = system.molecule(molid);
    let mut delta = vec![];
    for position in molecule.particles().position {
        delta.push(position + Vector3D::new(rng.gen(), rng.gen(), rng.gen()));
    }

    bencher.iter(|| cache.move_particles_cost(&system, molecule.indexes().collect(), &delta))
}

fn cache_move_all_rigid_molecules(bencher: &mut Bencher) {
    let mut system = utils::get_system("propane");
    let mut cache = EnergyCache::new();
    cache.init(&system);

    let mut rng = utils::get_rng([
        106, 32, 216, 89, 97, 75, 51, 16, 90, 137, 27, 66, 167, 233, 109, 177
    ]);

    for mut molecule in system.molecules_mut() {
        let delta = Vector3D::new(rng.gen(), rng.gen(), rng.gen());
        for position in molecule.particles_mut().position {
            *position += delta;
        }
    }

    bencher.iter(|| cache.move_all_rigid_molecules_cost(&system))
}


benchmark_group!(energy_computation, energy, forces, virial);
benchmark_group!(monte_carlo_cache, cache_move_particles, cache_move_all_rigid_molecules);

benchmark_main!(energy_computation, monte_carlo_cache);
