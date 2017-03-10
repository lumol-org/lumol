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
    let system = utils::get_system("propane");
    bencher.iter(||{
        let _ = system.potential_energy();
    })
}

fn forces(bencher: &mut Bencher) {
    let system = utils::get_system("propane");
    bencher.iter(||{
        let _ = system.forces();
    })
}

fn virial(bencher: &mut Bencher) {
    let system = utils::get_system("propane");
    bencher.iter(||{
        let _ = system.virial();
    })
}

fn cache_move_particles(bencher: &mut Bencher) {
    let system = utils::get_system("propane");
    let mut cache = EnergyCache::new();
    cache.init(&system);

    let mut rng = rand::weak_rng();

    let molecule = rng.choose(system.molecules()).unwrap();
    let mut delta = vec![];
    for i in molecule {
        let position = system[i].position;
        delta.push(position + Vector3D::new(rng.gen(), rng.gen(), rng.gen()));
    }

    bencher.iter(||{
        cache.move_particles_cost(&system, molecule.iter().collect(), &delta)
    })
}

fn cache_move_all_rigid_molecules(bencher: &mut Bencher) {
    let mut system = utils::get_system("propane");
    let mut cache = EnergyCache::new();
    cache.init(&system);

    let mut rng = rand::weak_rng();
    for molecule in system.molecules().to_owned() {
        let delta = Vector3D::new(rng.gen(), rng.gen(), rng.gen());
        for i in molecule {
            system[i].position += delta;
        }
    }

    bencher.iter(||{
        cache.move_all_rigid_molecules_cost(&system)
    })
}


benchmark_group!(energy_computation, energy, forces, virial);
benchmark_group!(monte_carlo_cache, cache_move_particles, cache_move_all_rigid_molecules);

benchmark_main!(energy_computation, monte_carlo_cache);
