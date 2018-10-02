// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

#[macro_use]
extern crate criterion;
extern crate rand;
extern crate lumol;

use criterion::Criterion;

use lumol::sys::EnergyCache;
use lumol::sys::compute::{MolecularVirial, AtomicVirial, PotentialEnergy, Forces, Compute};

mod utils;

fn energy_computation(c: &mut Criterion) {
    let system = utils::get_system("propane");
    c.bench_function("propane::energy", move |b| b.iter(|| {
        let _ = PotentialEnergy.compute(&system);
    }));

    let system = utils::get_system("propane");
    c.bench_function("propane::force", move |b| b.iter(|| {
        let _ = Forces.compute(&system);
    }));

    let system = utils::get_system("propane");
    c.bench_function("propane::atomic_virial", move |b| b.iter(|| {
        let _ = AtomicVirial.compute(&system);
    }));

    let system = utils::get_system("propane");
    c.bench_function("propane::molecular_virial", move |b| b.iter(|| {
        let _ = MolecularVirial.compute(&system);
    }));
}

fn monte_carlo_cache(c: &mut Criterion) {
    let system = utils::get_system("propane");
    let mut cache = EnergyCache::new();
    cache.init(&system);

    c.bench_function("propane::move_molecule_cost", move |b| b.iter_with_setup(
        || utils::move_rigid_molecule(&system),
        |(molid, positions)| cache.move_molecule_cost(&system, molid, &positions)
    ));

    let mut system = utils::get_system("propane");
    let mut cache = EnergyCache::new();
    cache.init(&system);

    c.bench_function("propane::move_all_molecules_cost", move |b| b.iter_with_setup(
        || utils::move_all_rigid_molecule(&mut system),
        |system| cache.move_all_molecules_cost(&system)
    ));
}

criterion_group!(propane, energy_computation, monte_carlo_cache);
criterion_main!(propane);
