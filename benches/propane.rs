// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license
use criterion::{Criterion, BatchSize, criterion_group, criterion_main};

use lumol::EnergyCache;
use lumol::compute::{MolecularVirial, AtomicVirial, PotentialEnergy, Forces, Compute};

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

    c.bench_function("propane::move_molecule_cost", move |b| b.iter_batched(
        || utils::move_rigid_molecule(&system),
        |(molid, positions)| cache.move_molecule_cost(&system, molid, &positions),
        BatchSize::SmallInput
    ));

    let system = utils::get_system("propane");
    let mut cache = EnergyCache::new();
    cache.init(&system);

    c.bench_function("propane::move_all_molecules_cost", move |b| b.iter_batched_ref(
        || utils::move_all_rigid_molecule(system.clone()),
        |system| cache.move_all_molecules_cost(system),
        BatchSize::SmallInput
    ));
}

criterion_group!(propane, energy_computation, monte_carlo_cache);
criterion_main!(propane);
