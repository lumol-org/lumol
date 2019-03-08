// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license
use criterion::{Criterion, BatchSize, criterion_group, criterion_main};

use lumol::compute::{MolecularVirial, AtomicVirial, PotentialEnergy, Forces, Compute};
use lumol::EnergyCache;
mod utils;

fn energy_computation(c: &mut Criterion) {
    let system = utils::get_system("argon");
    c.bench_function("argon::energy", move |b| b.iter(|| {
        let _ = PotentialEnergy.compute(&system);
    }));

    let system = utils::get_system("argon");
    c.bench_function("argon::force", move |b| b.iter(|| {
        let _ = Forces.compute(&system);
    }));

    let system = utils::get_system("argon");
    c.bench_function("argon::atomic_virial", move |b| b.iter(|| {
        let _ = AtomicVirial.compute(&system);
    }));

    let system = utils::get_system("argon");
    c.bench_function("argon::molecular_virial", move |b| b.iter(|| {
        let _ = MolecularVirial.compute(&system);
    }));
}

fn monte_carlo_cache(c: &mut Criterion) {
    let system = utils::get_system("argon");
    let mut cache = EnergyCache::new();
    cache.init(&system);

    c.bench_function("argon::move_molecule_cost", move |b| b.iter_batched(
        || utils::move_rigid_molecule(&system),
        |(molid, positions)| cache.move_molecule_cost(&system, molid, &positions),
        BatchSize::SmallInput
    ));

    let system = utils::get_system("argon");
    let mut cache = EnergyCache::new();
    cache.init(&system);

    c.bench_function("argon::move_all_molecules_cost", move |b| b.iter_batched_ref(
        || utils::move_all_rigid_molecule(system.clone()),
        |system| cache.move_all_molecules_cost(system),
        BatchSize::SmallInput
    ));
}

criterion_group!(argon, energy_computation, monte_carlo_cache);
criterion_main!(argon);
