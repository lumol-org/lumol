// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license
use criterion::{Criterion, BatchSize, criterion_group, criterion_main};

use lumol::{EnergyCache, Vector3D};
use lumol::{GlobalPotential, Ewald, SharedEwald, Wolf};

mod utils;

fn get_wolf() -> Wolf {
    Wolf::new(12.0)
}

fn get_ewald() -> SharedEwald {
    SharedEwald::new(Ewald::new(9.5, 7, None))
}

fn ewald_energy_computation(c: &mut Criterion) {
    let system = utils::get_system("nacl");
    let ewald = get_ewald();
    c.bench_function("nacl::ewald::energy", move |b| b.iter(|| {
        let _ = ewald.energy(&system);
    }));

    let system = utils::get_system("nacl");
    let ewald = get_ewald();
    c.bench_function("nacl::ewald::force", move |b| b.iter_batched_ref(
        || vec![Vector3D::zero(); system.size()],
        |forces| ewald.forces(&system, forces),
        BatchSize::SmallInput
    ));

    let system = utils::get_system("nacl");
    let ewald = get_ewald();
    c.bench_function("nacl::ewald::atomic_virial", move |b| b.iter(|| {
        let _ = ewald.atomic_virial(&system);
    }));

    let system = utils::get_system("nacl");
    let ewald = get_ewald();
    c.bench_function("nacl::ewald::molecular_virial", move |b| b.iter(|| {
        let _ = ewald.molecular_virial(&system);
    }));
}

fn ewald_monte_carlo_cache(c: &mut Criterion) {
    let mut system = utils::get_system("nacl");
    system.set_coulomb_potential(Box::new(get_ewald()));
    let mut cache = EnergyCache::new();
    cache.init(&system);

    c.bench_function("nacl::ewald::move_molecule_cost", move |b| b.iter_batched(
        || utils::move_rigid_molecule(&system),
        |(molid, positions)| cache.move_molecule_cost(&system, molid, &positions),
        BatchSize::SmallInput
    ));

    let mut system = utils::get_system("nacl");
    system.set_coulomb_potential(Box::new(get_ewald()));
    let mut cache = EnergyCache::new();
    cache.init(&system);

    c.bench_function("nacl::ewald::move_all_molecules_cost", move |b| b.iter_batched_ref(
        || utils::move_all_rigid_molecule(system.clone()),
        |system| cache.move_all_molecules_cost(system),
        BatchSize::SmallInput
    ));
}

fn wolf_energy_computation(c: &mut Criterion) {
    let system = utils::get_system("nacl");
    let wolf = get_wolf();
    c.bench_function("nacl::wolf::energy", move |b| b.iter(|| {
        let _ = wolf.energy(&system);
    }));

    let system = utils::get_system("nacl");
    let wolf = get_wolf();
    c.bench_function("nacl::wolf::force", move |b| b.iter_batched_ref(
        || vec![Vector3D::zero(); system.size()],
        |forces| wolf.forces(&system, forces),
        BatchSize::SmallInput
    ));

    let system = utils::get_system("nacl");
    let wolf = get_wolf();
    c.bench_function("nacl::wolf::atomic_virial", move |b| b.iter(|| {
        let _ = wolf.atomic_virial(&system);
    }));

    let system = utils::get_system("nacl");
    let wolf = get_wolf();
    c.bench_function("nacl::wolf::molecular_virial", move |b| b.iter(|| {
        let _ = wolf.molecular_virial(&system);
    }));
}

fn wolf_monte_carlo_cache(c: &mut Criterion) {
    let mut system = utils::get_system("nacl");
    system.set_coulomb_potential(Box::new(get_wolf()));
    let mut cache = EnergyCache::new();
    cache.init(&system);

    c.bench_function("nacl::wolf::move_molecule_cost", move |b| b.iter_batched(
        || utils::move_rigid_molecule(&system),
        |(molid, positions)| cache.move_molecule_cost(&system, molid, &positions),
        BatchSize::SmallInput
    ));

    let mut system = utils::get_system("nacl");
    system.set_coulomb_potential(Box::new(get_wolf()));
    let mut cache = EnergyCache::new();
    cache.init(&system);

    c.bench_function("nacl::wolf::move_all_molecules_cost", move |b| b.iter_batched_ref(
        || utils::move_all_rigid_molecule(system.clone()),
        |system| cache.move_all_molecules_cost(system),
        BatchSize::SmallInput
    ));
}

criterion_group!(ewald, ewald_energy_computation, ewald_monte_carlo_cache);
criterion_group!(wolf, wolf_energy_computation, wolf_monte_carlo_cache);

criterion_main!(ewald, wolf);
