// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 Lumol's contributors — BSD license

#[macro_use]
extern crate bencher;
extern crate lumol;

use bencher::Bencher;

use lumol::Logger;
use lumol::energy::{Ewald, Wolf, GlobalPotential};

use std::sync::{Once, ONCE_INIT};
pub static START: Once = ONCE_INIT;

mod utils;


fn energy_ewald(bencher: &mut Bencher) {
    START.call_once(|| {Logger::stdout();});
    let system = utils::get_system("nacl");
    let mut ewald = Ewald::new(9.5, 7);

    bencher.iter(||{
        let _ = ewald.energy(&system);
    })
}

fn forces_ewald(bencher: &mut Bencher) {
    START.call_once(|| {Logger::stdout();});
    let system = utils::get_system("nacl");
    let mut ewald = Ewald::new(9.5, 7);

    bencher.iter(||{
        let _ = ewald.forces(&system);
    })
}

fn virial_ewald(bencher: &mut Bencher) {
    START.call_once(|| {Logger::stdout();});
    let system = utils::get_system("nacl");
    let mut ewald = Ewald::new(9.5, 7);

    bencher.iter(||{
        let _ = ewald.virial(&system);
    })
}

fn energy_wolf(bencher: &mut Bencher) {
    START.call_once(|| {Logger::stdout();});
    let system = utils::get_system("nacl");
    let mut wolf = Wolf::new(12.0);

    bencher.iter(||{
        let _ = wolf.energy(&system);
    })
}

fn forces_wolf(bencher: &mut Bencher) {
    START.call_once(|| {Logger::stdout();});
    let system = utils::get_system("nacl");
    let mut wolf = Wolf::new(12.0);

    bencher.iter(||{
        let _ = wolf.forces(&system);
    })
}

fn virial_wolf(bencher: &mut Bencher) {
    START.call_once(|| {Logger::stdout();});
    let system = utils::get_system("nacl");
    let mut wolf = Wolf::new(12.0);

    bencher.iter(||{
        let _ = wolf.virial(&system);
    })
}

benchmark_group!(ewald, energy_ewald, forces_ewald, virial_ewald);
benchmark_group!(wolf, energy_wolf, forces_wolf, virial_wolf);
benchmark_main!(ewald, wolf);
