// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 Lumol's contributors — BSD license

#[macro_use]
extern crate bencher;
extern crate lumol;

use bencher::Bencher;

use lumol::Logger;
use lumol::energy::{Ewald, PairRestriction, CoulombicPotential, GlobalPotential};

use std::sync::{Once, ONCE_INIT};
pub static START: Once = ONCE_INIT;

mod utils;

fn get_ewald() -> Ewald {
    let mut ewald = Ewald::new(9.0, 7);
    ewald.set_restriction(PairRestriction::InterMolecular);
    ewald
}

fn energy_ewald(bencher: &mut Bencher) {
    START.call_once(|| {Logger::stdout();});
    let system = utils::get_system("water");
    let mut ewald = get_ewald();

    bencher.iter(||{
        let _ = ewald.energy(&system);
    })
}

fn forces_ewald(bencher: &mut Bencher) {
    START.call_once(|| {Logger::stdout();});
    let system = utils::get_system("water");
    let mut ewald = get_ewald();

    bencher.iter(||{
        let _ = ewald.forces(&system);
    })
}

fn virial_ewald(bencher: &mut Bencher) {
    START.call_once(|| {Logger::stdout();});
    let system = utils::get_system("water");
    let mut ewald = get_ewald();

    bencher.iter(||{
        let _ = ewald.virial(&system);
    })
}

benchmark_group!(ewald, energy_ewald, forces_ewald, virial_ewald);
benchmark_main!(ewald);
