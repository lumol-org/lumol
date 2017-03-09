// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 Lumol's contributors — BSD license

#[macro_use]
extern crate bencher;
extern crate lumol;
extern crate lumol_input;

use bencher::Bencher;
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

benchmark_group!(argon, energy, forces, virial);
benchmark_main!(argon);
