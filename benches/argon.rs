// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 Lumol's contributors — BSD license

#[macro_use]
extern crate bencher;
extern crate lumol;

use bencher::Bencher;

use lumol::units;
use lumol::energy::{PairInteraction, LennardJones};

mod utils;


fn energy(bencher: &mut Bencher) {
    let mut system = utils::get_system("argon");
    system.interactions_mut().add_pair("Ar", "Ar", PairInteraction::new(
        Box::new(LennardJones{
            epsilon: units::from(1.0, "kJ/mol").unwrap(),
            sigma: units::from(3.4, "A").unwrap(),
        }),
        units::from(10.0, "A").unwrap(),
    ));

    bencher.iter(||{
        let _ = system.potential_energy();
    })
}

fn forces(bencher: &mut Bencher) {
    let mut system = utils::get_system("argon");
    system.interactions_mut().add_pair("Ar", "Ar", PairInteraction::new(
        Box::new(LennardJones{
            epsilon: units::from(1.0, "kJ/mol").unwrap(),
            sigma: units::from(3.4, "A").unwrap(),
        }),
        units::from(10.0, "A").unwrap(),
    ));

    bencher.iter(||{
        let _ = system.forces();
    })
}

fn virial(bencher: &mut Bencher) {
    let mut system = utils::get_system("argon");
    system.interactions_mut().add_pair("Ar", "Ar", PairInteraction::new(
        Box::new(LennardJones{
            epsilon: units::from(1.0, "kJ/mol").unwrap(),
            sigma: units::from(3.4, "A").unwrap(),
        }),
        units::from(10.0, "A").unwrap(),
    ));

    bencher.iter(||{
        let _ = system.virial();
    })
}

benchmark_group!(argon, energy, forces, virial);
benchmark_main!(argon);
