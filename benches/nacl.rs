// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 Lumol's contributors — BSD license

#[macro_use]
extern crate bencher;
extern crate lumol;

use bencher::Bencher;

use lumol::Logger;
use lumol::sys::{System, Trajectory, UnitCell};
use lumol::energy::{Ewald, PairRestriction, CoulombicPotential, GlobalPotential};

use std::path::Path;
use std::sync::{Once, ONCE_INIT};
pub static START: Once = ONCE_INIT;

pub fn get_system() -> System {
    let data_dir = Path::new(file!()).parent().unwrap();
    let configuration = data_dir.join("data").join("nacl.xyz");
    let mut system = Trajectory::open(configuration)
                                .and_then(|mut traj| traj.read())
                                .unwrap();
    system.set_cell(UnitCell::cubic(11.2804));

    for particle in &mut system {
        let charge = match particle.name() {
            "Na" => 1.0,
            "Cl" => -1.0,
            _ => panic!("Missing charge value for {}", particle.name())
        };
        particle.charge = charge;
    }

    return system;
}

pub fn get_ewald() -> Ewald {
    let mut ewald = Ewald::new(5.5, 7);
    ewald.set_restriction(PairRestriction::InterMolecular);
    return ewald;
}

fn energy_ewald(bencher: &mut Bencher) {
    START.call_once(|| {Logger::stdout();});
    let system = get_system();
    let mut ewald = get_ewald();

    bencher.iter(||{
        let _ = ewald.energy(&system);
    })
}

fn forces_ewald(bencher: &mut Bencher) {
    START.call_once(|| {Logger::stdout();});
    let system = get_system();
    let mut ewald = get_ewald();

    bencher.iter(||{
        let _ = ewald.forces(&system);
    })
}

fn virial_ewald(bencher: &mut Bencher) {
    START.call_once(|| {Logger::stdout();});
    let system = get_system();
    let mut ewald = get_ewald();

    bencher.iter(||{
        let _ = ewald.virial(&system);
    })
}

benchmark_group!(ewald, energy_ewald, forces_ewald, virial_ewald);
benchmark_main!(ewald);
