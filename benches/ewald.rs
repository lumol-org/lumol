// Cymbalum, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license
#![feature(test)]
extern crate test;
extern crate cymbalum;
use cymbalum::*;

use std::path::Path;
use std::sync::{Once, ONCE_INIT};
pub static START: Once = ONCE_INIT;

pub fn get_system(name: &str) -> System {
    let data_dir = Path::new(file!()).parent().unwrap();
    let configuration = data_dir.join("data").join(name);
    let mut system = io::Trajectory::open(configuration)
                                     .and_then(|mut traj| traj.read())
                                     .unwrap();
    let cell = if name == "NaCl.xyz" {
        UnitCell::cubic(11.2804)
    } else if name == "water.xyz" {
        UnitCell::cubic(18.0)
    } else { unreachable!() };
    system.set_cell(cell);

    // Use a smaller system to have faster benches and more iteration
    if name == "water.xyz" {
        for _ in 0..50 {
            let last = system.size() - 1;
            system.remove_molecule_containing(last);
        }
    }

    for p in &mut system {
        let charge = match p.name() {
            "O" => -0.82,
            "H" => 0.41,
            "Na" => 1.0,
            "Cl" => -1.0,
            _ => panic!("Missing charge value for {}", p.name())
        };
        p.charge = charge;
    }

    return system;
}

pub fn get_ewald() -> Ewald {
    let mut ewald = Ewald::new(5.5, 7);
    ewald.set_restriction(PairRestriction::InterMolecular);
    return ewald;
}

mod nacl {
    use super::*;
    use cymbalum::*;
    use test::Bencher;

    #[bench]
    fn energy(bencher: &mut Bencher) {
        START.call_once(|| {Logger::stdout();});
        let system = get_system("NaCl.xyz");
        let ewald = get_ewald();

        bencher.iter(||{
            let _ = ewald.energy(&system);
        })
    }

    #[bench]
    fn forces(bencher: &mut Bencher) {
        START.call_once(|| {Logger::stdout();});
        let system = get_system("NaCl.xyz");
        let ewald = get_ewald();

        bencher.iter(||{
            let _ = ewald.forces(&system);
        })
    }
}


mod water {
    use super::*;
    use cymbalum::*;
    use test::Bencher;

    #[bench]
    fn energy(bencher: &mut Bencher) {
        START.call_once(|| {Logger::stdout();});
        let system = get_system("water.xyz");
        let ewald = get_ewald();

        bencher.iter(||{
            let _ = ewald.energy(&system);
        })
    }

    #[bench]
    fn forces(bencher: &mut Bencher) {
        START.call_once(|| {Logger::stdout();});
        let system = get_system("water.xyz");
        let ewald = get_ewald();

        bencher.iter(||{
            let _ = ewald.forces(&system);
        })
    }
}
