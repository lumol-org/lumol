// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

//! Testing energy computation for SPC/E water using data from
//! https://www.nist.gov/mml/csd/chemical-informatics-research-group/spce-water-reference-calculations-9%C3%A5-cutoff
//! https://www.nist.gov/mml/csd/chemical-informatics-research-group/spce-water-reference-calculations-10å-cutoff

extern crate lumol;
extern crate lumol_input as input;

use lumol::sys::{System, UnitCell};
use lumol::sys::Trajectory;
use lumol::energy::{PairInteraction, LennardJones, NullPotential};
use lumol::energy::{Ewald, SharedEwald, PairRestriction, CoulombicPotential};
use lumol::consts::K_BOLTZMANN;

use std::path::Path;
use std::fs::File;
use std::io::prelude::*;

pub fn get_system(path: &str, cutoff: f64) -> System {
    let path = Path::new(file!()).parent().unwrap()
                                 .join("data")
                                 .join("nist-spce")
                                 .join(path);
    let mut system = Trajectory::open(&path)
                                .and_then(|mut traj| traj.read())
                                .unwrap();

    let mut file = File::open(path).unwrap();
    let mut buffer = String::new();
    file.read_to_string(&mut buffer).unwrap();
    let line = buffer.lines().skip(1).next().unwrap();
    let mut splited = line.split_whitespace();
    assert_eq!(splited.next(), Some("cell:"));
    let a: f64 = splited.next().expect("Missing 'a' cell parameter")
                        .parse().expect("'a' cell parameter is not a float");
    let b: f64 = splited.next().expect("Missing 'b' cell parameter")
                        .parse().expect("'b' cell parameter is not a float");
    let c: f64 = splited.next().expect("Missing 'c' cell parameter")
                        .parse().expect("'c' cell parameter is not a float");

    system.cell = UnitCell::ortho(a, b, c);

    for i in 0..system.size() {
        if i % 3 == 0 {
            system.add_bond(i, i + 1);
            system.add_bond(i, i + 2);
        }
    }

    for particle in system.particles_mut() {
        particle.charge = match particle.name() {
            "H" => 0.42380,
            "O" => -2.0 * 0.42380,
            other => panic!("Unknown particle name: {}", other)
        }
    }

    let mut lj = PairInteraction::new(Box::new(LennardJones{
        epsilon: 78.19743111 * K_BOLTZMANN,
        sigma: 3.16555789
    }), cutoff);
    lj.enable_tail_corrections();
    system.add_pair_potential("O", "O", lj);

    system.add_pair_potential("O", "H",
        PairInteraction::new(Box::new(NullPotential), cutoff)
    );
    system.add_pair_potential("H", "H",
        PairInteraction::new(Box::new(NullPotential), cutoff)
    );

    let mut ewald = Ewald::new(cutoff, 5);
    ewald.set_alpha(5.6 / f64::min(f64::min(a, b), c));
    let mut ewald = SharedEwald::new(ewald);
    ewald.set_restriction(PairRestriction::InterMolecular);
    system.set_coulomb_potential(Box::new(ewald));

    return system;
}

mod cutoff_9 {
    use super::*;
    use lumol::consts::K_BOLTZMANN;

    #[test]
    fn nist1() {
        let system = get_system("spce-1.xyz", 9.0);

        let energy = system.potential_energy() / K_BOLTZMANN;
        let expected = -4.88608e5;
        assert!(f64::abs((energy - expected) / expected) < 1e-3);
    }

    #[test]
    fn nist2() {
        let system = get_system("spce-2.xyz", 9.0);

        let energy = system.potential_energy() / K_BOLTZMANN;
        let expected = -1.06602e6;
        assert!(f64::abs((energy - expected) / expected) < 1e-3);
    }

    #[test]
    fn nist3() {
        let system = get_system("spce-3.xyz", 9.0);

        let energy = system.potential_energy() / K_BOLTZMANN;
        let expected = -1.71488e6;
        assert!(f64::abs((energy - expected) / expected) < 1e-3);
    }

    #[test]
    fn nist4() {
        let system = get_system("spce-4.xyz", 9.0);

        let energy = system.potential_energy() / K_BOLTZMANN;
        let expected = -3.08010e6;
        assert!(f64::abs((energy - expected) / expected) < 1e-3);
    }
}

mod cutoff_10 {
    use super::*;
    use lumol::consts::K_BOLTZMANN;

    #[test]
    fn nist1() {
        let system = get_system("spce-1.xyz", 10.0);

        let energy = system.potential_energy() / K_BOLTZMANN;
        let expected = -4.88604e5;
        assert!(f64::abs((energy - expected) / expected) < 1e-3);
    }

    #[test]
    fn nist2() {
        let system = get_system("spce-2.xyz", 10.0);

        let energy = system.potential_energy() / K_BOLTZMANN;
        let expected = -1.06590e6;
        assert!(f64::abs((energy - expected) / expected) < 1e-3);
    }

    #[test]
    fn nist3() {
        let system = get_system("spce-3.xyz", 10.0);

        let energy = system.potential_energy() / K_BOLTZMANN;
        let expected = -1.71488e6;
        assert!(f64::abs((energy - expected) / expected) < 1e-3);
    }

    #[test]
    fn nist4() {
        let system = get_system("spce-4.xyz", 10.0);

        let energy = system.potential_energy() / K_BOLTZMANN;
        let expected = -3.20501e6;
        assert!(f64::abs((energy - expected) / expected) < 1e-3);
    }
}
