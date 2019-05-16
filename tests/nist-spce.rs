// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

#![allow(clippy::unreadable_literal, clippy::needless_return)]

//! Testing energy computation for SPC/E water using data from
//! https://www.nist.gov/mml/csd/chemical-informatics-research-group/spce-water-reference-calculations-9å-cutoff
//! https://www.nist.gov/mml/csd/chemical-informatics-research-group/spce-water-reference-calculations-10å-cutoff
//!
//! Additionaly, forces are compared on these configuration with LAMMPS
//! calculations.
use lumol::consts::K_BOLTZMANN;
use lumol::energy::{CoulombicPotential, Ewald, PairRestriction, SharedEwald};
use lumol::energy::{LennardJones, NullPotential, PairInteraction};
use lumol::sys::{System, UnitCell};
use lumol::sys::TrajectoryBuilder;
use lumol::types::Vector3D;
use lumol::units;

use std::fs::File;
use std::io::prelude::*;
use std::path::Path;

pub fn get_system(path: &str) -> System {
    let path = Path::new(file!()).parent().unwrap().join("data").join("nist-spce").join(path);
    let mut system = TrajectoryBuilder::new().open(&path).and_then(|mut traj| traj.read()).unwrap();

    let mut file = File::open(path).unwrap();
    let mut buffer = String::new();
    file.read_to_string(&mut buffer).unwrap();
    let line = buffer.lines().nth(1).unwrap();
    let mut splited = line.split_whitespace();
    assert_eq!(splited.next(), Some("cell:"));
    let a: f64 = splited.next()
                        .expect("Missing 'a' cell parameter")
                        .parse()
                        .expect("'a' cell parameter is not a float");
    let b: f64 = splited.next()
                        .expect("Missing 'b' cell parameter")
                        .parse()
                        .expect("'b' cell parameter is not a float");
    let c: f64 = splited.next()
                        .expect("Missing 'c' cell parameter")
                        .parse()
                        .expect("'c' cell parameter is not a float");

    system.cell = UnitCell::ortho(a, b, c);

    for i in 0..system.size() {
        if i % 3 == 0 {
            system.add_bond(i, i + 1);
            system.add_bond(i, i + 2);
        }
    }

    for particle in system.particles_mut() {
        match particle.name.as_ref() {
            "H" => *particle.charge = 0.42380,
            "O" => *particle.charge = -2.0 * 0.42380,
            other => panic!("Unknown particle name: {}", other),
        }
    }

    return system;
}

pub fn set_nist_interactions(system: &mut System, cutoff: f64) {
    let mut lj = PairInteraction::new(
        Box::new(LennardJones {
            epsilon: 78.19743111 * K_BOLTZMANN,
            sigma: 3.16555789,
        }),
        cutoff,
    );
    lj.enable_tail_corrections();
    system.set_pair_potential(("O", "O"), lj);
    system.set_pair_potential(("O", "H"), PairInteraction::new(Box::new(NullPotential), cutoff));
    system.set_pair_potential(("H", "H"), PairInteraction::new(Box::new(NullPotential), cutoff));

    let alpha = 5.6 / f64::min(f64::min(system.cell.a(), system.cell.b()), system.cell.c());
    let mut ewald = SharedEwald::new(Ewald::new(cutoff, 5, alpha));
    ewald.set_restriction(PairRestriction::InterMolecular);
    system.set_coulomb_potential(Box::new(ewald));
}


// For comparaison of forces, as LAMMPS does not allow to set all parameters
// to the same one used by NIST
pub fn set_lammps_interactions(system: &mut System, cutoff: f64, kmax: usize, alpha: f64) {
    let lj = PairInteraction::new(
        Box::new(LennardJones {
            epsilon: 78.19743111 * K_BOLTZMANN,
            sigma: 3.16555789,
        }),
        cutoff,
    );
    system.set_pair_potential(("O", "O"), lj);
    system.set_pair_potential(("O", "H"), PairInteraction::new(Box::new(NullPotential), cutoff));
    system.set_pair_potential(("H", "H"), PairInteraction::new(Box::new(NullPotential), cutoff));

    let mut ewald = SharedEwald::new(Ewald::new(cutoff, kmax, alpha));
    ewald.set_restriction(PairRestriction::InterMolecular);
    system.set_coulomb_potential(Box::new(ewald));
}

pub fn get_forces(path: &str) -> Vec<Vector3D> {
    let path = Path::new(file!()).parent().unwrap().join("data").join("nist-spce").join(path);

    let mut file = File::open(path).unwrap();
    let mut buffer = String::new();
    file.read_to_string(&mut buffer).unwrap();
    let mut lines = buffer.lines();

    let natoms: usize = lines.next().unwrap().parse().unwrap();
    lines.next();

    let mut forces = vec![Vector3D::new(0.0, 0.0, 0.0); natoms];
    for (i, line) in lines.enumerate() {
        let mut splitted = line.split_whitespace();
        let num: usize = splitted.next().unwrap().parse().unwrap();
        assert_eq!(num, i + 1);

        forces[i][0] = splitted.next().unwrap().parse().unwrap();
        forces[i][1] = splitted.next().unwrap().parse().unwrap();
        forces[i][2] = splitted.next().unwrap().parse().unwrap();
    }

    return forces;
}

mod cutoff_9 {
    use super::*;
    use lumol::consts::K_BOLTZMANN;

    #[test]
    fn nist1() {
        let mut system = get_system("spce-1.xyz");
        set_nist_interactions(&mut system, 9.0);

        let energy = system.potential_energy() / K_BOLTZMANN;
        let expected = -4.88608e5;
        assert!(f64::abs((energy - expected) / expected) < 1e-3);

        let evaluator = system.energy_evaluator();
        let energy = evaluator.pairs() / K_BOLTZMANN;
        let expected = 9.98560e4;
        assert!(f64::abs((energy - expected) / expected) < 1e-3);

        let energy = evaluator.pairs_tail() / K_BOLTZMANN;
        let expected = -1.12959e3;
        assert!(f64::abs((energy - expected) / expected) < 1e-3);

        let energy = evaluator.coulomb() / K_BOLTZMANN;
        let expected = -5.87334e5;
        assert!(f64::abs((energy - expected) / expected) < 1e-3);
    }

    #[test]
    fn nist1_forces() {
        let mut system = get_system("spce-1.xyz");
        set_lammps_interactions(&mut system, 9.0, 8, 0.364209);

        let mut forces = system.forces();
        for force in &mut forces {
            force[0] = units::to(force[0], "kcal/mol/A").unwrap();
            force[1] = units::to(force[1], "kcal/mol/A").unwrap();
            force[2] = units::to(force[2], "kcal/mol/A").unwrap();
        }

        let expected = get_forces("forces-9-1.xyz");
        assert_eq!(forces.len(), expected.len());

        for (force, expected) in forces.iter().zip(&expected) {
            let delta = force - expected;
            assert!(f64::abs(delta[0] / expected[0]) < 5e-3);
            assert!(f64::abs(delta[1] / expected[1]) < 5e-3);
            assert!(f64::abs(delta[2] / expected[2]) < 5e-3);
        }
    }

    #[test]
    fn nist2() {
        let mut system = get_system("spce-2.xyz");
        set_nist_interactions(&mut system, 9.0);

        let energy = system.potential_energy() / K_BOLTZMANN;
        let expected = -1.06602e6;
        assert!(f64::abs((energy - expected) / expected) < 1e-3);

        let evaluator = system.energy_evaluator();
        let energy = evaluator.pairs() / K_BOLTZMANN;
        let expected = 1.94941e5;
        assert!(f64::abs((energy - expected) / expected) < 1e-3);

        let energy = evaluator.pairs_tail() / K_BOLTZMANN;
        let expected = -4.51836e3;
        assert!(f64::abs((energy - expected) / expected) < 1e-3);

        let energy = evaluator.coulomb() / K_BOLTZMANN;
        let expected = -1.25645e6;
        assert!(f64::abs((energy - expected) / expected) < 1e-3);
    }

    #[test]
    fn nist2_forces() {
        let mut system = get_system("spce-2.xyz");
        set_lammps_interactions(&mut system, 9.0, 8, 0.370036);

        let mut forces = system.forces();
        for force in &mut forces {
            force[0] = units::to(force[0], "kcal/mol/A").unwrap();
            force[1] = units::to(force[1], "kcal/mol/A").unwrap();
            force[2] = units::to(force[2], "kcal/mol/A").unwrap();
        }

        let expected = get_forces("forces-9-2.xyz");
        assert_eq!(forces.len(), expected.len());

        for (force, expected) in forces.iter().zip(&expected) {
            let delta = force - expected;
            assert!(f64::abs(delta[0] / expected[0]) < 5e-3);
            assert!(f64::abs(delta[1] / expected[1]) < 5e-3);
            assert!(f64::abs(delta[2] / expected[2]) < 5e-3);
        }
    }

    #[test]
    fn nist3() {
        let mut system = get_system("spce-3.xyz");
        set_nist_interactions(&mut system, 9.0);

        let energy = system.potential_energy() / K_BOLTZMANN;
        let expected = -1.71488e6;
        assert!(f64::abs((energy - expected) / expected) < 1e-3);

        let evaluator = system.energy_evaluator();
        let energy = evaluator.pairs() / K_BOLTZMANN;
        let expected = 3.57106e5;
        assert!(f64::abs((energy - expected) / expected) < 1e-3);

        let energy = evaluator.pairs_tail() / K_BOLTZMANN;
        let expected = -1.01663e4;
        assert!(f64::abs((energy - expected) / expected) < 1e-3);

        let energy = evaluator.coulomb() / K_BOLTZMANN;
        let expected = -2.06205e6;
        assert!(f64::abs((energy - expected) / expected) < 1e-3);
    }

    #[test]
    fn nist3_forces() {
        let mut system = get_system("spce-3.xyz");
        set_lammps_interactions(&mut system, 9.0, 8, 0.373403);

        let mut forces = system.forces();
        for force in &mut forces {
            force[0] = units::to(force[0], "kcal/mol/A").unwrap();
            force[1] = units::to(force[1], "kcal/mol/A").unwrap();
            force[2] = units::to(force[2], "kcal/mol/A").unwrap();
        }

        let expected = get_forces("forces-9-3.xyz");
        assert_eq!(forces.len(), expected.len());

        for (force, expected) in forces.iter().zip(&expected) {
            let delta = force - expected;
            for i in 0..3 {
                // Dynamic tolerance depending on the exact value of the force
                let tol = match f64::abs(expected[i]) {
                    a if a < 1e-1 => 1e-1,
                    b if b < 1.0 => 5e-2,
                    _ => 1e-2,
                };
                assert!(f64::abs(delta[i] / expected[i]) < tol);
            }
        }
    }

    #[test]
    fn nist4() {
        let mut system = get_system("spce-4.xyz");
        set_nist_interactions(&mut system, 9.0);

        let energy = system.potential_energy() / K_BOLTZMANN;
        let expected = -3.08010e6;
        assert!(f64::abs((energy - expected) / expected) < 1e-3);

        let evaluator = system.energy_evaluator();
        let energy = evaluator.pairs() / K_BOLTZMANN;
        let expected = 4.53536e5;
        assert!(f64::abs((energy - expected) / expected) < 1e-3);

        let energy = evaluator.pairs_tail() / K_BOLTZMANN;
        let expected = -1.88265e4;
        assert!(f64::abs((energy - expected) / expected) < 1e-3);

        let energy = evaluator.coulomb() / K_BOLTZMANN;
        let expected = -3.51481e6;
        assert!(f64::abs((energy - expected) / expected) < 1e-3);
    }

    #[test]
    fn nist4_forces() {
        let mut system = get_system("spce-4.xyz");
        set_lammps_interactions(&mut system, 9.0, 12, 0.370914);

        let mut forces = system.forces();
        for force in &mut forces {
            force[0] = units::to(force[0], "kcal/mol/A").unwrap();
            force[1] = units::to(force[1], "kcal/mol/A").unwrap();
            force[2] = units::to(force[2], "kcal/mol/A").unwrap();
        }

        let expected = get_forces("forces-9-4.xyz");
        assert_eq!(forces.len(), expected.len());

        for (force, expected) in forces.iter().zip(&expected) {
            let delta = force - expected;
            for i in 0..3 {
                // Dynamic tolerance depending on the exact value of the force
                let tol = match f64::abs(expected[i]) {
                    a if a < 1e-1 => 1e-1,
                    b if b < 1.0 => 5e-2,
                    _ => 1e-2,
                };
                assert!(f64::abs(delta[i] / expected[i]) < tol);
            }
        }
    }
}

mod cutoff_10 {
    use super::*;
    use lumol::consts::K_BOLTZMANN;

    #[test]
    fn nist1() {
        let mut system = get_system("spce-1.xyz");
        set_nist_interactions(&mut system, 10.0);

        let energy = system.potential_energy() / K_BOLTZMANN;
        let expected = -4.88604e5;
        assert!(f64::abs((energy - expected) / expected) < 1e-3);

        let evaluator = system.energy_evaluator();
        let energy = evaluator.pairs() / K_BOLTZMANN;
        let expected = 9.95387e4;
        assert!(f64::abs((energy - expected) / expected) < 1e-3);

        let energy = evaluator.pairs_tail() / K_BOLTZMANN;
        let expected = -8.23715e2;
        assert!(f64::abs((energy - expected) / expected) < 1e-3);

        let energy = evaluator.coulomb() / K_BOLTZMANN;
        let expected = -5.87319e5;
        assert!(f64::abs((energy - expected) / expected) < 1e-3);
    }

    #[test]
    fn nist1_forces() {
        let mut system = get_system("spce-1.xyz");
        set_lammps_interactions(&mut system, 10.0, 7, 0.326983);

        let mut forces = system.forces();
        for force in &mut forces {
            force[0] = units::to(force[0], "kcal/mol/A").unwrap();
            force[1] = units::to(force[1], "kcal/mol/A").unwrap();
            force[2] = units::to(force[2], "kcal/mol/A").unwrap();
        }

        let expected = get_forces("forces-10-1.xyz");
        assert_eq!(forces.len(), expected.len());

        for (force, expected) in forces.iter().zip(&expected) {
            let delta = force - expected;
            assert!(f64::abs(delta[0] / expected[0]) < 5e-3);
            assert!(f64::abs(delta[1] / expected[1]) < 5e-3);
            assert!(f64::abs(delta[2] / expected[2]) < 5e-3);
        }
    }

    #[test]
    fn nist2() {
        let mut system = get_system("spce-2.xyz");
        set_nist_interactions(&mut system, 10.0);

        let energy = system.potential_energy() / K_BOLTZMANN;
        let expected = -1.06590e6;
        assert!(f64::abs((energy - expected) / expected) < 1e-3);

        let evaluator = system.energy_evaluator();
        let energy = evaluator.pairs() / K_BOLTZMANN;
        let expected = 1.93712e5;
        assert!(f64::abs((energy - expected) / expected) < 1e-3);

        let energy = evaluator.pairs_tail() / K_BOLTZMANN;
        let expected = -3.29486e3;
        assert!(f64::abs((energy - expected) / expected) < 1e-3);

        let energy = evaluator.coulomb() / K_BOLTZMANN;
        let expected = -1.25632e6;
        assert!(f64::abs((energy - expected) / expected) < 1e-3);
    }

    #[test]
    fn nist2_forces() {
        let mut system = get_system("spce-2.xyz");
        set_lammps_interactions(&mut system, 10.0, 8, 0.332241);

        let mut forces = system.forces();
        for force in &mut forces {
            force[0] = units::to(force[0], "kcal/mol/A").unwrap();
            force[1] = units::to(force[1], "kcal/mol/A").unwrap();
            force[2] = units::to(force[2], "kcal/mol/A").unwrap();
        }

        let expected = get_forces("forces-10-2.xyz");
        assert_eq!(forces.len(), expected.len());

        for (force, expected) in forces.iter().zip(&expected) {
            let delta = force - expected;
            assert!(f64::abs(delta[0] / expected[0]) < 5e-3);
            assert!(f64::abs(delta[1] / expected[1]) < 5e-3);
            assert!(f64::abs(delta[2] / expected[2]) < 5e-3);
        }
    }

    #[test]
    fn nist3() {
        let mut system = get_system("spce-3.xyz");
        set_nist_interactions(&mut system, 10.0);

        let energy = system.potential_energy() / K_BOLTZMANN;
        let expected = -1.71488e6;
        assert!(f64::abs((energy - expected) / expected) < 1e-3);

        let evaluator = system.energy_evaluator();
        let energy = evaluator.pairs() / K_BOLTZMANN;
        let expected = 3.54344e5;
        assert!(f64::abs((energy - expected) / expected) < 1e-3);

        let energy = evaluator.pairs_tail() / K_BOLTZMANN;
        let expected = -7.41343e3;
        assert!(f64::abs((energy - expected) / expected) < 1e-3);

        let energy = evaluator.coulomb() / K_BOLTZMANN;
        let expected = -2.06182e6;
        assert!(f64::abs((energy - expected) / expected) < 1e-3);
    }

    #[test]
    fn nist3_forces() {
        let mut system = get_system("spce-3.xyz");
        set_lammps_interactions(&mut system, 10.0, 8, 0.335278);

        let mut forces = system.forces();
        for force in &mut forces {
            force[0] = units::to(force[0], "kcal/mol/A").unwrap();
            force[1] = units::to(force[1], "kcal/mol/A").unwrap();
            force[2] = units::to(force[2], "kcal/mol/A").unwrap();
        }

        let expected = get_forces("forces-10-3.xyz");
        assert_eq!(forces.len(), expected.len());

        for (force, expected) in forces.iter().zip(&expected) {
            let delta = force - expected;
            for i in 0..3 {
                // Dynamic tolerance depending on the exact value of the force
                let tol = match f64::abs(expected[i]) {
                    a if a < 1e-1 => 1e-1,
                    b if b < 1.0 => 5e-2,
                    _ => 1e-2,
                };
                assert!(f64::abs(delta[i] / expected[i]) < tol);
            }
        }
    }

    #[test]
    fn nist4() {
        let mut system = get_system("spce-4.xyz");
        set_nist_interactions(&mut system, 10.0);

        let energy = system.potential_energy() / K_BOLTZMANN;
        let expected = -3.20501e6;
        assert!(f64::abs((energy - expected) / expected) < 1e-3);

        let evaluator = system.energy_evaluator();
        let energy = evaluator.pairs() / K_BOLTZMANN;
        let expected = 4.48593e5;
        assert!(f64::abs((energy - expected) / expected) < 1e-3);

        let energy = evaluator.pairs_tail() / K_BOLTZMANN;
        let expected = -1.37286e4;
        assert!(f64::abs((energy - expected) / expected) < 1e-3);

        let energy = evaluator.coulomb() / K_BOLTZMANN;
        let expected = -3.63987e6;
        assert!(f64::abs((energy - expected) / expected) < 1e-3);
    }

    #[test]
    fn nist4_forces() {
        let mut system = get_system("spce-4.xyz");
        set_lammps_interactions(&mut system, 10.0, 11, 0.333033);

        let mut forces = system.forces();
        for force in &mut forces {
            force[0] = units::to(force[0], "kcal/mol/A").unwrap();
            force[1] = units::to(force[1], "kcal/mol/A").unwrap();
            force[2] = units::to(force[2], "kcal/mol/A").unwrap();
        }

        let expected = get_forces("forces-10-4.xyz");
        assert_eq!(forces.len(), expected.len());

        for (force, expected) in forces.iter().zip(&expected) {
            let delta = force - expected;
            for i in 0..3 {
                // Dynamic tolerance depending on the exact value of the force
                let tol = match f64::abs(expected[i]) {
                    a if a < 1e-1 => 1e-1,
                    b if b < 1.0 => 5e-2,
                    _ => 1e-2,
                };
                assert!(f64::abs(delta[i] / expected[i]) < tol);
            }
        }
    }
}
