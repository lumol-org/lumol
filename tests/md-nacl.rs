// Cymbalum, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

//! Testing physical properties of a NaCl crystal
extern crate cymbalum;
use cymbalum::*;

use std::path::Path;
use std::sync::{Once, ONCE_INIT};
pub static START: Once = ONCE_INIT;

pub fn setup_system(potential: &str, file: &str) -> System {
    let data_dir = Path::new(file!()).parent().unwrap();

    let configuration = String::from("NaCl-") + file + ".xyz";
    let configuration = data_dir.join("data").join(configuration);
    let mut system = input::Trajectory::open(configuration)
                                        .and_then(|mut traj| traj.read())
                                        .unwrap();

    if file == "small" {
        system.set_cell(UnitCell::cubic(11.2804));
    } else if file == "big" {
        system.set_cell(UnitCell::cubic(22.5608));
    } else {
        unreachable!();
    }

    let potential = String::from("NaCl-") + potential + ".toml";
    let potential = data_dir.join("data").join(potential);
    input::read_interactions(&mut system, potential).unwrap();

    let mut velocities = BoltzmannVelocities::new(units::from(300.0, "K").unwrap());
    velocities.init(&mut system);
    return system;
}

mod wolf {
    use super::*;
    use cymbalum::*;

    #[test]
    fn constant_energy() {
        START.call_once(|| {Logger::stdout();});
        let mut system = setup_system("wolf", "small");
        let mut simulation = Simulation::new(Box::new(
            MolecularDynamics::new(units::from(1.0, "fs").unwrap())
        ));

        let e_initial = system.total_energy();
        simulation.run(&mut system, 1000);
        let e_final = system.total_energy();
        assert!(f64::abs((e_initial - e_final)/e_final) < 1e-6);
    }

    #[test]
    fn anisotropic_berendsen() {
        START.call_once(|| {Logger::stdout();});
        let mut system = setup_system("wolf", "small");
        let mut simulation = Simulation::new(Box::new(
            MolecularDynamics::from_integrator(Box::new(
                AnisoBerendsenBarostat::hydrostatic(
                    units::from(1.0, "fs").unwrap(),
                    units::from(5e4, "bar").unwrap(),
                    1000.0
                )
            ))
        ));

        simulation.run(&mut system, 10000);
        let pressure = units::from(5e4, "bar").unwrap();
        assert!(f64::abs(system.pressure() - pressure)/pressure < 1e-2);
    }
}

mod ewald {
    use super::*;
    use cymbalum::*;

    #[test]
    fn constant_energy() {
        START.call_once(|| {Logger::stdout();});
        let mut system = setup_system("ewald", "small");
        let mut simulation = Simulation::new(Box::new(
            MolecularDynamics::new(units::from(1.0, "fs").unwrap())
        ));

        let e_initial = system.total_energy();
        simulation.run(&mut system, 100);
        let e_final = system.total_energy();
        assert!(f64::abs((e_initial - e_final)/e_final) < 5e-3);
    }

    #[test]
    fn energy() {
        START.call_once(|| {Logger::stdout();});
        let system = setup_system("ewald", "big");
        let energy = units::to(system.total_energy(), "kcal/mol").unwrap();

        // Energy of this system given by LAMMPS in kcal/mol
        const LAMMPS_ENERGY: f64 = -48610.136;
        assert!(f64::abs((energy - LAMMPS_ENERGY)/LAMMPS_ENERGY) < 1e-3);
    }
}
