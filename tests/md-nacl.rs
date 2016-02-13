// Cymbalum, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

//! Testing physical properties of a NaCl crystal
extern crate cymbalum;
use self::cymbalum::*;

use std::path::Path;
use std::sync::{Once, ONCE_INIT};
pub static START: Once = ONCE_INIT;

pub fn setup_system(potential: &str) -> System {
    let data_dir = Path::new(file!()).parent().unwrap();
    let configuration = data_dir.join("data").join("NaCl.xyz");
    let mut system = System::from_file(configuration.to_str().unwrap()).unwrap();
    system.set_cell(UnitCell::cubic(11.2804));

    let potentials = data_dir.join("data").join(potential);
    input::read_interactions(&mut system, potentials).unwrap();

    let mut velocities = BoltzmanVelocities::new(units::from(300.0, "K").unwrap());
    velocities.init(&mut system);
    return system;
}

mod wolf {
    use super::*;
    use cymbalum::*;
    #[test]
    fn constant_energy() {
        START.call_once(|| {Logger::stdout();});
        let mut system = setup_system("NaCl-wolf.yml");
        let mut simulation = Simulation::new(
            MolecularDynamics::new(units::from(1.0, "fs").unwrap())
        );

        let e_initial = system.total_energy();
        simulation.run(&mut system, 1000);
        let e_final = system.total_energy();
        assert!(f64::abs((e_initial - e_final)/e_final) < 1e-6);
    }

    #[test]
    fn anisotropic_berendsen() {
        START.call_once(|| {Logger::stdout();});
        let mut system = setup_system("NaCl-wolf.yml");
        let mut simulation = Simulation::new(
            MolecularDynamics::from_integrator(
                AnisoBerendsenBarostat::hydrostatic(
                    units::from(1.0, "fs").unwrap(),
                    units::from(5e4, "bar").unwrap()
                )
            )
        );

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
        let mut system = setup_system("NaCl-ewald.yml");
        let mut simulation = Simulation::new(
            MolecularDynamics::new(units::from(1.0, "fs").unwrap())
        );

        let e_initial = system.total_energy();
        simulation.run(&mut system, 100);
        let e_final = system.total_energy();
        assert!(f64::abs((e_initial - e_final)/e_final) < 1e-4);
    }
}
