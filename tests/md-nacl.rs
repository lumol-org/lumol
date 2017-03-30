// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 Lumol's contributors — BSD license

//! Testing physical properties of a NaCl crystal
extern crate lumol;
extern crate lumol_input as input;
extern crate env_logger;

use std::sync::{Once, ONCE_INIT};
pub static START: Once = ONCE_INIT;

mod utils;

mod wolf {
    use START;
    use lumol::units;
    use input::Input;
    use std::path::Path;

    #[test]
    fn constant_energy() {
        START.call_once(|| {::env_logger::init().unwrap();});
        let path = Path::new(file!()).parent().unwrap()
                                     .join("data")
                                     .join("md-nacl")
                                     .join("nve-wolf-small.toml");
        let mut config = Input::new(path).unwrap().read().unwrap();

        let e_initial = config.system.total_energy();
        config.simulation.run(&mut config.system, config.nsteps);

        let e_final = config.system.total_energy();
        assert!(f64::abs((e_initial - e_final)/e_final) < 1e-6);
    }

    #[test]
    fn anisotropic_berendsen() {
        START.call_once(|| {::env_logger::init().unwrap();});
        let path = Path::new(file!()).parent().unwrap()
                                     .join("data")
                                     .join("md-nacl")
                                     .join("npt-wolf-small.toml");
        let mut config = Input::new(path).unwrap().read().unwrap();

        let collecter = ::utils::Collecter::new(9000);
        let temperatures = collecter.temperatures();
        let pressures = collecter.pressures();

        config.simulation.add_output(Box::new(collecter));
        config.simulation.run(&mut config.system, config.nsteps);

        let expected = units::from(50000.0, "bar").unwrap();
        let pressure = ::utils::mean(pressures.clone());
        assert!(f64::abs(pressure - expected) / expected < 1e-3);

        let expected = units::from(273.0, "K").unwrap();
        let temperature = ::utils::mean(temperatures.clone());
        assert!(f64::abs(temperature - expected) / expected < 1e-2);
    }
}

mod ewald {
    use START;
    use lumol::units;
    use input::Input;
    use std::path::Path;

    #[test]
    fn constant_energy() {
        START.call_once(|| {::env_logger::init().unwrap();});
        let path = Path::new(file!()).parent().unwrap()
                                     .join("data")
                                     .join("md-nacl")
                                     .join("nve-ewald-small.toml");
        let mut config = Input::new(path).unwrap().read().unwrap();

        let e_initial = config.system.total_energy();
        config.simulation.run(&mut config.system, config.nsteps);
        let e_final = config.system.total_energy();
        assert!(f64::abs((e_initial - e_final)/e_final) < 5e-3);
    }

    #[test]
    fn energy() {
        START.call_once(|| {::env_logger::init().unwrap();});
        let path = Path::new(file!()).parent().unwrap()
                                     .join("data")
                                     .join("md-nacl")
                                     .join("energy-ewald-big.toml");
        let system = Input::new(path).unwrap().read_system().unwrap();
        let energy = units::to(system.total_energy(), "kcal/mol").unwrap();

        // Energy of this system given by LAMMPS in kcal/mol
        const LAMMPS_ENERGY: f64 = -48610.136;
        assert!(f64::abs((energy - LAMMPS_ENERGY)/LAMMPS_ENERGY) < 1e-3);
    }
}
