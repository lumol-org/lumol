// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

//! Testing physical properties of a sodium chloride crystal
use std::sync::Once;
static START: Once = Once::new();

mod utils;

mod wolf {
    use crate::START;
    use lumol::input::Input;
    use lumol::units;
    use std::path::Path;

    #[test]
    fn constant_energy() {
        START.call_once(::env_logger::init);
        let path = Path::new(file!()).parent()
                                     .unwrap()
                                     .join("data")
                                     .join("md-nacl")
                                     .join("nve-wolf-small.toml");
        let mut config = Input::new(path).unwrap().read().unwrap();

        let e_initial = config.system.total_energy();
        config.simulation.run(&mut config.system, config.nsteps);

        let e_final = config.system.total_energy();
        assert!(f64::abs((e_initial - e_final) / e_final) < 1e-4);
    }

    #[test]
    fn anisotropic_berendsen() {
        START.call_once(::env_logger::init);
        let path = Path::new(file!()).parent()
                                     .unwrap()
                                     .join("data")
                                     .join("md-nacl")
                                     .join("npt-wolf-small.toml");
        let mut config = Input::new(path).unwrap().read().unwrap();

        let collector = crate::utils::Collector::starting_at(9000);
        let temperatures = collector.temperatures();
        let pressures = collector.pressures();

        config.simulation.add_output(Box::new(collector));
        config.simulation.run(&mut config.system, config.nsteps);

        let expected = units::from(50000.0, "bar").unwrap();
        let pressure = crate::utils::mean(pressures);
        assert!(f64::abs(pressure - expected) / expected < 2e-3);

        let expected = units::from(273.0, "K").unwrap();
        let temperature = crate::utils::mean(temperatures);
        assert!(f64::abs(temperature - expected) / expected < 1e-2);
    }
}

mod ewald {
    use crate::START;
    use lumol::input::Input;
    use std::path::Path;

    #[test]
    fn constant_energy() {
        START.call_once(::env_logger::init);
        let path = Path::new(file!()).parent()
                                     .unwrap()
                                     .join("data")
                                     .join("md-nacl")
                                     .join("nve-ewald-small.toml");
        let mut config = Input::new(path).unwrap().read().unwrap();

        let e_initial = config.system.total_energy();
        config.simulation.run(&mut config.system, config.nsteps);
        let e_final = config.system.total_energy();
        assert!(f64::abs((e_initial - e_final) / e_final) < 5e-3);
    }

    #[test]
    fn constant_energy_kspace() {
        START.call_once(::env_logger::init);
        let path = Path::new(file!()).parent()
                                     .unwrap()
                                     .join("data")
                                     .join("md-nacl")
                                     .join("nve-ewald-kspace.toml");
        let mut config = Input::new(path).unwrap().read().unwrap();

        let e_initial = config.system.total_energy();
        config.simulation.run(&mut config.system, config.nsteps);
        let e_final = config.system.total_energy();
        assert!(f64::abs((e_initial - e_final) / e_final) < 5e-3);
    }
}
