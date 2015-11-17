/*
 * Cymbalum, Molecular Simulation in Rust
 * Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
*/
//! Testing physical properties of a NaCl crystal
extern crate cymbalum;
use self::cymbalum::*;

use std::path::Path;
use std::sync::{Once, ONCE_INIT};
pub static START: Once = ONCE_INIT;

pub fn setup_universe(potential: &str) -> Universe {
    let data_dir = Path::new(file!()).parent().unwrap();
    let configuration = data_dir.join("data").join("NaCl.xyz");
    let mut universe = Universe::from_file(configuration.to_str().unwrap()).unwrap();
    universe.set_cell(UnitCell::cubic(11.2804));

    let potentials = data_dir.join("data").join(potential);
    input::read_interactions(&mut universe, potentials).unwrap();

    let mut velocities = BoltzmanVelocities::new(units::from(300.0, "K").unwrap());
    velocities.init(&mut universe);
    return universe;
}

mod wolf {
    use super::*;
    use cymbalum::*;
    #[test]
    fn constant_energy() {
        START.call_once(|| {Logger::stdout();});
        let mut universe = setup_universe("NaCl-wolf.yml");
        let mut simulation = Simulation::new(
            MolecularDynamics::new(units::from(1.0, "fs").unwrap())
        );

        let e_initial = universe.total_energy();
        simulation.run(&mut universe, 1000);
        let e_final = universe.total_energy();
        assert!(f64::abs((e_initial - e_final)/e_final) < 1e-6);
    }

    #[test]
    fn anisotropic_berendsen() {
        START.call_once(|| {Logger::stdout();});
        let mut universe = setup_universe("NaCl-wolf.yml");
        let mut simulation = Simulation::new(
            MolecularDynamics::from_integrator(
                AnisoBerendsenBarostat::hydrostatic(
                    units::from(1.0, "fs").unwrap(),
                    units::from(5e4, "bar").unwrap()
                )
            )
        );

        simulation.run(&mut universe, 10000);
        let pressure = units::from(5e4, "bar").unwrap();
        assert!(f64::abs(universe.pressure() - pressure)/pressure < 1e-2);
    }
}

mod ewald {
    use super::*;
    use cymbalum::*;
    #[test]
    fn constant_energy() {
        START.call_once(|| {Logger::stdout();});
        let mut universe = setup_universe("NaCl-ewald.yml");
        let mut simulation = Simulation::new(
            MolecularDynamics::new(units::from(1.0, "fs").unwrap())
        );

        let e_initial = universe.total_energy();
        simulation.run(&mut universe, 100);
        let e_final = universe.total_energy();
        assert!(f64::abs((e_initial - e_final)/e_final) < 1e-4);
    }
}
