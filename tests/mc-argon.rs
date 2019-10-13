// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license
use lumol::input::Input;
use lumol::units;

use std::path::Path;
use std::sync::Once;
static START: Once = Once::new();

mod utils;

#[test]
fn npt() {
    START.call_once(::env_logger::init);
    let path = Path::new(file!()).parent()
                                 .unwrap()
                                 .join("data")
                                 .join("mc-argon")
                                 .join("npt.toml");

    let mut config = Input::new(path).unwrap().read().unwrap();

    let collecter = utils::Collecter::starting_at(500);
    let pressures = collecter.pressures();

    config.simulation.add_output(Box::new(collecter));
    config.simulation.run(&mut config.system, config.nsteps);

    let expected = units::from(200.0, "bar").unwrap();
    let pressure = crate::utils::mean(pressures.clone());
    assert!(f64::abs(pressure - expected) / expected < 1e-2);
}
