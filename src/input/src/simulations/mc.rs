// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license
use toml::value::Table;
use std::path::PathBuf;

use lumol::sys::{read_molecule, molecule_type};
use lumol::sim::mc::*;
use lumol::units;

use error::{Error, Result};
use FromTomlWithData;
use extract;
use simulations::get_input_path;

impl FromTomlWithData for MonteCarlo {
    type Data = PathBuf;
    fn from_toml(config: &Table, root: PathBuf) -> Result<MonteCarlo> {
        let temperature = try!(extract::str("temperature", config, "Monte Carlo propagator"));
        let temperature = try!(units::from_str(temperature));

        let mut mc = MonteCarlo::new(temperature);

        let has_update_frequency = config.get("update_frequency").is_some();
        if has_update_frequency {
            let update_frequency =
                try!(extract::uint("update_frequency", config, "Monte Carlo propagator"));
            mc.set_amplitude_update_frequency(update_frequency);
        }

        let moves = try!(extract::slice("moves", config, "Monte Carlo propagator"));
        for mc_move in moves {
            let mc_move = try!(mc_move.as_table()
                .ok_or(Error::from("All moves must be tables in Monte Carlo")));

            let frequency = if mc_move.get("frequency").is_some() {
                try!(extract::number("frequency", mc_move, "Monte Carlo move"))
            } else {
                1.0
            };

            let target_acceptance = if mc_move.get("target_acceptance").is_some() {
                Some(try!(extract::number("target_acceptance", mc_move, "Monte Carlo move")))
            } else {
                None
            };

            let mc_move: Box<MCMove> = match try!(extract::typ(mc_move, "Monte Carlo move")) {
                "Translate" => Box::new(try!(Translate::from_toml(mc_move, root.clone()))),
                "Rotate" => Box::new(try!(Rotate::from_toml(mc_move, root.clone()))),
                "Resize" => Box::new(try!(Resize::from_toml(mc_move, root.clone()))),
                other => return Err(Error::from(format!("Unknown Monte Carlo move '{}'", other))),
            };

            match target_acceptance {
                Some(ta) => {
                    if !has_update_frequency {
                        return Err(Error::from(
                            "No 'update_frequency' found. Please specify \
                            'update_frequency' in combination with 'target_acceptance'"
                        ));
                    } else if ta < 0.0 || ta > 1.0 {
                        return Err(Error::from(
                            "'target_acceptance' has to be between 0.0 and 1.0"
                        ));
                    } else {
                        mc.add_move_with_acceptance(mc_move, frequency, ta)
                    }
                }
                None => mc.add(mc_move, frequency),
            }
        }
        return Ok(mc);
    }
}

impl FromTomlWithData for Translate {
    type Data = PathBuf;
    fn from_toml(config: &Table, root: PathBuf) -> Result<Translate> {
        let delta = try!(extract::str("delta", config, "Translate move"));
        let delta = try!(units::from_str(delta));

        if config.get("molecule").is_some() {
            let molfile = try!(extract::str("molecule", config, "Translate move"));
            let molfile = get_input_path(root, molfile);
            let (molecule, atoms) = try!(read_molecule(molfile));
            let moltype = molecule_type(&molecule, &atoms);
            Ok(Translate::with_moltype(delta, moltype))
        } else {
            Ok(Translate::new(delta))
        }
    }
}

impl FromTomlWithData for Rotate {
    type Data = PathBuf;
    fn from_toml(config: &Table, root: PathBuf) -> Result<Rotate> {
        let delta = try!(extract::str("delta", config, "Rotate move"));
        let delta = try!(units::from_str(delta));

        if config.get("molecule").is_some() {
            let molfile = try!(extract::str("molecule", config, "Rotate move"));
            let molfile = get_input_path(root, molfile);
            let (molecule, atoms) = try!(read_molecule(molfile));
            let moltype = molecule_type(&molecule, &atoms);
            Ok(Rotate::with_moltype(delta, moltype))
        } else {
            Ok(Rotate::new(delta))
        }
    }
}

impl FromTomlWithData for Resize {
    type Data = PathBuf;
    fn from_toml(config: &Table, _: PathBuf) -> Result<Resize> {
        let pressure = try!(extract::str("pressure", config, "Resize move"));
        let pressure = try!(units::from_str(pressure));

        let delta = try!(extract::str("delta", config, "Resize move"));
        let delta = try!(units::from_str(delta));

        Ok(Resize::new(pressure, delta))
    }
}
