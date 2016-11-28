// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license
use toml::Table;
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
        let temperature = try!(extract::str("temperature", config, "Monte-Carlo propagator"));
        let temperature = try!(units::from_str(temperature));

        let mut mc = MonteCarlo::new(temperature);

        let moves = try!(extract::slice("moves", config, "Monte-Carlo propagator"));
        for mc_move in moves {
            let mc_move = try!(mc_move.as_table().ok_or(
                Error::from("All moves must be tables in Monte-Carlo")
            ));

            let frequency = if mc_move.get("frequency").is_some() {
                try!(extract::number("frequency", mc_move, "Monte-Carlo move"))
            } else {
                1.0
            };

            let mc_move: Box<MCMove> = match try!(extract::typ(mc_move, "Monte-Carlo move")) {
                "Translate" => Box::new(try!(Translate::from_toml(mc_move, root.clone()))),
                "Rotate" => Box::new(try!(Rotate::from_toml(mc_move, root.clone()))),
                other => return Err(Error::from(
                    format!("Unknown Monte-Carlo move '{}'", other)
                ))
            };

            mc.add(mc_move, frequency);
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
            let molfile = try!(extract::str("molecule", config, "Translate move"));
            let molfile = get_input_path(root, molfile);
            let (molecule, atoms) = try!(read_molecule(molfile));
            let moltype = molecule_type(&molecule, &atoms);
            Ok(Rotate::with_moltype(delta, moltype))
        } else {
            Ok(Rotate::new(delta))
        }
    }
}
