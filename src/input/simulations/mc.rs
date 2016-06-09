// Cymbalum, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license
use toml::Table;
use input::error::{Error, Result};
use input::FromToml;
use input::read_molecule;
use system::moltype;

use simulation::mc::*;
use units;

impl FromToml for MonteCarlo {
    fn from_toml(config: &Table) -> Result<MonteCarlo> {
        let temperature = extract_str!("temperature", config as "Monte-Carlo propagator");
        let temperature = try!(units::from_str(temperature));

        let mut mc = MonteCarlo::new(temperature);

        let moves = extract_slice!("moves", config as "Monte-Carlo propagator");
        for mc_move in moves {
            let mc_move = try!(mc_move.as_table().ok_or(
                Error::from("All moves must be tables in Monte-Carlo")
            ));

            let frequency = if mc_move.get("frequency").is_some() {
                extract_number!("frequency", mc_move as "Monte-Carlo move")
            } else {
                1.0
            };

            let mc_move: Box<MCMove> = match extract_type!(mc_move) {
                "Translate" => Box::new(try!(Translate::from_toml(mc_move))),
                "Rotate" => Box::new(try!(Rotate::from_toml(mc_move))),
                other => return Err(Error::from(
                    format!("Unknown Monte-Carlo move '{}'", other)
                ))
            };

            mc.add(mc_move, frequency);
        }
        return Ok(mc);
    }
}

impl FromToml for Translate {
    fn from_toml(config: &Table) -> Result<Translate> {
        let delta = extract_str!("delta", config as "Translate move");
        let delta = try!(units::from_str(delta));

        if config.get("molecule").is_some() {
            let molfile = extract_str!("molecule", config as "Translate move");
            let (molecule, atoms) = try!(read_molecule(molfile));
            let moltype = moltype(&molecule, &atoms);
            Ok(Translate::with_moltype(delta, moltype))
        } else {
            Ok(Translate::new(delta))
        }
    }
}

impl FromToml for Rotate {
    fn from_toml(config: &Table) -> Result<Rotate> {
        let delta = extract_str!("delta", config as "Rotate move");
        let delta = try!(units::from_str(delta));

        if config.get("molecule").is_some() {
            let molfile = extract_str!("molecule", config as "Translate move");
            let (molecule, atoms) = try!(read_molecule(molfile));
            let moltype = moltype(&molecule, &atoms);
            Ok(Rotate::with_moltype(delta, moltype))
        } else {
            Ok(Rotate::new(delta))
        }
    }
}

#[cfg(test)]
mod tests {
    use std::path::Path;
    use input::read_config;
    use input::testing::bad_inputs;

    #[test]
    fn mc() {
        let path = Path::new(file!()).parent().unwrap()
                                     .join("data")
                                     .join("mc.toml");
        assert!(read_config(&path).is_ok());
    }

    #[test]
    fn bad_mc() {
        for path in bad_inputs("simulations", "mc") {
            assert!(read_config(path).is_err());
        }
    }
}
