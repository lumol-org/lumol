// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license
use toml::Table;

use lumol::chfl::read_molecule;
use lumol::system::moltype;
use lumol::simulation::mc::*;
use lumol::units;

use error::{Error, Result};
use FromToml;
use extract;

impl FromToml for MonteCarlo {
    fn from_toml(config: &Table) -> Result<MonteCarlo> {
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
        let delta = try!(extract::str("delta", config, "Translate move"));
        let delta = try!(units::from_str(delta));

        if config.get("molecule").is_some() {
            let molfile = try!(extract::str("molecule", config, "Translate move"));
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
        let delta = try!(extract::str("delta", config, "Rotate move"));
        let delta = try!(units::from_str(delta));

        if config.get("molecule").is_some() {
            let molfile = try!(extract::str("molecule", config, "Translate move"));
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
    use Input;
    use testing::bad_inputs;

    #[test]
    fn mc() {
        let path = Path::new(file!()).parent().unwrap()
                                     .join("data")
                                     .join("mc.toml");
        let input = Input::new(path).unwrap();
        assert!(input.read().is_ok());
    }

    #[test]
    fn bad_mc() {
        for path in bad_inputs("simulations", "mc") {
            assert!(Input::new(path).and_then(|input| input.read()).is_err());
        }
    }
}
