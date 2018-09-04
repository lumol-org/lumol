// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license
use std::path::PathBuf;
use toml::value::Table;

use lumol::sim::mc::*;
use lumol::sys::{molecule_type, read_molecule};
use lumol::units;

use FromTomlWithData;
use error::{Error, Result};
use extract;
use simulations::get_input_path;

impl FromTomlWithData for MonteCarlo {
    type Data = PathBuf;
    fn from_toml(config: &Table, root: PathBuf) -> Result<MonteCarlo> {
        let temperature = extract::str("temperature", config, "Monte Carlo propagator")?;
        let temperature = units::from_str(temperature)?;

        let mut mc = MonteCarlo::new(temperature);

        let has_update_frequency = config.get("update_frequency").is_some();
        if has_update_frequency {
            let update_frequency = extract::uint("update_frequency", config, "Monte Carlo propagator")?;
            mc.set_amplitude_update_frequency(update_frequency);
        }

        let moves = extract::slice("moves", config, "Monte Carlo propagator")?;
        for mc_move in moves {
            let mc_move = mc_move.as_table().ok_or(
                Error::from("All moves must be tables in Monte Carlo")
            )?;

            let frequency = if mc_move.get("frequency").is_some() {
                extract::number("frequency", mc_move, "Monte Carlo move")?
            } else {
                1.0
            };

            let target_acceptance = if mc_move.get("target_acceptance").is_some() {
                Some(extract::number("target_acceptance", mc_move, "Monte Carlo move")?)
            } else {
                None
            };

            let mc_move: Box<MCMove> = match extract::typ(mc_move, "Monte Carlo move")? {
                "Translate" => Box::new(Translate::from_toml(mc_move, root.clone())?),
                "Rotate" => Box::new(Rotate::from_toml(mc_move, root.clone())?),
                "Resize" => Box::new(Resize::from_toml(mc_move, root.clone())?),
                other => return Err(Error::from(format!("Unknown Monte Carlo move '{}'", other))),
            };

            match target_acceptance {
                Some(ta) => {
                    if !has_update_frequency {
                        return Err(Error::from(
                            "No 'update_frequency' found. Please specify \
                             'update_frequency' in combination with 'target_acceptance'",
                        ));
                    } else if ta < 0.0 || ta > 1.0 {
                        return Err(
                            Error::from("'target_acceptance' has to be between 0.0 and 1.0"),
                        );
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
        let delta = extract::str("delta", config, "Translate move")?;
        let delta = units::from_str(delta)?;

        if config.get("molecule").is_some() {
            let molfile = extract::str("molecule", config, "Translate move")?;
            let molfile = get_input_path(root, molfile);
            let (molecule, atoms) = read_molecule(molfile)?;
            let moltype = molecule_type(&molecule, atoms.as_slice());
            Ok(Translate::with_moltype(delta, moltype))
        } else {
            Ok(Translate::new(delta))
        }
    }
}

impl FromTomlWithData for Rotate {
    type Data = PathBuf;
    fn from_toml(config: &Table, root: PathBuf) -> Result<Rotate> {
        let delta = extract::str("delta", config, "Rotate move")?;
        let delta = units::from_str(delta)?;

        if config.get("molecule").is_some() {
            let molfile = extract::str("molecule", config, "Rotate move")?;
            let molfile = get_input_path(root, molfile);
            let (molecule, atoms) = read_molecule(molfile)?;
            let moltype = molecule_type(&molecule, atoms.as_slice());
            Ok(Rotate::with_moltype(delta, moltype))
        } else {
            Ok(Rotate::new(delta))
        }
    }
}

impl FromTomlWithData for Resize {
    type Data = PathBuf;
    fn from_toml(config: &Table, _: PathBuf) -> Result<Resize> {
        let pressure = extract::str("pressure", config, "Resize move")?;
        let pressure = units::from_str(pressure)?;

        let delta = extract::str("delta", config, "Resize move")?;
        let delta = units::from_str(delta)?;

        Ok(Resize::new(pressure, delta))
    }
}
