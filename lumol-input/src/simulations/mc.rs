// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license
#![allow(clippy::wildcard_imports)]

use std::path::PathBuf;
use toml::value::Table;

use lumol_sim::mc::*;
use lumol_core::read_molecule;
use lumol_core::units;

use crate::{Error, FromTomlWithData};
use crate::extract;
use crate::simulations::get_input_path;

impl FromTomlWithData for MonteCarlo {
    type Data = PathBuf;
    fn from_toml(config: &Table, root: PathBuf) -> Result<MonteCarlo, Error> {
        let temperature = extract::str("temperature", config, "Monte Carlo propagator")?;
        let temperature = units::from_str(temperature)?;
        let has_update_frequency = config.get("update_frequency").is_some();

        let mut builder = MonteCarloBuilder::new(temperature);
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

            let mc_move: Box<dyn MCMove> = match extract::typ(mc_move, "Monte Carlo move")? {
                "Translate" => Box::new(Translate::from_toml(mc_move, root.clone())?),
                "Rotate" => Box::new(Rotate::from_toml(mc_move, root.clone())?),
                "Resize" => Box::new(Resize::from_toml(mc_move, root.clone())?),
                other => return Err(Error::from(format!("unknown Monte Carlo move '{other}'"))),
            };

            match target_acceptance {
                Some(ta) => {
                    if !has_update_frequency {
                        return Err(Error::from(
                            "No 'update_frequency' found. Please specify 'update_frequency' in combination with 'target_acceptance'",
                        ));
                    } else if !(0.0..=1.0).contains(&ta) {
                        return Err(
                            Error::from("'target_acceptance' has to be between 0.0 and 1.0"),
                        );
                    }

                    builder.add(mc_move, frequency, ta);
                }
                None => builder.add(mc_move, frequency, None),
            }
        }

        let mut mc = builder.finish();
        if has_update_frequency {
            let update_frequency = extract::uint("update_frequency", config, "Monte Carlo propagator")?;
            mc.set_amplitude_update_frequency(update_frequency);
        }

        return Ok(mc);
    }
}

impl FromTomlWithData for Translate {
    type Data = PathBuf;
    fn from_toml(config: &Table, root: PathBuf) -> Result<Translate, Error> {
        let delta = extract::str("delta", config, "Translate move")?;
        let delta = units::from_str(delta)?;

        if config.get("molecule").is_some() {
            let molfile = extract::str("molecule", config, "Translate move")?;
            let molfile = get_input_path(root, molfile);
            let hash = read_molecule(molfile)?.as_ref().hash();
            Ok(Translate::new(delta, hash))
        } else {
            Ok(Translate::new(delta, None))
        }
    }
}

impl FromTomlWithData for Rotate {
    type Data = PathBuf;
    fn from_toml(config: &Table, root: PathBuf) -> Result<Rotate, Error> {
        let delta = extract::str("delta", config, "Rotate move")?;
        let delta = units::from_str(delta)?;

        if config.get("molecule").is_some() {
            let molfile = extract::str("molecule", config, "Rotate move")?;
            let molfile = get_input_path(root, molfile);
            let hash = read_molecule(molfile)?.as_ref().hash();
            Ok(Rotate::new(delta, hash))
        } else {
            Ok(Rotate::new(delta, None))
        }
    }
}

impl FromTomlWithData for Resize {
    type Data = PathBuf;
    fn from_toml(config: &Table, _: PathBuf) -> Result<Resize, Error> {
        let pressure = extract::str("pressure", config, "Resize move")?;
        let pressure = units::from_str(pressure)?;

        let delta = extract::str("delta", config, "Resize move")?;
        let delta = units::from_str(delta)?;

        Ok(Resize::new(pressure, delta))
    }
}
