// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license
use toml::value::Table;

use lumol::sim::min::*;
use lumol::units;

use FromToml;
use error::{Error, Result};
use extract;

impl FromToml for Minimization {
    fn from_toml(config: &Table) -> Result<Minimization> {
        let minimizer = extract::table("minimizer", config, "minimization propagator")?;

        let minimizer: Box<Minimizer> = match extract::typ(minimizer, "minimizer")? {
            "SteepestDescent" => Box::new(SteepestDescent::from_toml(minimizer)?),
            other => return Err(Error::from(format!("Unknown minimizer '{}'", other))),
        };

        if let Some(tolerance) = config.get("tolerance") {
            let tolerance = tolerance.as_table().ok_or(
                Error::from("'tolerance' must be a table in minimization propagator")
            )?;
            let tolerance = Tolerance::from_toml(tolerance)?;
            Ok(Minimization::with_tolerance(minimizer, tolerance))
        } else {
            Ok(Minimization::new(minimizer))
        }
    }
}

impl FromToml for Tolerance {
    fn from_toml(config: &Table) -> Result<Tolerance> {
        let energy = extract::str("energy", config, "minimization tolerance")?;
        let force2 = extract::str("force2", config, "minimization tolerance")?;

        Ok(Tolerance {
            energy: units::from_str(energy)?,
            force2: units::from_str(force2)?,
        })
    }
}


impl FromToml for SteepestDescent {
    fn from_toml(_: &Table) -> Result<SteepestDescent> {
        Ok(SteepestDescent::new())
    }
}
