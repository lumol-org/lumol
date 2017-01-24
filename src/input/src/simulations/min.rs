// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 Lumol's contributors — BSD license
use toml::Table;

use lumol::sim::min::*;
use lumol::units;

use error::{Error, Result};
use FromToml;
use extract;

impl FromToml for Minimization {
    fn from_toml(config: &Table) -> Result<Minimization> {
        let minimizer = try!(extract::table("minimizer", config, "minimization propagator"));

        let minimizer: Box<Minimizer> = match try!(extract::typ(minimizer, "minimizer")) {
            "SteepestDescent" => Box::new(try!(
                SteepestDescent::from_toml(minimizer)
            )),
            other => return Err(Error::from(
                format!("Unknown minimizer '{}'", other)
            ))
        };

        if let Some(tolerance) = config.get("tolerance") {
            let tolerance = try!(tolerance.as_table().ok_or(Error::from(
                "'tolerance' must be a table in minimization propagator"
            )));
            let tolerance = try!(Tolerance::from_toml(tolerance));
            Ok(Minimization::with_tolerance(minimizer, tolerance))
        } else {
            Ok(Minimization::new(minimizer))
        }
    }
}

impl FromToml for Tolerance {
    fn from_toml(config: &Table) -> Result<Tolerance> {
        let energy = try!(extract::str("energy", config, "minimization tolerance"));
        let force2 = try!(extract::str("force2", config, "minimization tolerance"));

        Ok(Tolerance {
            energy: try!(units::from_str(energy)),
            force2: try!(units::from_str(force2)),
        })

    }
}


impl FromToml for SteepestDescent {
    fn from_toml(_: &Table) -> Result<SteepestDescent> {
        Ok(SteepestDescent::new())
    }
}
