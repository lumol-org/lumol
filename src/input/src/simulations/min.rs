// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license
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

        if let Some(criteria) = config.get("criteria") {
            let criteria = try!(criteria.as_table().ok_or(Error::from(
                "'criteria' must be a table in minimization propagator"
            )));
            let criteria = try!(EnergyCriteria::from_toml(criteria));
            Ok(Minimization::with_criteria(minimizer, criteria))
        } else {
            Ok(Minimization::new(minimizer))
        }
    }
}

impl FromToml for EnergyCriteria {
    fn from_toml(config: &Table) -> Result<EnergyCriteria> {
        let energy = try!(extract::str("energy", config, "minimization criteria"));
        let force2 = try!(extract::str("force2", config, "minimization criteria"));

        Ok(EnergyCriteria {
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
