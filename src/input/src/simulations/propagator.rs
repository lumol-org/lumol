// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license
use lumol::sim::{Propagator, MolecularDynamics, MonteCarlo};

use error::{Error, Result};
use {FromToml, FromTomlWithData};
use extract;
use super::Input;

impl Input {
    /// Get the the simulation propagator. This is an internal function, public
    /// because of the code organization.
    // TODO: use restricted privacy here
    #[doc(hidden)]
    pub fn read_propagator(&self) -> Result<Box<Propagator>> {
        let config = try!(self.simulation_table());
        let propagator = try!(extract::table("propagator", config, "simulation"));
        match try!(extract::typ(propagator, "propagator")) {
            "MolecularDynamics" => Ok(Box::new(try!(
                MolecularDynamics::from_toml(propagator)
            ))),
            "MonteCarlo" => Ok(Box::new(try!(
                MonteCarlo::from_toml(propagator, self.path.clone())
            ))),
            other => Err(Error::from(
                format!("Unknown propagator type '{}'", other)
            ))
        }
    }
}
