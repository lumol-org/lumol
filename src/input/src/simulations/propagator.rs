// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license
use lumol::simulation::{Propagator, MolecularDynamics, MonteCarlo};

use error::{Error, Result};
use FromToml;
use extract;
use super::Input;

impl Input {
    /// Get the the simulation propagator. This is an internal function, public
    /// because of the code organization.
    // TODO: use restricted privacy here
    pub fn read_propagator(&self) -> Result<Box<Propagator>> {
        let config = try!(self.simulation_table());
        let propagator = try!(extract::table("propagator", config, "simulation"));
        match try!(extract::typ(propagator, "propagator")) {
            "MolecularDynamics" => Ok(Box::new(try!(
                MolecularDynamics::from_toml(propagator)
            ))),
            "MonteCarlo" => Ok(Box::new(try!(
                MonteCarlo::from_toml(propagator)
            ))),
            other => Err(Error::from(
                format!("Unknown propagator type '{}'", other)
            ))
        }
    }
}


#[cfg(test)]
mod tests {
    use Input;
    use testing::bad_inputs;

    #[test]
    fn bad_propagators() {
        for path in bad_inputs("simulations", "propagator") {
            assert!(Input::new(path).and_then(|input| input.read()).is_err());
        }
    }
}
