// Cymbalum, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license
use toml::Table;
use input::error::{Error, Result};
use input::FromToml;

use simulation::{Propagator, MolecularDynamics};

pub fn read_propagator(config: &Table) -> Result<Box<Propagator>> {
    let propagator = extract_table!("propagator", config as "simulation");
    match extract_type!(propagator) {
        "MolecularDynamics" => Ok(Box::new(try!(MolecularDynamics::from_toml(propagator)))),
        other => Err(Error::from(
            format!("Unknown propagator type '{}'", other)
        ))
    }
}


#[cfg(test)]
mod tests {
    use input::read_config;
    use input::testing::bad_inputs;

    #[test]
    fn bad_propagators() {
        for path in bad_inputs("simulations", "propagator") {
            assert!(read_config(path).is_err());
        }
    }
}
