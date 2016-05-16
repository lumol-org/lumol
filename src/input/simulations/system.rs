// Cymbalum, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license
use toml::{Table, Value};
use chemfiles;

use input::error::{Error, Result};
use input::Trajectory;
use input::read_interactions;
use input::interactions::read_interactions_toml;

use units;
use system::*;

pub fn read_system(config: &Table) -> Result<System> {
    let config = try!(system_table(config));

    let file = extract_str!("file", config as "system input");
    let mut trajectory = try!(Trajectory::open(file));

    if let Some(cell) = try!(read_cell(config)) {
        try!(trajectory.as_chemfiles().set_cell(cell));
    }

    if config.get("topology").is_some() {
        let topology = extract_str!("topology", config as "system input");
        try!(trajectory.as_chemfiles().set_topology_file(topology));
    }

    let guess_bonds = if let Some(guess_bonds) = config.get("guess_bonds") {
        try!(guess_bonds.as_bool().ok_or(
            Error::from("'guess_bonds' should be a boolean value in system")
        ))
    } else { false };

    let mut system = if guess_bonds {
        try!(trajectory.read_guess_bonds())
    } else {
        try!(trajectory.read())
    };

    try!(read_potentials(&mut system, config));

    if let Some(velocities) = config.get("velocities") {
        let velocities = try!(velocities.as_table().ok_or(
            Error::from("'velocities' must be a table in system input")
        ));
        try!(init_velocities(&mut system, velocities));
    }

    Ok(system)
}

fn system_table(config: &Table) -> Result<&Table> {
    let systems = extract_slice!("systems", config as "input file");
    if systems.len() != 1 {
        return Err(Error::from(
            "Only one system is supported in the input"
        ));
    }

    let system = try!(systems[0].as_table().ok_or(
        Error::from("Systems should be tables")
    ));

    return Ok(system);
}

fn read_cell(config: &Table) -> Result<Option<chemfiles::UnitCell>> {
    if let Some(cell) = config.get("cell") {
        match *cell {
            Value::Array(ref cell) => {
                if cell.len() == 3 {
                    let a = try!(get_cell_number(&cell[0]));
                    let b = try!(get_cell_number(&cell[1]));
                    let c = try!(get_cell_number(&cell[2]));

                    Ok(Some(try!(chemfiles::UnitCell::new(a, b, c))))
                } else if cell.len() == 6 {
                    let a = try!(get_cell_number(&cell[0]));
                    let b = try!(get_cell_number(&cell[1]));
                    let c = try!(get_cell_number(&cell[2]));
                    let alpha = try!(get_cell_number(&cell[3]));
                    let beta  = try!(get_cell_number(&cell[4]));
                    let gamma = try!(get_cell_number(&cell[5]));

                    Ok(Some(try!(chemfiles::UnitCell::triclinic(a, b, c, alpha, beta, gamma))))
                } else {
                    Err(Error::from("'cell' array must have a size of 3 or 6"))
                }
            },
            Value::Integer(cell) => {
                let cell = cell as f64;
                Ok(Some(try!(chemfiles::UnitCell::new(cell, cell, cell))))
            },
            Value::Float(cell) => {
                Ok(Some(try!(chemfiles::UnitCell::new(cell, cell, cell))))
            },
            _ => Err(Error::from("'cell' must be a number or an array"))
        }
    } else {
        Ok(None)
    }
}

fn get_cell_number(value: &Value) -> Result<f64> {
    if let Some(value) = value.as_integer() {
        Ok(value as f64)
    } else if let Some(value) = value.as_float() {
        Ok(value)
    } else {
        Err(Error::from("values must be numbers in 'cell' array"))
    }
}

fn init_velocities(system: &mut System, config: &Table) -> Result<()> {
    if config.get("init").is_some() {
        let temperature = extract_str!("init", config as "velocities initializer");
        let temperature = try!(units::from_str(temperature));
        let mut velocities = BoltzmanVelocities::new(temperature);
        velocities.init(system);
    } else {
        warn!("'velocities' key does nothing in this input file");
    }

    Ok(())
}

fn read_potentials(system: &mut System, config: &Table) -> Result<()> {
    if let Some(potentials) = config.get("potentials") {
        if let Some(potentials) = potentials.as_str() {
            try!(read_interactions(system, potentials));
        } else if let Some(potentials) = potentials.as_table() {
            try!(read_interactions_toml(system, potentials));
        } else {
            return Err(Error::from("'potentials' must be a string or a table"))
        }
    } else {
        warn!("No potentials found in input file");
    }
    Ok(())
}


#[cfg(test)]
mod tests {
    use std::path::Path;
    use input::read_config;
    use input::testing::bad_inputs;
    use system::{Particle, UnitCell};

    #[test]
    fn system() {
        let path = Path::new(file!()).parent().unwrap()
                                     .join("data")
                                     .join("system.toml");
        let config = read_config(&path).unwrap();
        assert_eq!(*config.system.cell(), UnitCell::cubic(20.0));
        assert_approx_eq!(config.system.temperature(), 300.0, 1e-12);

        let mut system = config.system;
        let last = system.size();
        system.add_particle(Particle::new("A"));
        system.add_particle(Particle::new("B"));
        assert_eq!(system.pair_potentials(last, last + 1).len(), 10);
    }

    #[test]
    fn inline_potentials() {
        let path = Path::new(file!()).parent().unwrap()
                                     .join("data")
                                     .join("potentials.toml");
        let config = read_config(&path).unwrap();
        assert_eq!(*config.system.cell(), UnitCell::cubic(20.0));
        assert_approx_eq!(config.system.temperature(), 300.0, 1e-12);
        assert_eq!(config.system.pair_potentials(0, 1).len(), 2);
    }

    #[test]
    fn bad_systems() {
        for path in bad_inputs("simulations", "system") {
            println!("{:?}",  path);
            assert!(read_config(path).is_err());
        }
    }
}
