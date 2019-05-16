// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

//! Molecular dynamics simulation of a crystal of sodium chloride, reading system and
//! potentials from files.
use lumol::{TrajectoryBuilder, UnitCell};
use lumol::units;

use lumol::sim::Simulation;
use lumol::sim::md::{MolecularDynamics, RescaleThermostat};
use lumol::sim::{BoltzmannVelocities, InitVelocities};

use lumol::input::InteractionsInput;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Read the system fromt the `data/nacl.xyz` file
    let mut system = TrajectoryBuilder::new().open("data/nacl.xyz")?
                                             .read()?;
    // Set the unit cell, as there is no unit cell data in XYZ files
    system.cell = UnitCell::cubic(units::from(22.5608, "A")?);
    // Read the interactions from the `data/nacl.toml` TOML file
    let input = InteractionsInput::new("data/nacl.toml")?;
    input.read(&mut system)?;

    let mut velocities = BoltzmannVelocities::new(units::from(300.0, "K")?);
    velocities.init(&mut system);

    let mut md = MolecularDynamics::new(units::from(1.0, "fs")?);
    // Use a velocity rescaling thermostat
    md.set_thermostat(Box::new(RescaleThermostat::new(units::from(300.0, "K")?)));

    let mut simulation = Simulation::new(Box::new(md));
    simulation.run(&mut system, 1000);

    Ok(())
}
