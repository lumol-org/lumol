// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

//! Molecular dynamics simulation of a crystal of sodium chloride, reading system and
//! potentials from files.
extern crate lumol;
extern crate lumol_input as input;

use lumol::sys::{UnitCell, TrajectoryBuilder};
use lumol::sys::veloc::{BoltzmannVelocities, InitVelocities};
use lumol::sim::Simulation;
use lumol::sim::md::{MolecularDynamics, RescaleThermostat};
use lumol::units;

use input::InteractionsInput;

fn main() {
    // Read the system fromt the `data/nacl.xyz` file
    let mut system = TrajectoryBuilder::new()
                                       .open("data/nacl.xyz")
                                       .and_then(|mut traj| traj.read())
                                       .unwrap();
    // Set the unit cell, as there is no unit cell data in XYZ files
    system.cell = UnitCell::cubic(units::from(22.5608, "A").unwrap());
    // Read the interactions from the `data/NaCl.toml` TOML file
    let input = InteractionsInput::new("data/NaCl.toml").unwrap();
    input.read(&mut system).unwrap();

    let mut velocities = BoltzmannVelocities::new(units::from(300.0, "K").unwrap());
    velocities.init(&mut system);

    let mut md = MolecularDynamics::new(units::from(1.0, "fs").unwrap());
    // Use a velocity rescaling thermostat
    md.set_thermostat(Box::new(
        RescaleThermostat::new(units::from(300.0, "K").unwrap())
    ));

    let mut simulation = Simulation::new(Box::new(md));
    simulation.run(&mut system, 1000);
}
