// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

//! Monte Carlo simulation of a Xenon crystal melt.
extern crate lumol;

use lumol::energy::{LennardJones, PairInteraction};
use lumol::out::TrajectoryOutput;
use lumol::sim::Simulation;
use lumol::sim::mc::{MonteCarlo, Translate};
use lumol::sys::{TrajectoryBuilder, UnitCell};
use lumol::units;

fn main() -> Result<(), Box<std::error::Error>> {
    let mut system = TrajectoryBuilder::new().open("data/xenon.xyz")?
                                             .read()?;
    system.cell = UnitCell::cubic(units::from(21.65, "A")?);

    let lj = Box::new(LennardJones {
        sigma: units::from(4.57, "A")?,
        epsilon: units::from(1.87, "kJ/mol")?,
    });
    system.add_pair_potential(("Xe", "Xe"), PairInteraction::new(lj, 12.0));

    // Create a Monte Carlo propagator
    let mut mc = MonteCarlo::new(units::from(500.0, "K")?);
    // Add the `Translate` move with 0.5 A amplitude and 1.0 frequency
    mc.add(Box::new(Translate::new(units::from(0.5, "A")?, None)), 1.0);
    let mut simulation = Simulation::new(Box::new(mc));

    let trajectory_out = Box::new(TrajectoryOutput::new("trajectory.xyz")?);
    simulation.add_output_with_frequency(trajectory_out, 50);

    simulation.run(&mut system, 20_000);

    return Ok(());
}
