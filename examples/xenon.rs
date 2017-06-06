// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

//! Monte Carlo simulation of a Xenon crystal melt.
extern crate lumol;

use lumol::sys::{TrajectoryBuilder, UnitCell};
use lumol::energy::{PairInteraction, LennardJones};
use lumol::sim::Simulation;
use lumol::sim::mc::{MonteCarlo, Translate};
use lumol::out::TrajectoryOutput;
use lumol::units;

fn main() {
    let mut system = TrajectoryBuilder::new()
                                       .open("data/xenon.xyz")
                                       .and_then(|mut traj| traj.read())
                                       .unwrap();
    system.cell = UnitCell::cubic(units::from(21.65, "A").unwrap());

    let lj = Box::new(LennardJones{
        sigma: units::from(4.57, "A").unwrap(),
        epsilon: units::from(1.87, "kJ/mol").unwrap()
    });
    system.add_pair_potential("Xe", "Xe", PairInteraction::new(lj, 12.0));

    // Create a Monte Carlo propagator
    let mut mc = MonteCarlo::new(units::from(500.0, "K").unwrap());
    // Add the `Translate` move with 0.5 A amplitude and 1.0 frequency
    mc.add(
        Box::new(Translate::new(units::from(0.5, "A").unwrap())),
        1.0
    );
    let mut simulation = Simulation::new(Box::new(mc));

    let trajectory_out = Box::new(TrajectoryOutput::new("trajectory.xyz").unwrap());
    simulation.add_output_with_frequency(trajectory_out, 50);

    simulation.run(&mut system, 20000);
}
