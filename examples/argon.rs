// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 Lumol's contributors — BSD license

//! Molecular dynamics simulation of an Argon crystal melt.
//!
//! In this example, we do everything by hand, from the system setup to the
//! simulation run.
extern crate lumol;

use lumol::sys::{System, UnitCell, Particle};
use lumol::sys::veloc::{BoltzmannVelocities, InitVelocities};
use lumol::energy::{LennardJones, PairInteraction};
use lumol::sim::{Simulation, MolecularDynamics};
use lumol::out::{TrajectoryOutput, EnergyOutput};
use lumol::types::Vector3D;
use lumol::units;

fn main() {
    let mut system = System::from_cell(UnitCell::cubic(17.0));

    // Create a cubic crystal of Argon by hand.
    for i in 0..5 {
        for j in 0..5 {
            for k in 0..5 {
                let mut part = Particle::new("Ar");
                part.position = Vector3D::new(
                        i as f64 * 3.4,
                        j as f64 * 3.4,
                        k as f64 * 3.4
                );
                system.add_particle(part);
            }
        }
    }

    let lj = Box::new(LennardJones {
        sigma: units::from(3.4, "A").unwrap(),
        epsilon: units::from(1.0, "kJ/mol").unwrap()
    });
    system.interactions_mut().add_pair("Ar", "Ar",
        PairInteraction::new(lj, units::from(8.5, "A").unwrap())
    );

    let mut velocities = BoltzmannVelocities::new(units::from(300.0, "K").unwrap());
    velocities.seed(129);
    velocities.init(&mut system);

    let mut simulation = Simulation::new(Box::new(
        MolecularDynamics::new(units::from(1.0, "fs").unwrap())
    ));

    let trajectory_out = Box::new(TrajectoryOutput::new("trajectory.xyz").unwrap());
    // Write the trajectory to `trajectory.xyz` every 10 steps
    simulation.add_output_with_frequency(trajectory_out, 10);

    let energy_out = Box::new(EnergyOutput::new("energy.dat").unwrap());
    // Write the energy to `energy.dat` every step
    simulation.add_output(energy_out);

    simulation.run(&mut system, 5000);

    println!("All done!")
}
