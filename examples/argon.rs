// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license
#![allow(clippy::cast_lossless)]

//! Molecular dynamics simulation of an Argon crystal melt.
//!
//! In this example, we do everything by hand, from the system setup to the
//! simulation run.
use lumol::{Particle, Molecule, System, UnitCell, Vector3D};
use lumol::energy::{LennardJones, PairInteraction};
use lumol::units;

use lumol::sim::output::{EnergyOutput, TrajectoryOutput};
use lumol::sim::{MolecularDynamics, Simulation};
use lumol::sim::{BoltzmannVelocities, InitVelocities};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut system = System::with_cell(UnitCell::cubic(17.0));

    // Create a cubic crystal of Argon by hand.
    for i in 0..5 {
        for j in 0..5 {
            for k in 0..5 {
                let position = Vector3D::new(i as f64 * 3.4, j as f64 * 3.4, k as f64 * 3.4);
                let particle = Particle::with_position("Ar", position);
                system.add_molecule(Molecule::new(particle));
            }
        }
    }

    let lj = Box::new(LennardJones {
        sigma: units::from(3.4, "A")?,
        epsilon: units::from(1.0, "kJ/mol")?,
    });
    system.set_pair_potential(
        ("Ar", "Ar"),
        PairInteraction::new(lj, units::from(8.5, "A")?),
    );

    let mut velocities = BoltzmannVelocities::new(units::from(300.0, "K")?);
    velocities.seed(129);
    velocities.init(&mut system);

    let md = MolecularDynamics::new(units::from(1.0, "fs")?);
    let mut simulation = Simulation::new(Box::new(md));

    let trajectory_out = Box::new(TrajectoryOutput::new("trajectory.xyz")?);
    // Write the trajectory to `trajectory.xyz` every 10 steps
    simulation.add_output_with_frequency(trajectory_out, 10);

    let energy_out = Box::new(EnergyOutput::new("energy.dat")?);
    // Write the energy to `energy.dat` every step
    simulation.add_output(energy_out);

    simulation.run(&mut system, 5000);

    println!("All done!");

    Ok(())
}
