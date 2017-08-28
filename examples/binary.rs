// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

//! Monte Carlo simulation of a binary mixture of H20 and CO2.
extern crate lumol;
extern crate lumol_input as input;

use lumol::sys::{Molecule, Particle, ParticleVec, TrajectoryBuilder, UnitCell};
use lumol::sys::{read_molecule, molecule_type};
use lumol::sim::Simulation;
use lumol::sim::mc::{MonteCarlo, Translate, Rotate};
use lumol::units;

use input::InteractionsInput;

fn main() {
    let mut system = TrajectoryBuilder::new().open("data/binary.xyz")
                                             .and_then(|mut traj| traj.read())
                                             .unwrap();
    // Add bonds in the system
    for i in 0..system.molecules().len() / 3 {
        system.add_bond(3 * i,     3 * i + 1);
        system.add_bond(3 * i + 1, 3 * i + 2);
    }

    system.cell = UnitCell::cubic(25.0);
    let input = InteractionsInput::new("data/binary.toml").unwrap();
    input.read(&mut system).unwrap();

    let co2 = {
        // We can read files to get molecule type
        let (molecule, atoms) = read_molecule("data/CO2.xyz").unwrap();
        molecule_type(&molecule, atoms.as_slice())
    };
    let h2o = {
        // Or define a new molecule by hand
        let mut molecule = Molecule::new(0);
        molecule.merge_with(Molecule::new(1));
        molecule.merge_with(Molecule::new(2));

        molecule.add_bond(0, 1);
        molecule.add_bond(1, 2);

        let mut atoms = ParticleVec::new();
        atoms.push(Particle::new("H"));
        atoms.push(Particle::new("O"));
        atoms.push(Particle::new("H"));

        molecule_type(&molecule, atoms.as_slice())
    };

    let mut mc = MonteCarlo::new(units::from(500.0, "K").unwrap());

    // Use the molecular types of CO2 and H2O to specify different probabilities
    mc.add(Box::new(Translate::with_moltype(units::from(0.5, "A").unwrap(), co2)), 1.0);
    mc.add(Box::new(Rotate::with_moltype(units::from(10.0, "deg").unwrap(), co2)), 1.0);

    mc.add(Box::new(Translate::with_moltype(units::from(10.0, "A").unwrap(), h2o)), 2.0);
    mc.add(Box::new(Rotate::with_moltype(units::from(20.0, "deg").unwrap(), h2o)), 2.0);

    let mut simulation = Simulation::new(Box::new(mc));
    simulation.run(&mut system, 200_000_000);
}
