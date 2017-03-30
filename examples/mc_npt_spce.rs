// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 Lumol's contributors — BSD license
extern crate lumol;
extern crate lumol_input as input;

use lumol::sys::{UnitCell, System, Trajectory};
use lumol::units;
use lumol::sim::Simulation;
use lumol::energy::{PairInteraction, PairRestriction};
use lumol::energy::{LennardJones, NullPotential};
use lumol::energy::{Ewald, SharedEwald};
use lumol::sim::mc::{MonteCarlo, Translate, Rotate, Resize};
use lumol::out::{EnergyOutput, PropertiesOutput, TrajectoryOutput};

use std::path::Path;

fn get_system() -> System {
    let data_dir = Path::new(file!()).parent().unwrap().join("data");
    let configuration = data_dir.join("spce.xyz");
    let mut system = Trajectory::open(configuration)
                                .and_then(|mut traj| traj.read_guess_bonds())
                                .unwrap();
    system.set_cell(UnitCell::cubic(20.0));

    // Add intermolecular interactions
    let lj = Box::new(LennardJones{
        sigma: units::from(3.16555789, "A").unwrap(),
        epsilon: units::from(0.65017254, "kJ/mol").unwrap()
    });
    let mut pairs = PairInteraction::new(lj, 6.0);
    // Restrict interactions to act only between different molecules.
    pairs.set_restriction(PairRestriction::InterMolecular);
    pairs.enable_tail_corrections();
    system.interactions_mut().add_pair("O", "O", pairs);

    // Add bonds: we use fixed bond lengths assuming
    // the equilibrium bond length. This means that both
    // energy as well as virial for the bond potential are
    // zero. Hence, we use a `NullPotential`.
    let bond = Box::new(NullPotential{});
    system.interactions_mut().add_bond("O", "H", bond);
    let angle = Box::new(NullPotential{});
    system.interactions_mut().add_angle("H", "O", "H", angle);

    // Set charges
    let h = 0.42380;
    let o = - 2.0 * h;
    for particle in &mut system {
        match particle.name() {
            "H" => particle.charge = h,
            "O" => particle.charge = o,
            _ => panic!("Unknown particle name in charge setting."),
        }
    }

    // Add ewald summation method.
    let ewald = Ewald::new(6.0, 7);
    system.interactions_mut().set_coulomb(Box::new(SharedEwald::new(ewald)));

    // Check if bonds are guessed correctly.
    assert_eq!(system.size() as f64 / system.molecules().len() as f64, 3.0);

    system
}

fn main() {
    let mut system = get_system();

    let temperature = units::from(300.0, "K").unwrap();
    let pressure = units::from(1.0, "bar").unwrap();

    let mut mc = MonteCarlo::new(temperature);

    // Build move set for NPT
    let delta_trans = units::from(50.0, "A").unwrap();
    let delta_rot = units::from(20.0, "deg").unwrap();
    let delta_vol = units::from(0.05 * system.volume(), "A^3").unwrap();

    // We strive for 50% acceptances - but that's arbitrary.
    // You should try different values to make your system run
    // efficiently.
    // Also, amplitudes are limited (180° for angles and the largest
    // cut-off in the system for translations). If you don't limit
    // amplitudes, they may explode.
    // I.e. for a dilute gas, translation amplitudes grow infinitely
    // since the system has so few particles that almost all moves are
    // accepted no matter what the amplitude will be.
    mc.add_move_with_acceptance(
        Box::new(Translate::new(delta_trans)), 45.0, 0.5);
    mc.add_move_with_acceptance(
        Box::new(Rotate::new(delta_rot)), 45.0, 0.5);
    mc.add_move_with_acceptance(
        Box::new(Resize::new(pressure, delta_vol)), 2.0, 0.5);
    mc.set_amplitude_update_frequency(500);

    // Setup simulation.
    let mut simulation = Simulation::new(Box::new(mc));

    // Add output.
    simulation.add_output_with_frequency(
        Box::new(PropertiesOutput::new("npt_spce_prp.dat").unwrap()), 500);
    simulation.add_output_with_frequency(
        Box::new(EnergyOutput::new("npt_spce_ener.dat").unwrap()), 500);
    simulation.add_output_with_frequency(
        Box::new(TrajectoryOutput::new("npt_spce_conf.xyz").unwrap()), 10000);

    // Often, a simulation is described using `MC cycles`.
    // We define a `cycle` to contain `nmols+2` moves.
    // (See frequencies of moves:
    // nmols/2 translations + nmols/2 rotations + 2 resize moves)
    let moves_per_cycle = 92;
    let cycles = 10000;

    // Some output and start of the simulation.
    println!("Simuation of 90 spce molecules.");
    let ener_init = units::to(system.total_energy(), "kJ/mol").unwrap();
    println!("Initial energy: {}", ener_init);
    println!("Simulating {} cycles with a total of {} moves",
        cycles, cycles * moves_per_cycle);
    println!("  running ....");
    simulation.run(&mut system, cycles * moves_per_cycle);
    println!("Done.");
}
