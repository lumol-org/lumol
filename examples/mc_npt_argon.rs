// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 Lumol's contributors — BSD license

//! Testing physical properties of a Lennard-Jones Argon using Monte-Carlo
//! simulation.
extern crate lumol;
extern crate lumol_input as input;

use lumol::sys::{System, Trajectory, UnitCell};
use lumol::energy::{LennardJones, PairInteraction};
use lumol::sim::Simulation;
use lumol::sim::mc::{MonteCarlo, Translate, Resize};
use lumol::units;
use lumol::out::{EnergyOutput, PropertiesOutput};

use std::path::Path;

fn get_system() -> System {
    let data_dir = Path::new(file!()).parent().unwrap();
    let configuration = data_dir.join("data").join("argon.xyz");
    let mut system = Trajectory::open(configuration)
                               .and_then(|mut traj| traj.read())
                               .unwrap();

    system.set_cell(UnitCell::cubic(31.0));

    // add LJ cut + tail corrections
    let lj = Box::new(LennardJones{
        sigma: units::from(3.405, "A").unwrap(),
        epsilon: units::from(1.0, "kJ/mol").unwrap()
    });

    let mut pairs = PairInteraction::new(lj, 10.215);
    pairs.enable_tail_corrections();
    system.interactions_mut().add_pair("Ar", "Ar", pairs);
    system
}


fn main() {
    let mut system = get_system();

    // To compare to NIST NVT data
    // State (T=0.85, rho=0.76)
    // rc = 10.215 A
    // T  = 102.231042155 K
    // V  = 25306.1891827 A^3
    // p  = 20.1586274874 bar
    // L  = 29.3590670058 A

    let temperature = units::from(102.231042155, "K").unwrap();
    let pressure = units::from(20.1586274874, "bar").unwrap();
    let mut mc = MonteCarlo::new(temperature);

    // Build move set for NPT
    let delta_trans = units::from(3.0, "A").unwrap();
    let delta_vol = units::from(0.05 * system.volume(), "A^3").unwrap();

    mc.add_move_with_acceptance(
        Box::new(Translate::new(delta_trans)), 500.0, 0.5);
    mc.add_move_with_acceptance(
        Box::new(Resize::new(pressure, delta_vol)), 2.0, 0.5);
    mc.set_amplitude_update_frequency(200);

    let mut simulation = Simulation::new(Box::new(mc));
    simulation.add_output_with_frequency(
        Box::new(PropertiesOutput::new("npt_argon_prp.dat").unwrap()), 500);
    simulation.add_output_with_frequency(
        Box::new(EnergyOutput::new("npt_argon_ener.dat").unwrap()), 500);

    println!("Starting simulation.");

    // run simulation
    simulation.run(&mut system, 2_000_000);
}
