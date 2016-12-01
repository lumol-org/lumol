// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

//! Testing physical properties of a Lennard-Jones Argon 
//! Monte-Carlo simulation
extern crate lumol;
extern crate lumol_input as input;

use lumol::Logger;
use lumol::sys::{System, Trajectory, UnitCell};
use lumol::energy::{LennardJones, PairInteraction};
use lumol::sim::Simulation;
use lumol::sim::mc::{MonteCarlo, Translate, Resize};
use lumol::units;

use std::sync::{Once, ONCE_INIT};
static START: Once = ONCE_INIT;

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
    // cut off: rc = 4.0 * sigma
    let mut pairs = PairInteraction::new(lj, 13.62);
    pairs.enable_tail_corrections();
    system.interactions_mut().add_pair("Ar", "Ar", pairs);
    return system;
}

#[test]
fn mc_npt_argon() {
    START.call_once(|| {Logger::stdout();});
    let mut system = get_system();
    
    // rng: set this here in case we want to change the seed
    //let mut rng = Box::new(rand::XorShiftRng::new_unseeded());
    //rng.reseed([2015u32, 42u32, 3u32, 12u32]);

    // set thermodynamic state (T,p) 
    // in LJ units, this is:
    // T* = 1.0
    // p* = 0.03
    // We use this state to compare to [1]
    
    let temperature = units::from(120.0, "K").unwrap();
    let pressure = units::from(12.6, "bar").unwrap();
    //let mut mc = MonteCarlo::from_rng(temperature, rng); // This is not working, why?
    let mut mc = MonteCarlo::new(temperature);

    // Build move set for NPT
    let delta_trans = units::from(3.0, "A").unwrap();
    let delta_vol = units::from(0.01 * system.volume(), "A^3").unwrap();

    // on average, for every particle move, perform two volume moves
    mc.add_move_with_acceptance(
        Box::new(Translate::new(delta_trans)), 500.0, 0.3);
    mc.add_move_with_acceptance(
        Box::new(Resize::new(pressure, delta_vol)), 2.0, 0.3);
    // if a move was called 500 times, perform an update
    mc.set_amplitude_update_frequency(500);
        
    let mut simulation = Simulation::new(Box::new(mc));
    
    // run simulation for at least 500 cycles (so that a volume update is performed)
    simulation.run(&mut system, 500000);

    // [1] Lotfi et al.: 
    // rho* = N/V*sigma^3 = 0.70179(37)
    // u*   = U/N/epsilon = -4.9018(25)
    // from NPT MD simulations with:
    // N = 1372
    // rc = 5.7*sigma (lower for higher densities) + tail corrections
    // 5000 timesteps equilibration
    // 55_000 timesteps production   
}
