extern crate lumol_tutorial_potential;
use lumol_tutorial_potential::Mie;

extern crate lumol;
extern crate lumol_input as input;
use input::Input;
use lumol::energy::PairInteraction;
use lumol::units;

// This function runs a small Monte-Carlo simulation for Argon using a Mie potential.
fn main() {
    // read configuration from input file located at the data/ folder
    let config = Input::new("data/argon.toml").unwrap().read().unwrap();
    let mut system = config.system;
    let mut simulation = config.simulation;

    // build a potential
    let mie = Mie::new(units::from(3.405, "A").unwrap(),
                       units::from(1.0, "kJ/mol").unwrap(),
                       12.0,
                       6.0);
    // use the potential with a cut off radius of 3 * sigma
    let mut pairs = PairInteraction::new(Box::new(mie), units::from(3.0 * 3.405, "A").unwrap());
    // use tail corrections to account for our truncation
    pairs.enable_tail_corrections();
    // finally add this interaction to the system for Argon atoms
    system.add_pair_potential("Ar", "Ar", pairs);

    // report the initial system energy
    let ener_init = units::to(system.total_energy(), "kJ/mol").unwrap();
    println!("initial energy     : {}", ener_init);

    // run the simulation
    simulation.run(&mut system, config.nsteps);

    // print some final information
    println!("final energy       : {}",
             units::to(system.total_energy(), "kJ/mol").unwrap());
    println!("It worked. Hooray!")
}
