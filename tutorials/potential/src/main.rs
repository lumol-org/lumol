use lumol_tutorial_potential::Mie;

use lumol::input::Input;
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
    let mut interaction = PairInteraction::new(Box::new(mie), units::from(3.0 * 3.405, "A").unwrap());
    // use tail corrections to account for our truncation
    interaction.enable_tail_corrections();
    // finally use this interaction for Argon atoms
    system.set_pair_potential(("Ar", "Ar"), interaction);

    // report the initial system energy
    let initial_energy = units::to(system.total_energy(), "kJ/mol").unwrap();
    println!("initial energy     : {}", initial_energy);

    // run the simulation
    simulation.run(&mut system, config.nsteps);

    // print some final information
    let final_energy = units::to(system.total_energy(), "kJ/mol").unwrap();
    println!("final energy       : {}", final_energy);
    println!("It worked. Hooray!")
}
