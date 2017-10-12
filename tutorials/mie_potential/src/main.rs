extern crate mie_potential;
use mie_potential::Mie;

extern crate lumol;
use lumol::sys::{System, UnitCell, Particle};
use lumol::energy::PairInteraction;
use lumol::sim::{Simulation, MonteCarlo};
use lumol::sim::mc::Translate;
use lumol::out::{EnergyOutput, PropertiesOutput};
use lumol::types::Vector3D;
use lumol::units;

// helper function to build a system made of 125 particles
// interacting via a Mie potential (using LJ exponents)
fn get_system() -> System {
    let mut system = System::with_cell(UnitCell::cubic(28.0));
    for i in 0..5 {
        for j in 0..5 {
            for k in 0..5 {
                let mut part = Particle::new("Ar");
                part.position = Vector3D::new(i as f64 * 3.4, j as f64 * 3.4, k as f64 * 3.4);
                system.add_particle(part);
            }
        }
    }
    println!("Added {} atoms to the system!", system.size());

    // build a potential
    let mie = Mie::new(units::from(3.405, "A").unwrap(),
                       units::from(1.0, "kJ/mol").unwrap(),
                       12.0,
                       6.0);
    // use the potential with a cut off radius of 3 * sigma
    let mut pairs = PairInteraction::new(Box::new(mie), 3.0 * 3.405);
    // use tail corrections to account for our truncation
    pairs.enable_tail_corrections();
    // finally add this interaction to the system for Argon atoms
    system.add_pair_potential("Ar", "Ar", pairs);
    system
}

// This function runs a small Monte-Carlo simulation for Argon using a Mie potential.
fn main() {
    let mut system = get_system();

    // Create the propagator
    let temperature = units::from(108.0, "K").unwrap();
    let mut mc = MonteCarlo::new(temperature);
    
    // Add moves to the propagator
    // delta_trans describes how far particles will be translated (initially)
    let delta_trans = units::from(0.15, "A").unwrap();
    // a move with an target acceptance will change delta_trans so that the desired acceptance is achieved
    mc.add_move_with_acceptance(Box::new(Translate::new(delta_trans)), 500.0, 0.5);
    mc.set_amplitude_update_frequency(1000);

    // create the simulation object using the propagator
    let mut simulation = Simulation::new(Box::new(mc));
    // add outputs
    simulation.add_output_with_frequency(
        Box::new(PropertiesOutput::new("mie_properties.dat").unwrap()), 125);
    simulation.add_output_with_frequency(Box::new(EnergyOutput::new("mie_energies.dat").unwrap()),
                                         125);
    
    // we have 125 particles so a "cycle" will be 125 particle displacements
    let steps = 125;
    let cycles = 1_000;
    println!("Simulating {} cycles with a total of {} steps",
             cycles,
             cycles * steps);

    // report the initial system energy
    let ener_init = units::to(system.total_energy(), "kJ/mol").unwrap();
    println!("initial energy     : {}", ener_init);

    // run the simulation
    simulation.run(&mut system, cycles * steps);
    
    // print some final information
    println!("final energy       : {}", units::to(system.total_energy(), "kJ/mol").unwrap());
    println!("It worked. Hooray!")
}