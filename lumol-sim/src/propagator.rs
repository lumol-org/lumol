// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

//! A propagator is responsible for updating the system during a simulation
use lumol_core::{System, DegreesOfFreedom};

/// Possible temperature computation strategies. Different propagators needs
/// different ways to compute the temperature: Monte Carlo temperature is a
/// constant of the simulation, whereas for molecular dynamics we use the
/// instantaneous velocities.
#[derive(Clone, PartialEq, Debug)]
pub enum TemperatureStrategy {
    /// No specific strategy, use whatever strategy was already in use.
    None,
    /// Use the instantaneous velocities to compute the temperature
    Velocities,
    /// Use a fixed external temperature
    External(f64),
}

/// The propagator trait is the main algorithm of a simulation, i.e. the one
/// which update the system. The main function here is `propagate`, which
/// should propagate the simulation for one step.
pub trait Propagator {
    /// Get the temperature computation strategy for this propagator.
    ///
    /// This function is called once at thr beginning of the simulation
    fn temperature_strategy(&self) -> TemperatureStrategy;

    /// Get the number of degrees of freedom simulated by this propagator
    ///
    /// This function is called once at thr beginning of the simulation
    fn degrees_of_freedom(&self, system: &System) -> DegreesOfFreedom;

    /// Setup code, preparing all the meta-information needed about the
    /// simulation.
    fn setup(&mut self, _: &System) {}

    /// Propagate the system for one simulation step.
    fn propagate(&mut self, system: &mut System);

    /// Finish the simulation, and maybe output some information about it
    fn finish(&mut self, _: &System) {}
}
