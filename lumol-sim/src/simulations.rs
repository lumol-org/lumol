// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license
use log::{info, warn};

use lumol_core::{System, DegreesOfFreedom, Vector3D};

use crate::output::Output;
use crate::propagator::{Propagator, TemperatureStrategy};

/// Writing an output at a given frequency
struct OutputFrequency {
    /// The output to use
    output: Box<dyn Output>,
    /// The frequency. `output` will be used every time the system step matches
    /// this frequency.
    frequency: u64,
}

impl OutputFrequency {
    pub fn new(output: Box<dyn Output>) -> OutputFrequency {
        OutputFrequency {
            frequency: 1,
            output: output,
        }
    }

    pub fn with_frequency(output: Box<dyn Output>, frequency: u64) -> OutputFrequency {
        OutputFrequency {
            frequency: frequency,
            output: output,
        }
    }
}

impl Output for OutputFrequency {
    fn setup(&mut self, system: &System) {
        self.output.setup(system);
    }

    fn write(&mut self, system: &System) {
        if system.step % self.frequency == 0 {
            self.output.write(system);
        }
    }

    fn finish(&mut self, system: &System) {
        self.output.finish(system);
    }
}

/// The Simulation struct holds all the needed algorithms for running the
/// simulation. It should be use together with a `System` to perform the
/// simulation.
pub struct Simulation {
    propagator: Box<dyn Propagator>,
    outputs: Vec<OutputFrequency>,
}

impl Simulation {
    /// Create a new Simulation from a Propagator.
    pub fn new(propagator: Box<dyn Propagator>) -> Simulation {
        Simulation {
            propagator: propagator,
            outputs: Vec::new(),
        }
    }

    /// Run the simulation on System for `nsteps` steps.
    pub fn run(&mut self, system: &mut System, nsteps: usize) {
        match self.propagator.temperature_strategy() {
            TemperatureStrategy::External(temperature) => {
                system.simulated_temperature(Some(temperature));
            }
            TemperatureStrategy::Velocities => system.simulated_temperature(None),
            TemperatureStrategy::None => {}
        }

        system.simulated_degrees_of_freedom = self.propagator.degrees_of_freedom(system);
        match system.simulated_degrees_of_freedom {
            DegreesOfFreedom::Molecules => info!("All molecules are treated as rigid bodies"),
            DegreesOfFreedom::Particles => info!("All particles are allowed to move in this simulation"),
            DegreesOfFreedom::Frozen(n) => info!("{} degrees of freedom are frozen", n),
        }

        system.check();
        self.setup(system);
        for i in 0..nsteps {
            self.propagator.propagate(system);
            system.step += 1;
            for output in &mut self.outputs {
                output.write(system);
            }

            if i % 10_000 == 0 {
                sanity_check(system);
            }
        }
        self.finish(system);
    }

    /// Add a new `Output` algorithm in the outputs list
    pub fn add_output(&mut self, output: Box<dyn Output>) {
        self.outputs.push(OutputFrequency::new(output));
    }

    /// Add a new `Output` algorithm in the outputs list, which will be used
    /// at the given frequency. The output will be used every time the system
    /// step matches this frequency.
    pub fn add_output_with_frequency(&mut self, output: Box<dyn Output>, frequency: u64) {
        self.outputs.push(OutputFrequency::with_frequency(output, frequency));
    }

    fn setup(&mut self, system: &mut System) {
        self.propagator.setup(system);
        for output in &mut self.outputs {
            output.setup(system);
        }
    }

    fn finish(&mut self, system: &mut System) {
        self.propagator.finish(system);
        for output in &mut self.outputs {
            output.finish(system);
        }
    }
}

/// Perform some sanity checks on the system
fn sanity_check(system: &System) {
    for position in system.particles().position {
        // The value of 1e6 A should be a good enough threshold. Even with
        // big boxes (100 A), and going through the boxes multiple time,
        // the particles positions should stay bellow this point.
        if any(position, |x| x.abs() > 1e6) {
            warn!(
                "Some particles have moved very far from the origin, \
                 the simulation might be exploding"
            );
            // we don't want to spam the output, so we return early if a
            // problem was found
            return;
        }
    }

    for velocity in system.particles().velocity {
        // Velocity threshold is 1000 A / fs
        if any(velocity, |x| x.abs() > 1000.0) {
            warn!(
                "Some particles have a very high velocity, \
                 the simulation might be exploding"
            );
            return;
        }
    }
}

fn any<F: Fn(f64) -> bool>(vector: &Vector3D, function: F) -> bool {
    function(vector[0]) || function(vector[1]) || function(vector[2])
}
