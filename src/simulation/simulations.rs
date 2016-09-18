// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

use system::System;

use super::Propagator;
use super::Output;

/// Writting an output at a given frequency
struct OutputFrequency {
    /// The output to use
    output: Box<Output>,
    /// The frequency. `output` will be used everytime the system step matches
    /// this frequency.
    frequency: u64,
}

impl OutputFrequency {
    pub fn new(output: Box<Output>) -> OutputFrequency {
        OutputFrequency{
            frequency: 1,
            output: output,
        }
    }

    pub fn with_frequency(output: Box<Output>, frequency: u64) -> OutputFrequency {
        OutputFrequency{
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
        if system.step() % self.frequency == 0 {
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
    propagator: Box<Propagator>,
    outputs: Vec<OutputFrequency>
}

impl Simulation {
    /// Create a new Simulation from a Propagator.
    pub fn new(propagator: Box<Propagator>) -> Simulation {
        Simulation {
            propagator: propagator,
            outputs: Vec::new(),
        }
    }

    /// Run the simulation on System for `nsteps` steps.
    pub fn run(&mut self, system: &mut System, nsteps: usize) {
        self.setup(system);
        for _ in 0..nsteps {
            self.propagator.propagate(system);
            system.increment_step();
            for output in &mut self.outputs {
                output.write(system);
            }
        }
        self.finish(system);
    }

    /// Add a new `Output` algorithm in the outputs list
    pub fn add_output(&mut self, output: Box<Output>) {
        self.outputs.push(OutputFrequency::new(output));
    }

    /// Add a new `Output` algorithm in the outputs list, which will be used
    /// at the given frequency. The output will be used everytime the system
    ///  step matches this frequency.
    pub fn add_output_with_frequency(&mut self, output: Box<Output>, frequency: u64) {
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

    #[cfg(test)]
    pub fn outputs_len(&self) -> usize {
        self.outputs.len()
    }
}
