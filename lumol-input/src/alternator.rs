// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

use lumol_sim::md::Control;
use lumol_core::System;

/// Helper struct that can wrap an algorithm to make it run only a fraction of
/// the times it is called.
#[derive(Debug)]
pub struct Alternator<T> {
    every: u64,
    count: u64,
    inner: T,
}

impl<T> Alternator<T> {
    /// Wrap the algorithm `base` to call it only every `every` time.
    pub fn new(every: u64, inner: T) -> Alternator<T> {
        Alternator {
            every: every,
            inner: inner,
            count: 0,
        }
    }

    /// Check if this is the appropriate time to run the algorithm.
    pub fn can_run(&mut self) -> bool {
        self.count += 1;
        self.count % self.every == 0
    }
}

impl<T> AsRef<T> for Alternator<T> {
    fn as_ref(&self) -> &T {
        &self.inner
    }
}

impl<T> AsMut<T> for Alternator<T> {
    fn as_mut(&mut self) -> &mut T {
        &mut self.inner
    }
}

impl<T: Control> Control for Alternator<T> {
    fn control(&mut self, system: &mut System) {
        if self.can_run() {
            self.as_mut().control(system);
        }
    }
}
