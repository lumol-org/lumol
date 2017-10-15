// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license
#![allow(dead_code)]

use lumol::out::Output;
use lumol::sys::System;

use std::rc::Rc;
use std::sync::RwLock;

pub type SharedVec = Rc<RwLock<Vec<f64>>>;

/// Collect pressure and temperature of a simulation after a starting step
pub struct Collecter {
    start: u64,
    pressures: SharedVec,
    temperatures: SharedVec,
}

impl Collecter {
    pub fn new(start: u64) -> Collecter {
        Collecter {
            start: start,
            pressures: Rc::new(RwLock::new(Vec::new())),
            temperatures: Rc::new(RwLock::new(Vec::new())),
        }
    }

    pub fn temperatures(&self) -> SharedVec {
        self.temperatures.clone()
    }

    pub fn pressures(&self) -> SharedVec {
        self.pressures.clone()
    }
}

impl Output for Collecter {
    fn write(&mut self, system: &System) {
        if system.step() < self.start {
            return
        }

        self.pressures.write().unwrap().push(system.pressure());
        self.temperatures.write().unwrap().push(system.temperature());
    }
}

pub fn mean(data: SharedVec) -> f64 {
    let data = data.read().unwrap();
    data.iter().sum::<f64>() / data.len() as f64
}
