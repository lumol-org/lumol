// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

//! Type and algorithms for simulations
//!
//! The main stuct is [`Simulation`], containing all the data to run a single
//! simulation. A given `System` can be used with multiple [`Simulation`], for
//! example starting with an energy minimization, and then running some
//! Molecular dynamics;  
//!
//! [`Simulation`]: struct.Simulation.html

pub use lumol_sim::*;
