// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license
#![allow(dead_code, clippy::needless_return)]

use lumol::{System, TrajectoryBuilder, Vector3D};
use lumol::input::InteractionsInput;
use std::path::Path;

use rand::{SeedableRng, Rng};
use rand_xorshift::XorShiftRng;

pub fn get_system(name: &str) -> System {
    let data = Path::new(file!()).parent().unwrap().join("..").join("data");

    let system = data.join(String::from(name) + ".pdb");
    let mut system = TrajectoryBuilder::new().open(system)
                                             .and_then(|mut trajectory| trajectory.read())
                                             .unwrap();

    let interactions = data.join(String::from(name) + ".toml");
    InteractionsInput::new(interactions).and_then(|input| input.read(&mut system))
                                        .unwrap();

    return system;
}

pub fn get_rng() -> XorShiftRng {
    XorShiftRng::from_seed([145, 59, 58, 50, 238, 182, 97, 28, 107, 149, 227, 40, 90, 109, 196, 129])
}

pub fn move_rigid_molecule(system: &System) -> (usize, Vec<Vector3D>) {
    let mut rng = get_rng();

    let molid = rng.gen_range(0..system.molecules().count());
    let molecule = system.molecule(molid);
    let delta = Vector3D::new(rng.gen(), rng.gen(), rng.gen());
    let mut positions = Vec::new();
    for position in molecule.particles().position {
        positions.push(position + delta);
    }

    return (molid, positions);
}

pub fn move_all_rigid_molecule(mut system: System) -> System {
    let mut rng = get_rng();

    for mut molecule in system.molecules_mut() {
        let delta = Vector3D::new(rng.gen(), rng.gen(), rng.gen());
        for position in molecule.particles_mut().position {
            *position += delta;
        }
    }

    return system;
}
