use sys::{System, Molecule};
use types::Vector3D;

/// Computes energy between a particle and the surrounding system.
fn trial_external_energy(system: &System, trialpos: &Vector3D, pid: usize) -> f64 {
    let mut energy = 0.0;
    let cell = system.cell();
    // maybe rewrite this loop to go over molecules
    for i in 0..system.size() {
        if system.are_in_same_molecule(i, pid) { break };
        let r = cell.distance(trialpos, &system[i].position);
        for potential in system.pair_potentials(pid, i) {
            let info = potential.restriction().information(-1);
            if !info.excluded {
                energy += info.scaling * potential.energy(r);
            }
        }
    }
    // TODO: add global interactions (short ranged only?)
    // Figure out how to update cache accordingly
    energy
}

/// Returns non-covalent energy contribution of the (partially) grown molecule.
pub fn trial_non_covalent_energy(system: &System, partialpos: &[Vector3D], trialpos: &Vector3D, molecule: &Molecule) -> f64 {
    let mut energy = 0.0;

    let length = partialpos.len();
    // only one particle? nothing to do.
    if length == 1 { return energy };
    // index of trial: length + start:
    // o -- o -- O <- trial
    // 0    1    2  == length
    let j = length + molecule.start(); // to increase readability in loop
    // iterate over partially grown molecule
    for (pos_i, i) in partialpos.iter().zip(molecule.iter()) {
        let r = system.cell().distance(pos_i, trialpos);
        let distance = system.bond_distance(i, j); // distance within molecule
        for potential in system.pair_potentials(i, j) {
            let info = potential.restriction().information(distance);
            if !info.excluded {
                energy += info.scaling * potential.energy(r);
            }
        }
    }

    // TODO ELECTROSTATICSSSAAAAAAGGGGGG
    energy += trial_external_energy(&system, &trialpos, j);
    energy
}


