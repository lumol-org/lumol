// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

//! Algorithm to compute physical properties of a System

use std::f64::consts::PI;

use rayon::prelude::*;
use soa_derive::soa_zip;
use log_once::warn_once;

use crate::consts::K_BOLTZMANN;
use crate::{Matrix3, Vector3D};
use crate::{System, DegreesOfFreedom};

use crate::utils::ThreadLocalVec;

/// The `Compute` trait allow to compute properties of a system, without
/// modifying this system. The `Output` type is the type of the computed
/// property.
pub trait Compute {
    /// The data type of the property
    type Output;
    /// Compute the property
    fn compute(&self, system: &System) -> Self::Output;
}

/// Compute all the forces acting on the system, and return a vector of
/// force acting on each particles
pub struct Forces;
impl Compute for Forces {
    type Output = Vec<Vector3D>;
    fn compute(&self, system: &System) -> Vec<Vector3D> {
        let natoms = system.size();
        let thread_local_forces = ThreadLocalVec::with_size(natoms);

        (0..natoms).into_par_iter().for_each(|i| {
            let mut forces = thread_local_forces.borrow_mut();
            let mut force_i = Vector3D::zero();
            for j in (i + 1)..system.size() {
                let path = system.bond_path(i, j);
                let d = system.nearest_image(i, j);
                let dn = d.normalized();
                let r = d.norm();
                if let Some(potential) = system.pair_potential(i, j) {
                    let info = potential.restriction().information(path);
                    if !info.excluded {
                        let force = info.scaling * potential.force(r) * dn;
                        force_i += force;
                        forces[j] -= force;
                    }
                }
            }
            forces[i] += force_i;
        });

        // At this point all the forces are computed, but the results are
        // scattered across all thread local Vecs, here we gather them.
        let mut forces = vec![Vector3D::zero(); natoms];
        thread_local_forces.sum_into(&mut forces);

        for molecule in system.molecules() {
            for bond in molecule.bonds() {
                let (i, j) = (bond.i(), bond.j());
                let d = system.nearest_image(i, j);
                let dn = d.normalized();
                let r = d.norm();
                if let Some(potential) = system.bond_potential(i, j) {
                    let force = potential.force(r) * dn;
                    forces[i] += force;
                    forces[j] -= force;
                }
            }

            for angle in molecule.angles() {
                let (i, j, k) = (angle.i(), angle.j(), angle.k());
                let (theta, d1, d2, d3) = system.angle_and_derivatives(i, j, k);
                if let Some(potential) = system.angle_potential(i, j, k) {
                    let force = potential.force(theta);
                    forces[i] += force * d1;
                    forces[j] += force * d2;
                    forces[k] += force * d3;
                }
            }

            for dihedral in molecule.dihedrals() {
                let (i, j, k, m) = (dihedral.i(), dihedral.j(), dihedral.k(), dihedral.m());
                let (phi, d1, d2, d3, d4) = system.dihedral_and_derivatives(i, j, k, m);
                if let Some(potential) = system.dihedral_potential(i, j, k, m) {
                    let force = potential.force(phi);
                    forces[i] += force * d1;
                    forces[j] += force * d2;
                    forces[k] += force * d3;
                    forces[m] += force * d4;
                }
            }
        }

        if let Some(coulomb) = system.coulomb_potential() {
            coulomb.forces(system, &mut forces);
        }

        for global in system.global_potentials() {
            global.forces(system, &mut forces);
        }
        return forces;
    }
}

/// Compute the potential energy of the system
pub struct PotentialEnergy;
impl Compute for PotentialEnergy {
    type Output = f64;
    fn compute(&self, system: &System) -> f64 {
        let evaluator = system.energy_evaluator();

        let mut energy = evaluator.pairs();
        energy += evaluator.pairs_tail();
        energy += evaluator.bonds();
        energy += evaluator.angles();
        energy += evaluator.dihedrals();
        energy += evaluator.coulomb();
        energy += evaluator.global();

        assert!(energy.is_finite(), "Potential energy is infinite!");
        return energy;
    }
}

/// Compute the kinetic energy of the system
///
/// $$ K = \sum_i m_i \vec v_i \cdot \vec v_i $$
pub struct KineticEnergy;
impl Compute for KineticEnergy {
    type Output = f64;
    fn compute(&self, system: &System) -> f64 {
        let mut energy = 0.0;
        for (&mass, velocity) in soa_zip!(system.particles(), [mass, velocity]) {
            energy += 0.5 * mass * velocity.norm2();
        }
        assert!(energy.is_finite(), "Kinetic energy is infinite!");
        return energy;
    }
}

/// Compute the total energy of the system
pub struct TotalEnergy;
impl Compute for TotalEnergy {
    type Output = f64;
    fn compute(&self, system: &System) -> f64 {
        let kinetic = KineticEnergy.compute(system);
        let potential = PotentialEnergy.compute(system);
        return kinetic + potential;
    }
}

/// Compute the instantaneous temperature of the system
///
/// $$ T = \frac {2}{k_B N_f} \sum_i m_i \vec v_i \cdot \vec v_i $$
///
/// where $N_f$ is the number of degrees of freedom in the system, $k_B$ is the
/// Boltzman constant, $m_i$ the mass of particle $i$ and $\vec v_i$ the
/// velocity of particle $i$.
pub struct Temperature;
impl Compute for Temperature {
    type Output = f64;
    fn compute(&self, system: &System) -> f64 {
        let kinetic = KineticEnergy.compute(system);
        let dof = system.degrees_of_freedom() as f64;
        return 2.0 * kinetic / (dof * K_BOLTZMANN);
    }
}

/// Compute the volume of the system
pub struct Volume;
impl Compute for Volume {
    type Output = f64;
    #[inline]
    fn compute(&self, system: &System) -> f64 {
        return system.cell.volume();
    }
}


/// Compute the virial tensor of the system using the atomic definition.
///
/// $$ \underline{W} = \sum_i \vec r_i \otimes \vec f_i - \underline H
///    \frac{\partial U}{\partial \underline H} $$
///
/// where $\underline{H}$ is the unit cell matrix, $\vec f_i$ the force acting
/// on the atom $i$ and $\vec r_i$ the position of the atom $i$
///
/// If all the interactions are pair interactions, this definition reduces to
/// $$ \underline{W} = \sum_i \sum_{j > i} \vec r_{ij} \otimes \vec f_{ij} $$
pub struct AtomicVirial;
impl Compute for AtomicVirial {
    type Output = Matrix3;
    fn compute(&self, system: &System) -> Matrix3 {
        assert!(!system.cell.is_infinite(), "Can not compute virial for infinite cell");

        // Pair potentials contributions
        let pair_virials = (0..system.size()).into_par_iter().map(|i| {
            let mut local_virial = Matrix3::zero();
            for j in (i + 1)..system.size() {
                let path = system.bond_path(i, j);
                if let Some(potential) = system.pair_potential(i, j) {
                    let info = potential.restriction().information(path);
                    if !info.excluded {
                        let d = system.nearest_image(i, j);
                        local_virial += info.scaling * potential.virial(&d);
                    }
                }
            }
            return local_virial;
        });
        let mut virial = pair_virials.sum();

        // Tail correction for pair potentials contribution
        let volume = system.cell.volume();
        let composition = system.composition();
        for (i, ni) in composition.all_particles() {
            for (j, nj) in composition.all_particles() {
                let two_pi_density = 2.0 * PI * (ni as f64) * (nj as f64) / volume;
                if let Some(potential) = system.interactions().pair((i, j)) {
                    virial += two_pi_density * potential.tail_virial();
                }
            }
        }

        // Bond potentials contributions
        for molecule in system.molecules() {
            for bond in molecule.bonds() {
                let (i, j) = (bond.i(), bond.j());
                let r = system.nearest_image(i, j);
                if let Some(potential) = system.bond_potential(i, j) {
                    virial += potential.virial(&r);
                }
            }
        }

        // Angles and dihedrals potentials do not contribute as they only have
        // an angular part (see DL_POLY 4 manual page 18, or Smith, W., 1993,
        // CCP5 Information Quarterly, 39, 14. 18, 21, 24).

        if let Some(coulomb) = system.coulomb_potential() {
            virial += coulomb.atomic_virial(system);
        }

        for global in system.global_potentials() {
            virial += global.atomic_virial(system);
        }

        return virial;
    }
}

/// Compute the virial tensor of the system using the molecular definition
///
/// This differs from the [`AtomicVirial`](struct.AtomicVirial.html) when using
/// rigid molecules, as it will contains the right contributions of the forces
/// maintaining the molecules rigid without needing to compute them.
///
/// $$ \underline{W} = \sum_i \vec r_i \otimes \vec f_i - \underline H
///    \frac{\partial U}{\partial \underline H} $$
///
/// where $\underline{H}$ is the unit cell matrix, $\vec f_i$ the force acting
/// on the atom $i$ (comprising the force needed to keep the molecules rigid)
/// and $\vec r_i$ the position of the atom $i$
///
/// If all the interactions are pair interaction, this definition reduces to
///
/// $$ \underline{W} = \sum_i \sum_{j > i} \sum_{a \in i} \sum_{b \in i}
///    \frac{\vec r_{ab} \otimes \vec f_{ab}}{r_{ab}^2} \vec r_{ab} \cdot \vec
///    r_{ij} $$
///
/// where $i$ and $j$ run over all the molecules in the system, while $a$ and
/// $b$ run over all the particles in these molecules
pub struct MolecularVirial;
impl Compute for MolecularVirial {
    type Output = Matrix3;
    fn compute(&self, system: &System) -> Matrix3 {
        assert!(!system.cell.is_infinite(), "Can not compute virial for infinite cell");

        // Pair potentials contributions, using the molecular virial definition
        // This is defined in Allen & Tildesley in equations 2.54; 2.61; 2.63.
        let pair_virials = system.molecules().enumerate().par_bridge().map(|(i, molecule_i)| {
            let mut local_virial = Matrix3::zero();
            let ri = molecule_i.center_of_mass();

            for molecule_j in system.molecules().skip(i + 1) {
                let rj = molecule_j.center_of_mass();
                let mut r_ij = ri - rj;
                system.cell.vector_image(&mut r_ij);

                for part_a in molecule_i.indexes() {
                    for part_b in molecule_j.indexes() {
                        let path = system.bond_path(part_a, part_b);
                        let r_ab = system.nearest_image(part_a, part_b);
                        if let Some(potential) = system.pair_potential(part_a, part_b) {
                            let info = potential.restriction().information(path);
                            if !info.excluded {
                                let w_ab = info.scaling * potential.virial(&r_ab);
                                local_virial += w_ab * (r_ab * r_ij) / r_ab.norm2();
                            }
                        }
                     }
                 }
            }
            return local_virial;
        });
        let mut virial = pair_virials.sum();

        // Tail correction for pair potentials contribution
        let volume = system.cell.volume();
        let composition = system.composition();
        for (i, ni) in composition.all_particles() {
            for (j, nj) in composition.all_particles() {
                let two_pi_density = 2.0 * PI * (ni as f64) * (nj as f64) / volume;
                if let Some(potential) = system.interactions().pair((i, j)) {
                    virial += two_pi_density * potential.tail_virial();
                }
            }
        }

        // Bond potentials contributions
        for molecule in system.molecules() {
            for bond in molecule.bonds() {
                let (i, j) = (bond.i(), bond.j());
                let r = system.nearest_image(i, j);
                if let Some(potential) = system.bond_potential(i, j) {
                    let w = potential.virial(&r);
                    if w.norm() > 1e-30 {
                        // Use the same sorting as interactions
                        let name_i = &system.particles().name[i];
                        let name_j = &system.particles().name[j];
                        let (name_i, name_j) = if name_i < name_j {
                            (name_i, name_j)
                        } else {
                            (name_j, name_i)
                        };
                        warn_once!(
                            "Ignoring non null bond potential ({}, {}) in molecular virial",
                            name_i, name_j
                        );
                    }
                }
            }
        }

        // Angles and dihedrals potentials do not contribute as they only have
        // an angular part (see DL_POLY 4 manual page 18, or Smith, W., 1993,
        // CCP5 Information Quarterly, 39, 14. 18, 21, 24).

        if let Some(coulomb) = system.coulomb_potential() {
            virial += coulomb.molecular_virial(system);
        }

        for global in system.global_potentials() {
            virial += global.molecular_virial(system);
        }

        return virial;
    }
}

/// Compute the virial tensor of the system, picking between [`AtomicVirial`]
/// and [`MolecularVirial`] depending on the number of degrees of freedom
/// simulated on the system.
///
/// [`AtomicVirial`]: struct.AtomicVirial.html
/// [`MolecularVirial`]: struct.MolecularVirial.html
pub struct Virial;
impl Compute for Virial {
    type Output = Matrix3;
    fn compute(&self, system: &System) -> Matrix3 {
        match system.simulated_degrees_of_freedom {
            DegreesOfFreedom::Molecules => MolecularVirial.compute(system),
            DegreesOfFreedom::Particles | DegreesOfFreedom::Frozen(_) => AtomicVirial.compute(system),
        }
    }
}

/// Compute the pressure of the system using the virial definition, at a given
/// temperature.
///
/// $$ p = \frac{N_f k_B T}{3 V} + \frac{Tr(\underline{W})}{3V} $$
///
/// where $N_f$ is the number of degrees of freedom in the system, $k_B$ is the
/// Boltzman constant, $T$ the temperature, $V$ the simulation volume, $Tr$ is
/// the matricial trace, and $\underline{W}$ the [`Virial`].
///
/// [`Virial`]: struct.Virial.html
pub struct PressureAtTemperature {
    /// Temperature for the pressure computation
    pub temperature: f64,
}

impl Compute for PressureAtTemperature {
    type Output = f64;
    fn compute(&self, system: &System) -> f64 {
        assert!(!system.cell.is_infinite(), "Can not compute pressure for infinite cell");
        assert!(self.temperature >= 0.0);
        let virial = system.virial().trace();
        let volume = system.volume();
        let dof = system.degrees_of_freedom() as f64;
        return (dof * K_BOLTZMANN * self.temperature + virial) / (3.0 * volume);
    }
}

/// Compute the pressure of the system using the virial definition.
///
/// $$ p = \sum_i m_i \vec v_i \cdot \vec v_i + \frac{Tr(\underline{W})}{3V} $$
///
/// where $m_i$ is the mass of particle $i$, $\vec v_i$ the velocity of particle
/// $i$, $V$ the simulation volume, $Tr$ is the matricial trace, and
/// $\underline{W}$ the [`Virial`].
///
/// [`Virial`]: struct.Virial.html
pub struct Pressure;
impl Compute for Pressure {
    type Output = f64;
    fn compute(&self, system: &System) -> f64 {
        let pressure = PressureAtTemperature {
            temperature: system.temperature(),
        };
        return pressure.compute(system);
    }
}

/// Compute the stress tensor of the system from the virial definition, at the
/// given temperature.
///
/// $$ \underline{\sigma} = \frac{1}{V} \left(\frac{N_f}{3} k_B T \space
///    \underline{1} + \underline{W} \right) $$
///
/// where $N_f$ is the number of degrees of freedom in the system, $k_B$ is the
/// Boltzman constant, $T$ the temperature, $V$ the simulation volume, $Tr$ is
/// the matricial trace, and $\underline{W}$ the [`Virial`].
pub struct StressAtTemperature {
    /// Temperature for the stress tensor computation
    pub temperature: f64,
}

impl Compute for StressAtTemperature {
    type Output = Matrix3;
    fn compute(&self, system: &System) -> Matrix3 {
        assert!(self.temperature >= 0.0);
        assert!(!system.cell.is_infinite(), "Can not compute stress for infinite cell");
        let virial = system.virial();
        let volume = system.volume();
        let dof = system.degrees_of_freedom() as f64;
        let kinetic = dof / 3.0 * K_BOLTZMANN * self.temperature * Matrix3::one();
        return (kinetic + virial) / volume;
    }
}

/// Compute the stress tensor of the system from the virial definition
///
/// $$ \underline{\sigma} = \frac{1}{V} \left( \sum_i m_i \vec v_i \otimes \vec v_i + \underline{W} \right) $$
///
/// where $N_f$ is the number of degrees of freedom in the system, $k_B$ is the
/// Boltzman constant, $T$ the instantaneous system temperature, $V$ the
/// simulation volume, $Tr$ is the matricial trace, and $\underline{W}$ the
/// [`Virial`].
pub struct Stress;
impl Compute for Stress {
    type Output = Matrix3;
    fn compute(&self, system: &System) -> Matrix3 {
        assert!(!system.cell.is_infinite(), "Can not compute stress for infinite cell");

        let mut kinetic = Matrix3::zero();
        for (&mass, velocity) in soa_zip!(system.particles(), [mass, velocity]) {
            kinetic += mass * velocity.tensorial(velocity);
        }

        let volume = system.volume();
        let virial = system.virial();
        return (kinetic + virial) / volume;
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::System;
    use crate::consts::K_BOLTZMANN;
    use crate::{Harmonic, NullPotential, PairInteraction};
    use crate::utils::system_from_xyz;
    use crate::units;

    use approx::assert_ulps_eq;

    fn test_pairs_system() -> System {
        let mut system = system_from_xyz(
            "2
            cell: 10.0
            F 0.0 0.0 0.0 -0.007225222699367925 -0.002405756495275919  0.0026065109398392215
            F 1.3 0.0 0.0  0.001179633958023287  0.003525262341736351 -0.0004132774783154952
            ",
        );

        let mut interaction = PairInteraction::new(
            Box::new(Harmonic {
                k: units::from(300.0, "kJ/mol/A^2").unwrap(),
                x0: units::from(1.2, "A").unwrap(),
            }),
            5.0,
        );
        interaction.enable_tail_corrections();
        system.set_pair_potential(("F", "F"), interaction);

        // unused interaction to check that we do handle this right
        system.set_pair_potential(("H", "O"), PairInteraction::new(Box::new(NullPotential), 0.0));


        return system;
    }

    fn test_molecular_system() -> System {
        let mut system = system_from_xyz(
            "4
            cell: 10.0
            F 0.0 0.0 0.0
            F 1.0 0.0 0.0
            F 1.0 1.0 0.0
            F 2.0 1.0 0.0
            ",
        );
        assert!(system.add_bond(0, 1).is_empty());
        assert!(system.add_bond(1, 2).is_empty());
        assert!(system.add_bond(2, 3).is_empty());
        assert_eq!(system.molecules().count(), 1);
        assert_eq!(system.molecule(0).bonds().len(), 3);

        system.set_pair_potential(("F", "F"), PairInteraction::new(Box::new(NullPotential), 0.0));

        system.set_bond_potential(
            ("F", "F"),
            Box::new(Harmonic {
                k: units::from(100.0, "kJ/mol/A^2").unwrap(),
                x0: units::from(2.0, "A").unwrap(),
            }),
        );

        system.set_angle_potential(
            ("F", "F", "F"),
            Box::new(Harmonic {
                k: units::from(100.0, "kJ/mol/deg^2").unwrap(),
                x0: units::from(88.0, "deg").unwrap(),
            }),
        );

        system.set_dihedral_potential(
            ("F", "F", "F", "F"),
            Box::new(Harmonic {
                k: units::from(100.0, "kJ/mol/deg^2").unwrap(),
                x0: units::from(185.0, "deg").unwrap(),
            }),
        );

        // unused interaction to check that we do handle this right
        system.set_pair_potential(("H", "O"), PairInteraction::new(Box::new(NullPotential), 0.0));

        return system;
    }

    #[test]
    fn forces_pairs() {
        let system = &test_pairs_system();
        let res = Forces.compute(system);

        let forces_tot = res[0] + res[1];
        assert_eq!(forces_tot, Vector3D::zero());

        let force = units::from(30.0, "kJ/mol/A").unwrap();
        assert_ulps_eq!(res[0][0], force);
        assert_ulps_eq!(res[0][1], 0.0);
        assert_ulps_eq!(res[0][1], 0.0);

        assert_ulps_eq!(res[1][0], -force);
        assert_ulps_eq!(res[1][1], 0.0);
        assert_ulps_eq!(res[1][1], 0.0);
    }

    #[test]
    fn force_molecular() {
        let system = test_molecular_system();
        let res = Forces.compute(&system);
        let forces_tot = res[0] + res[1] + res[2] + res[3];
        assert_ulps_eq!(forces_tot.norm2(), 0.0);
    }

    #[test]
    #[allow(clippy::unreadable_literal)]
    fn energy_pairs() {
        let system = &test_pairs_system();
        let kinetic = KineticEnergy.compute(system);
        let potential = PotentialEnergy.compute(system);
        let total = TotalEnergy.compute(system);

        assert_eq!(kinetic + potential, total);
        assert_ulps_eq!(kinetic, 0.0007483016557453698);

        assert_eq!(kinetic, system.kinetic_energy());
        assert_eq!(potential, system.potential_energy());
        assert_eq!(total, system.total_energy());
    }

    #[test]
    fn energy_molecular() {
        let system = test_molecular_system();
        assert_ulps_eq!(PotentialEnergy.compute(&system), units::from(1800.0, "kJ/mol").unwrap());
    }

    #[test]
    fn temperature() {
        let system = &test_pairs_system();
        let temperature = Temperature.compute(system);
        assert_ulps_eq!(temperature, 300.0);
        assert_eq!(temperature, system.temperature());
    }

    #[test]
    fn volume() {
        let system = &test_pairs_system();
        let volume = Volume.compute(system);
        assert_eq!(volume, 1000.0);
        assert_eq!(volume, system.volume());
    }

    #[test]
    #[should_panic(expected="Can not compute virial for infinite cell")]
    fn virial_infinite_cell() {
        let _ = Virial.compute(&System::new());
    }

    #[test]
    fn virial_pairs() {
        let system = &test_pairs_system();
        let virial = Virial.compute(system);

        let mut expected = Matrix3::zero();
        let force = units::from(30.0, "kJ/mol/A^2").unwrap();
        expected[0][0] = -force * 1.3;

        assert_ulps_eq!(virial, expected);
        assert_eq!(virial, system.virial());
    }

    #[test]
    fn virial_molecular() {
        let system = &test_molecular_system();
        let virial = Virial.compute(system);

        let mut expected = Matrix3::zero();
        let w = units::from(100.0, "kJ/mol/A").unwrap();
        expected[0][0] = 2.0 * w;
        expected[1][1] = 1.0 * w;

        assert_ulps_eq!(virial, expected);
        assert_eq!(virial, system.virial());
    }

    #[test]
    #[should_panic(expected="assertion failed: self.temperature >= 0.0")]
    fn pressure_at_temperature_negative_temperature() {
        let system = &test_pairs_system();
        let pressure = PressureAtTemperature { temperature: -4.0 };
        let _ = pressure.compute(system);
    }

    #[test]
    #[should_panic(expected="Can not compute pressure for infinite cell")]
    fn pressure_at_temperature_infinite_cell() {
        let pressure = PressureAtTemperature { temperature: -4.0 };
        let _ = pressure.compute(&System::new());
    }

    #[test]
    fn pressure_at_temperature() {
        let system = &mut test_pairs_system();

        let temperature = 550.0;
        let force = units::from(30.0, "kJ/mol/A").unwrap();
        let virial = -force * 1.3;
        let natoms = 2.0;
        let volume = 1000.0;

        // Direct computation
        let expected = natoms * K_BOLTZMANN * temperature / volume + virial / (3.0 * volume);
        let pressure = PressureAtTemperature {
            temperature: temperature,
        };
        let pressure = pressure.compute(system);
        assert_ulps_eq!(pressure, expected);

        // Computation with the real system temperature
        let pressure = PressureAtTemperature {
            temperature: system.temperature(),
        };
        let pressure = pressure.compute(system);
        assert_eq!(pressure, system.pressure());

        // Computation with an external temperature for the system
        let pressure = PressureAtTemperature {
            temperature: temperature,
        };
        let pressure = pressure.compute(system);
        system.simulated_temperature(Some(temperature));
        assert_eq!(pressure, system.pressure());
    }

    #[test]
    #[should_panic(expected="assertion failed: self.temperature >= 0.0")]
    fn stress_at_temperature_negative_temperature() {
        let system = &test_pairs_system();
        let stress = StressAtTemperature { temperature: -4.0 };
        let _ = stress.compute(system);
    }

    #[test]
    #[should_panic(expected="Can not compute stress for infinite cell")]
    fn stress_at_temperature_infinite_cell() {
        let stress = StressAtTemperature { temperature: 300.0 };
        let _ = stress.compute(&System::new());
    }

    #[test]
    fn stress_at_temperature() {
        let system = &mut test_pairs_system();
        let temperature = 550.0;
        let stress = StressAtTemperature {
            temperature: temperature,
        };
        let stress = stress.compute(system);
        let pressure = PressureAtTemperature {
            temperature: temperature,
        };
        let pressure = pressure.compute(system);

        // tail corrections are smaller than 1e-9
        let trace = (stress[0][0] + stress[1][1] + stress[2][2]) / 3.0;
        assert_ulps_eq!(trace, pressure);

        system.simulated_temperature(Some(temperature));
        assert_eq!(stress, system.stress());
    }

    #[test]
    #[should_panic(expected="Can not compute stress for infinite cell")]
    fn stress_infinite_cell() {
        let _ = Stress.compute(&System::new());
    }

    #[test]
    fn stress() {
        let system = &test_pairs_system();
        let stress = Stress.compute(system);
        let pressure = Pressure.compute(system);

        let trace = stress.trace() / 3.0;
        assert_ulps_eq!(trace, pressure);
        assert_eq!(stress, system.stress());
    }

    #[test]
    #[should_panic(expected="Can not compute pressure for infinite cell")]
    fn pressure_infinite_cell() {
        let _ = Pressure.compute(&System::new());
    }

    #[test]
    fn pressure() {
        let system = &test_pairs_system();
        let pressure = Pressure.compute(system);

        let force = units::from(30.0, "kJ/mol/A").unwrap();
        let virial = -force * 1.3;
        let natoms = 2.0;
        let temperature = 300.0;
        let volume = 1000.0;
        let expected = natoms * K_BOLTZMANN * temperature / volume + virial / (3.0 * volume);

        assert_ulps_eq!(pressure, expected);
        assert_eq!(pressure, system.pressure());
    }
}
