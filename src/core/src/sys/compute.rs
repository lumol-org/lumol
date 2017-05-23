// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

//! Computing properties of a system

use std::f64::consts::PI;

use consts::K_BOLTZMANN;
use types::{Matrix3, Vector3D, Zero, One};
use sys::System;

use parallel::prelude::*;
use parallel::ThreadLocalStore;

/// The compute trait allow to compute properties of a system, without
/// modifying this system. The Output type is the type of the computed
/// property.
pub trait Compute {
    /// The data type of the property
    type Output;
    /// Compute the property
    fn compute(&self, system: &System) -> Self::Output;
}

/******************************************************************************/
/// Compute all the forces acting on the system, and return a vector of
/// force acting on each particles
pub struct Forces;
impl Compute for Forces {
    type Output = Vec<Vector3D>;
    fn compute(&self, system: &System) -> Vec<Vector3D> {
        let natoms = system.size();
        let thread_forces_store = ThreadLocalStore::new(|| vec![Vector3D::zero(); natoms]);

        (0..natoms).into_par_iter().for_each(|i| {

            let mut thread_forces = thread_forces_store.borrow_mut();

            for j in (i+1)..system.size() {
                let distance = system.bond_distance(i, j);
                let d = system.nearest_image(i, j);
                let dn = d.normalized();
                let r = d.norm();
                for potential in system.pair_potentials(i, j) {
                    let info = potential.restriction().information(distance);
                    if !info.excluded {
                        let force = info.scaling * potential.force(r) * dn;
                        thread_forces[i] += force;
                        thread_forces[j] -= force;
                    }
                }
            }
        });

        // At this point all the forces are computed, but the
        // results are scattered across all thread local Vecs,
        // here we gather them.
        let mut forces = vec![Vector3D::zero(); natoms];
        thread_forces_store.sum_local_values(&mut forces);

        for molecule in system.molecules() {
            for bond in molecule.bonds() {
                let (i, j) = (bond.i(), bond.j());
                let d = system.nearest_image(i, j);
                let dn = d.normalized();
                let r = d.norm();
                for potential in system.bond_potentials(i, j) {
                    let force = potential.force(r) * dn;
                    forces[i] += force;
                    forces[j] -= force;
                }
            }

            for angle in molecule.angles() {
                let (i, j, k) = (angle.i(), angle.j(), angle.k());
                let (theta, d1, d2, d3) = system.angle_and_derivatives(i, j, k);
                for potential in system.angle_potentials(i, j, k) {
                    let force = potential.force(theta);
                    forces[i] += force * d1;
                    forces[j] += force * d2;
                    forces[k] += force * d3;
                }
            }

            for dihedral in molecule.dihedrals() {
                let (i, j, k, m) = (dihedral.i(), dihedral.j(), dihedral.k(), dihedral.m());
                let (phi, d1, d2, d3, d4) = system.dihedral_and_derivatives(i, j, k, m);
                for potential in system.dihedral_potentials(i, j, k, m) {
                    let force = potential.force(phi);
                    forces[i] += force * d1;
                    forces[j] += force * d2;
                    forces[k] += force * d3;
                    forces[m] += force * d4;
                }
            }
        }

        if let Some(coulomb) = system.interactions().coulomb() {
            let coulomb_forces = coulomb.forces(system);
            debug_assert_eq!(coulomb_forces.len(), natoms, "Wrong `forces` size in coulomb potentials");
            for (i, force) in coulomb_forces.iter().enumerate() {
                forces[i] += force;
            }
        }

        for global in system.interactions().globals() {
            let global_forces = global.forces(system);
            debug_assert_eq!(global_forces.len(), natoms, "Wrong `forces` size in global potentials");
            for (i, force) in global_forces.iter().enumerate() {
                forces[i] += force;
            }
        }
        return forces;
    }
}

/******************************************************************************/
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

/******************************************************************************/
/// Compute the kinetic energy of the system
pub struct KineticEnergy;
impl Compute for KineticEnergy {
    type Output = f64;
    fn compute(&self, system: &System) -> f64 {
        let mut energy = 0.0;
        for particle in system {
            energy += 0.5 * particle.mass * particle.velocity.norm2();
        }
        assert!(energy.is_finite(), "Kinetic energy is infinite!");
        return energy;
    }
}

/******************************************************************************/
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

/******************************************************************************/
/// Compute the instantaneous temperature of the system
pub struct Temperature;
impl Compute for Temperature {
    type Output = f64;
    fn compute(&self, system: &System) -> f64 {
        let kinetic = KineticEnergy.compute(system);
        let natoms = system.size() as f64;
        return 1.0/K_BOLTZMANN * 2.0 * kinetic/(3.0 * natoms);
    }
}

/******************************************************************************/
/// Compute the volume of the system
pub struct Volume;
impl Compute for Volume {
    type Output = f64;
    #[inline]
    fn compute(&self, system: &System) -> f64 {
        return system.cell().volume();
    }
}

/******************************************************************************/
/// Compute the virial tensor of the system, defined by
/// $$ W = \sum_i \sum_{j > i} \vec r_{ij} \otimes \vec f_{ij} $$
pub struct Virial;
impl Compute for Virial {
    type Output = Matrix3;
    fn compute(&self, system: &System) -> Matrix3 {
        assert!(!system.cell().is_infinite(), "Can not compute virial for infinite cell");

        let mut virial = (0..system.size()).par_map(|i| {
            let mut local_virial = Matrix3::zero();

            for j in (i+1)..system.size() {
                let distance = system.bond_distance(i, j);
                for potential in system.pair_potentials(i, j) {
                    let info = potential.restriction().information(distance);
                    if !info.excluded {
                        let d = system.nearest_image(i, j);
                        local_virial += info.scaling * potential.virial(&d);
                    }
                }
            }

            local_virial
        }).sum();

        let volume = system.cell().volume();
        let composition = system.composition();
        for i in system.particle_kinds() {
            let ni = composition[i] as f64;
            for j in system.particle_kinds() {
                let nj = composition[j] as f64;
                let potentials = system.interactions().pairs(i, j);
                for potential in potentials {
                    virial += 2.0 * PI * ni * nj * potential.tail_virial() / volume;
                }
            }
        }

        // TODO: implement virial computations for molecular potentials
        // (angles & dihedrals)

        if let Some(coulomb) = system.interactions().coulomb() {
            virial += coulomb.virial(system);
        }

        for global in system.interactions().globals() {
            virial += global.virial(system);
        }

        return virial;
    }
}

/******************************************************************************/
/// Compute the pressure of the system from the virial equation, at the given
/// temperature. This pressure is given by the following formula:
/// $$ p = \frac{N k_B T}{V} + \frac{1}{3V} \sum_i \vec f_i \cdot \vec r_i $$
pub struct PressureAtTemperature {
    /// Temperature for the pressure computation
    pub temperature: f64
}

impl Compute for PressureAtTemperature {
    type Output = f64;
    fn compute(&self, system: &System) -> f64 {
        assert!(!system.cell().is_infinite(), "Can not compute pressure for infinite cell");
        assert!(self.temperature >= 0.0);
        let virial_tensor = system.virial();
        let virial = virial_tensor.trace();
        let volume = system.cell().volume();
        let natoms = system.size() as f64;
        return natoms * K_BOLTZMANN * self.temperature / volume + virial / (3.0 * volume);
    }
}

/******************************************************************************/
/// Compute the stress tensor of the system from the virial equation, at the
/// given temperature. The stress tensor is defined by
/// $$ \sigma = \sigma = \frac{1}{V} (\sum_i m_i v_i \otimes v_i + \sum_i \sum_{j > i} \vec r_{ij} \otimes \vec f_{ij}) $$
/// but here the kinetic energy term is replaced by it average at the given
/// temperature.
pub struct StressAtTemperature {
    /// Temperature for the stress tensor computation
    pub temperature: f64
}

impl Compute for StressAtTemperature {
    type Output = Matrix3;
    fn compute(&self, system: &System) -> Matrix3 {
        assert!(self.temperature >= 0.0);
        assert!(!system.cell().is_infinite(), "Can not compute stress for infinite cell");
        let virial = system.virial();
        let volume = system.cell().volume();
        let natoms = system.size() as f64;
        let kinetic = natoms * K_BOLTZMANN * self.temperature * Matrix3::one();
        return (kinetic + virial) / volume;
    }
}

/******************************************************************************/
/// Compute the stress tensor of the system, defined by:
/// $$ \sigma = \frac{1}{V} (\sum_i m_i v_i \otimes v_i + \sum_i \sum_{j > i} \vec r_{ij} \otimes \vec f_{ij}) $$
pub struct Stress;
impl Compute for Stress {
    type Output = Matrix3;
    fn compute(&self, system: &System) -> Matrix3 {
        assert!(!system.cell().is_infinite(), "Can not compute stress for infinite cell");
        let mut kinetic = Matrix3::zero();
        for particle in system.iter() {
            let velocity = &particle.velocity;
            kinetic += particle.mass * velocity.tensorial(velocity);
        }

        let volume = system.cell().volume();
        let virial = system.virial();
        return (kinetic + virial) / volume;
    }
}

/******************************************************************************/
/// Compute the virial pressure of the system. This pressure is given by the
/// following formula:
/// $$ p = \frac{N k_B T}{V} + \frac{1}{3V} \sum_i \vec f_i \cdot \vec r_i $$
pub struct Pressure;
impl Compute for Pressure {
    type Output = f64;
    fn compute(&self, system: &System) -> f64 {
        return PressureAtTemperature{temperature: system.temperature()}.compute(system);
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use types::*;
    use sys::{System, Particle, UnitCell};
    use sys::veloc::{InitVelocities, BoltzmannVelocities};
    use energy::{Harmonic, NullPotential, PairInteraction};
    use consts::K_BOLTZMANN;
    use utils::unit_from;

    fn test_pairs_system() -> System {
        let mut system = System::from_cell(UnitCell::cubic(10.0));;
        system.add_particle(Particle::new("F"));
        system[0].position = Vector3D::zero();
        system.add_particle(Particle::new("F"));
        system[1].position = Vector3D::new(1.3, 0.0, 0.0);

        let mut interaction = PairInteraction::new(Box::new(Harmonic{
            k: unit_from(300.0, "kJ/mol/A^2"),
            x0: unit_from(1.2, "A")
        }), 5.0);
        interaction.enable_tail_corrections();
        system.interactions_mut().add_pair("F", "F", interaction);

        /// unused interaction to check that we do handle this right
        system.interactions_mut().add_pair("H", "O",
            PairInteraction::new(Box::new(NullPotential), 0.0)
        );

        let mut velocities = BoltzmannVelocities::new(unit_from(300.0, "K"));
        velocities.init(&mut system);
        return system;
    }

    fn test_molecular_system() -> System {
        let mut system = System::from_cell(UnitCell::cubic(10.0));;
        system.add_particle(Particle::new("F"));
        system[0].position = Vector3D::new(0.0, 0.0, 0.0);
        system.add_particle(Particle::new("F"));
        system[1].position = Vector3D::new(1.0, 0.0, 0.0);
        system.add_particle(Particle::new("F"));
        system[2].position = Vector3D::new(1.0, 1.0, 0.0);
        system.add_particle(Particle::new("F"));
        system[3].position = Vector3D::new(2.0, 1.0, 0.0);

        assert!(system.add_bond(0, 1).is_empty());
        assert!(system.add_bond(1, 2).is_empty());
        assert!(system.add_bond(2, 3).is_empty());

        system.interactions_mut().add_pair("F", "F",
            PairInteraction::new(Box::new(NullPotential), 0.0)
        );

        system.interactions_mut().add_bond("F", "F",
            Box::new(Harmonic{
                k: unit_from(100.0, "kJ/mol/A^2"),
                x0: unit_from(2.0, "A")
        }));

        system.interactions_mut().add_angle("F", "F", "F",
            Box::new(Harmonic{
                k: unit_from(100.0, "kJ/mol/deg^2"),
                x0: unit_from(88.0, "deg")
        }));

        system.interactions_mut().add_dihedral("F", "F", "F", "F",
            Box::new(Harmonic{
                k: unit_from(100.0, "kJ/mol/deg^2"),
                x0: unit_from(185.0, "deg")
        }));

        /// unused interaction to check that we do handle this right
        system.interactions_mut().add_pair("H", "O",
            PairInteraction::new(Box::new(NullPotential), 0.0)
        );

        return system;
    }

    #[test]
    fn forces_pairs() {
        let system = &test_pairs_system();
        let res = Forces.compute(system);

        let forces_tot = res[0] + res[1];
        assert_eq!(forces_tot, Vector3D::zero());

        let force = unit_from(30.0, "kJ/mol/A");
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
    fn energy_pairs() {
        let system = &test_pairs_system();
        let kinetic = KineticEnergy.compute(system);
        let potential = PotentialEnergy.compute(system);
        let total = TotalEnergy.compute(system);

        assert_eq!(kinetic + potential, total);
        assert_eq!(kinetic, 0.0007483016557453699);

        assert_eq!(kinetic, system.kinetic_energy());
        assert_eq!(potential, system.potential_energy());
        assert_eq!(total, system.total_energy());
    }

    #[test]
    fn energy_molecular() {
        let system = test_molecular_system();
        assert_ulps_eq!(PotentialEnergy.compute(&system), unit_from(1800.0, "kJ/mol"));
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
    #[should_panic]
    fn virial_infinite_cell() {
        let _ = Virial.compute(&System::new());
    }

    #[test]
    fn virial() {
        let system = &test_pairs_system();
        let virial = Virial.compute(system);

        let mut expected = Matrix3::zero();
        let force = unit_from(30.0, "kJ/mol/A");
        expected[(0, 0)] = - force * 1.3;

        assert_ulps_eq!(virial, expected);
        assert_eq!(virial, system.virial());
    }

    #[test]
    #[should_panic]
    fn pressure_at_temperature_negative_temperature() {
        let system = &test_pairs_system();
        let pressure = PressureAtTemperature{temperature: -4.0};
        let _ = pressure.compute(system);
    }

    #[test]
    #[should_panic]
    fn pressure_at_temperature_infinite_cell() {
        let pressure = PressureAtTemperature{temperature: -4.0};
        let _ = pressure.compute(&System::new());
    }

    #[test]
    fn pressure_at_temperature() {
        let system = &mut test_pairs_system();

        let temperature = 550.0;
        let force = unit_from(30.0, "kJ/mol/A");
        let virial = -force * 1.3;
        let natoms = 2.0;
        let volume = 1000.0;

        // Direct computation
        let expected = natoms * K_BOLTZMANN * temperature / volume + virial / (3.0 * volume);
        let pressure = PressureAtTemperature{temperature: temperature}.compute(system);
        assert_ulps_eq!(pressure, expected);

        // Computation with the real system temperature
        let pressure = PressureAtTemperature{temperature: system.temperature()};
        let pressure = pressure.compute(system);
        assert_eq!(pressure, system.pressure());

        // Computation with an external temperature for the system
        let pressure = PressureAtTemperature{temperature: temperature}.compute(system);
        system.external_temperature(Some(temperature));
        assert_eq!(pressure, system.pressure());
    }

    #[test]
    #[should_panic]
    fn stress_at_temperature_negative_temperature() {
        let system = &test_pairs_system();
        let stress = StressAtTemperature{temperature: -4.0};
        let _ = stress.compute(system);
    }

    #[test]
    #[should_panic]
    fn stress_at_temperature_infinite_cell() {
        let stress = StressAtTemperature{temperature: 300.0};
        let _ = stress.compute(&System::new());
    }

    #[test]
    fn stress_at_temperature() {
        let system = &mut test_pairs_system();
        let temperature = 550.0;
        let stress = StressAtTemperature{temperature: temperature}.compute(system);
        let pressure = PressureAtTemperature{temperature: temperature}.compute(system);

        // tail corrections are smaller than 1e-9
        let trace = (stress[(0, 0)] + stress[(1, 1)] + stress[(2, 2)]) / 3.0;
        assert_ulps_eq!(trace, pressure);

        system.external_temperature(Some(temperature));
        assert_eq!(stress, system.stress());
    }

    #[test]
    #[should_panic]
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
    #[should_panic]
    fn pressure_infinite_cell() {
        let _ = Pressure.compute(&System::new());
    }

    #[test]
    fn pressure() {
        let system = &test_pairs_system();
        let pressure = Pressure.compute(system);

        let force = unit_from(30.0, "kJ/mol/A");
        let virial = - force * 1.3;
        let natoms = 2.0;
        let temperature = 300.0;
        let volume = 1000.0;
        let expected = natoms * K_BOLTZMANN * temperature / volume + virial / (3.0 * volume);

        assert_ulps_eq!(pressure, expected);
        assert_eq!(pressure, system.pressure());
    }
}
