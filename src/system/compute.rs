// Cymbalum, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

//! Computing properties of a system
use constants::K_BOLTZMANN;
use types::{Matrix3, Vector3D, Zero, One};
use system::System;
use potentials::{Potential, Virial as VirialTrait};

/// The compute trait allow to compute properties of a system, whithout
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
        let mut res = vec![Vector3D::zero(); natoms];

        for i in 0..system.size() {
            for j in (i+1)..system.size() {
                let d = system.nearest_image(i, j);
                let dn = d.normalized();
                let r = d.norm();
                for potential in system.pair_potentials(i, j) {
                    let info = potential.restriction().informations(system, i, j);
                    if !info.excluded {
                        let force = info.scaling * potential.force(r) * dn;
                        res[i] += force;
                        res[j] -= force;
                    }
                }
            }
        }

        for molecule in system.molecules() {
            for bond in molecule.bonds() {
                let (i, j) = (bond.i(), bond.j());
                let d = system.nearest_image(i, j);
                let dn = d.normalized();
                let r = d.norm();
                for potential in system.bond_potentials(i, j) {
                    let force = potential.force(r) * dn;
                    res[i] += force;
                    res[j] -= force;
                }
            }

            for angle in molecule.angles() {
                let (i, j, k) = (angle.i(), angle.j(), angle.k());
                let (theta, d1, d2, d3) = system.angle_and_derivatives(i, j, k);
                for potential in system.angle_potentials(i, j, k) {
                    let force = potential.force(theta);
                    res[i] += force * d1;
                    res[j] += force * d2;
                    res[k] += force * d3;
                }
            }

            for dihedral in molecule.dihedrals() {
                let (i, j, k, m) = (dihedral.i(), dihedral.j(), dihedral.k(), dihedral.m());
                let (phi, d1, d2, d3, d4) = system.dihedral_and_derivatives(i, j, k, m);
                for potential in system.dihedral_potentials(i, j, k, m) {
                    let force = potential.force(phi);
                    res[i] += force * d1;
                    res[j] += force * d2;
                    res[k] += force * d3;
                    res[m] += force * d4;
                }
            }
        }

        if let Some(coulomb) = system.interactions().coulomb() {
            let forces = coulomb.borrow_mut().forces(system);
            debug_assert!(forces.len() == natoms, "Wrong `forces` size in coulomb potentials");
            for (i, force) in forces.iter().enumerate() {
                res[i] += force;
            }
        }

        for global in system.interactions().globals() {
            let forces = global.borrow_mut().forces(system);
            debug_assert!(forces.len() == natoms, "Wrong `forces` size in global potentials");
            for (i, force) in forces.iter().enumerate() {
                res[i] += force;
            }
        }
        return res;
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
/// Compute the instananeous temperature of the system
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
        let mut res = Matrix3::zero();
        for i in 0..system.size() {
            for j in (i+1)..system.size() {
                for potential in system.pair_potentials(i, j) {
                    let info = potential.restriction().informations(system, i, j);
                    if !info.excluded {
                        let d = system.nearest_image(i, j);
                        res += info.scaling * potential.virial(&d);
                    }
                }
            }
        }

        // FIXME: implement virial computations for molecular potentials
        // (angles & dihedrals)

        if let Some(coulomb) = system.interactions().coulomb() {
            res += coulomb.borrow_mut().virial(system);
        }

        for global in system.interactions().globals() {
            res += global.borrow_mut().virial(system);
        }

        return res;
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
        assert!(self.temperature >= 0.0);
        let virial_tensor = system.virial();
        let virial = virial_tensor.trace();
        let volume = system.volume();
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
        let virial = system.virial();
        let volume = system.volume();
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
        let mut kinetic = Matrix3::zero();
        for particle in system.iter() {
            let velocity = &particle.velocity;
            kinetic += particle.mass * velocity.tensorial(velocity);
        }

        let virial = Virial.compute(system);
        let volume = Volume.compute(system);
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
    use system::{System, Particle, UnitCell};
    use system::{InitVelocities, BoltzmannVelocities};
    use potentials::{Harmonic, NullPotential, PairInteraction};
    use constants::K_BOLTZMANN;
    use utils::unit_from;

    const EPS : f64 = 1e-8;

    fn test_pairs_system() -> System {
        let mut system = System::from_cell(UnitCell::cubic(10.0));;
        system.add_particle(Particle::new("F"));
        system[0].position = Vector3D::zero();
        system.add_particle(Particle::new("F"));
        system[1].position = Vector3D::new(1.3, 0.0, 0.0);

        let lj = Box::new(Harmonic{
            k: unit_from(300.0, "kJ/mol/A^2"),
            x0: unit_from(1.2, "A")
        });
        system.interactions_mut().add_pair("F", "F",
            PairInteraction::new(lj, 5.0)
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

        return system;
    }

    #[test]
    fn forces_pairs() {
        let system = &test_pairs_system();
        let res = Forces.compute(system);

        let forces_tot = res[0] + res[1];
        assert_eq!(forces_tot, Vector3D::zero());

        let force = unit_from(30.0, "kJ/mol/A");
        assert_approx_eq!(res[0][0], force, EPS);
        assert_approx_eq!(res[0][1], 0.0, EPS);
        assert_approx_eq!(res[0][1], 0.0, EPS);

        assert_approx_eq!(res[1][0], -force, EPS);
        assert_approx_eq!(res[1][1], 0.0, EPS);
        assert_approx_eq!(res[1][1], 0.0, EPS);
    }

    #[test]
    fn force_molecular() {
        let system = test_molecular_system();
        let res = Forces.compute(&system);
        let forces_tot = res[0] + res[1] + res[2] + res[3];
        assert_approx_eq!(forces_tot.norm2(), 0.0, 1e-12);
    }

    #[test]
    fn energy_pairs() {
        let system = &test_pairs_system();
        let kinetic = KineticEnergy.compute(system);
        let potential = PotentialEnergy.compute(system);
        let total = TotalEnergy.compute(system);

        assert_eq!(kinetic + potential, total);
        assert_eq!(kinetic, 0.0007483016557453699);
        assert_approx_eq!(potential, 1.5e-4, EPS);

        assert_eq!(kinetic, system.kinetic_energy());
        assert_eq!(potential, system.potential_energy());
        assert_eq!(total, system.total_energy());
    }

    #[test]
    fn energy_molecular() {
        let system = test_molecular_system();
        assert_approx_eq!(PotentialEnergy.compute(&system), unit_from(1800.0, "kJ/mol"));
    }

    #[test]
    fn temperature() {
        let system = &test_pairs_system();
        let temperature = Temperature.compute(system);
        assert_approx_eq!(temperature, 300.0, 1e-9);
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
    fn virial() {
        let system = &test_pairs_system();
        let virial = Virial.compute(system);

        let mut res = Matrix3::zero();
        let force = unit_from(30.0, "kJ/mol/A");
        res[(0, 0)] = - force * 1.3;

        for i in 0..3 {
            for j in 0..3 {
                assert_approx_eq!(virial[(i, j)], res[(i, j)], 1e-9);
            }
        }
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
    fn pressure_at_temperature() {
        let system = &test_pairs_system();
        let pressure = PressureAtTemperature{temperature: 550.0}.compute(system);

        let force = unit_from(30.0, "kJ/mol/A");
        let virial = -force * 1.3;
        let natoms = 2.0;
        let temperature = 550.0;
        let volume = 1000.0;
        let expected = natoms * K_BOLTZMANN * temperature / volume + virial / (3.0 * volume);

        assert_approx_eq!(pressure, expected, 1e-9);
        assert_eq!(pressure, system.pressure_at_temperature(550.0));

        let pressure = PressureAtTemperature{temperature: system.temperature()};
        let exact = pressure.compute(system);
        assert_eq!(exact, system.pressure());
    }

    #[test]
    #[should_panic]
    fn stress_at_temperature_negative_temperature() {
        let system = &test_pairs_system();
        let stress = StressAtTemperature{temperature: -4.0};
        let _ = stress.compute(system);
    }

    #[test]
    fn stress_at_temperature() {
        let system = &test_pairs_system();
        let stress = StressAtTemperature{temperature: 550.0}.compute(system);
        let pressure = PressureAtTemperature{temperature: 550.0}.compute(system);

        let trace = (stress[(0, 0)] + stress[(1, 1)] + stress[(2, 2)]) / 3.0;
        assert_approx_eq!(trace, pressure, 1e-9);
        assert_eq!(stress, system.stress_at_temperature(550.0));
    }

    #[test]
    fn stress() {
        let system = &test_pairs_system();
        let stress = Stress.compute(system);
        let pressure = Pressure.compute(system);

        let trace = stress.trace() / 3.0;
        assert_approx_eq!(trace, pressure, 1e-9);
        assert_eq!(stress, system.stress());
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

        assert_approx_eq!(pressure, expected, 1e-9);
        assert_eq!(pressure, system.pressure());
    }
}
