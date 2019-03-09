// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

//! Global potential are potentials acting on the whole system at once
//!
//! They can be coulombic potentials, or external provided potential function
//! for example.
use crate::PairRestriction;
use crate::Configuration;
use crate::{Matrix3, Vector3D};

/// A potential acting on the whole [System][System] at once.
///
/// This kind of potential is passed all the data in a system and returns the
/// energy/forces/virial contribution for this potential. It has three main
/// usages:
///
/// - Potentials that does not fit the [pair/bond/angle/dihedral
///   framework][Potential] for potentials: the coulombic potential solver is
///   the more important example;
/// - External potentials like an external electric field or a hard wall;
/// - Potentials coming from external software: DFT computations, metadynamic
///   drivers, or other.
///
/// [System]: ../sys/struct.System.html
/// [Potential]: index.html#potentials
///
/// # Examples
///
/// ```
/// use lumol_core::energy::{GlobalPotential, GlobalCache};
/// use lumol_core::types::{Vector3D, Matrix3};
/// use lumol_core::sys::{System, Configuration, Particle, Molecule, UnitCell};
///
/// /// Shift the energy of all the particles by a given delta.
/// #[derive(Clone)]
/// struct ShiftAll {
///     delta: f64,
/// }
///
/// impl GlobalPotential for ShiftAll {
///     fn cutoff(&self) -> Option<f64> {
///         None
///     }
///
///     fn energy(&self, configuration: &Configuration) -> f64 {
///         // shift all particles by delta
///         self.delta * configuration.size() as f64
///     }
///
///     fn forces(&self, configuration: &Configuration, forces: &mut [Vector3D]) {
///         // this potential does not changes the forces
///     }
///
///     fn atomic_virial(&self, configuration: &Configuration) -> Matrix3 {
///         // the virial is null as all the forces are null
///         Matrix3::zero()
///     }
/// }
///
/// // Not implementing `GlobalCache` means that we will not be able to use
/// // `ShiftAll` in Monte Carlo simulations.
/// impl GlobalCache for ShiftAll {
///     fn move_molecule_cost(&self, _: &Configuration, _: usize, _: &[Vector3D]) -> f64 {
///         unimplemented!()
///     }
///
///     fn update(&self) {
///         unimplemented!()
///     }
/// }
///
/// // A simple test
/// let mut system = System::with_cell(UnitCell::cubic(10.0));
/// system.add_molecule(Molecule::new(Particle::new("Ar")));
/// system.add_molecule(Molecule::new(Particle::new("Ar")));
///
/// system.add_global_potential(Box::new(ShiftAll{delta: 1.0}));
///
/// assert_eq!(system.potential_energy(), 2.0);
/// assert_eq!(system.forces(), vec![Vector3D::zero(); 2]);
/// assert_eq!(system.virial(), Matrix3::zero());
/// ```
pub trait GlobalPotential: GlobalCache + BoxCloneGlobal + Send + Sync {
    /// Return the cut off radius.
    fn cutoff(&self) -> Option<f64>;

    /// Compute the energetic contribution of this potential
    fn energy(&self, configuration: &Configuration) -> f64;

    /// Compute the force contribution of this potential. This function should
    /// return a vector containing the force acting on each particle in the
    /// configuration.
    fn forces(&self, configuration: &Configuration, forces: &mut [Vector3D]);

    /// Compute the total virial contribution of this potential, using the
    /// atomic virial definition
    fn atomic_virial(&self, configuration: &Configuration) -> Matrix3;

    /// Compute the total virial contribution of this potential, using the
    /// molecular virial definition. This default to `atomic_virial`.
    fn molecular_virial(&self, configuration: &Configuration) -> Matrix3 {
        return self.atomic_virial(configuration);
    }
}

impl_box_clone!(GlobalPotential, BoxCloneGlobal, box_clone_gobal);

/// Energetic cache for global potentials.
///
/// This trait provide all the functions needed by [`EnergyCache`][EnergyCache]
/// to compute partial energy updates in Monte Carlo simulations. You can use
/// a `panic!`ing implementation for all methods if you never need to use a
/// given [`GlobalPotential`][GlobalPotential] in Monte Carlo simulations.
///
/// All methods take a non-mutable `&self` receiver, which means you may want
/// to wrap the implemntation in `RwLock` or `Mutex` to allow for inner
/// mutability while still implementing `Send + Sync`.
///
/// [EnergyCache]: ../sys/struct.EnergyCache.html
/// [GlobalPotential]: trait.GlobalPotential.html
///
/// # Examples
///
/// ```
/// use lumol_core::energy::{GlobalPotential, GlobalCache};
/// use lumol_core::types::{Vector3D, Matrix3};
/// use lumol_core::sys::Configuration;
///
/// /// Shift the energy of all the particles by a given delta.
/// #[derive(Clone)]
/// struct ShiftAll {
///     delta: f64,
/// }
///
/// impl GlobalPotential for ShiftAll {
///     fn cutoff(&self) -> Option<f64> {
///         None
///     }
///
///     fn energy(&self, configuration: &Configuration) -> f64 {
///         // shift all particles by delta
///         self.delta * configuration.size() as f64
///     }
///
///     fn forces(&self, _: &Configuration, _: &mut [Vector3D]) {
///         // this potential does not changes the forces
///     }
///
///     fn atomic_virial(&self, _: &Configuration) -> Matrix3 {
///         // the virial is null as all the forces are null
///         Matrix3::zero()
///     }
/// }
///
/// impl GlobalCache for ShiftAll {
///     fn move_molecule_cost(&self, _: &Configuration, _: usize, _: &[Vector3D]) -> f64 {
///         // The cost of moving particles is null, because all the particles
///         // get the same energy shift whatever there position are.
///         return 0.0
///     }
///
///     fn update(&self) {
///         // We are not storing anything in the ShiftAll struct, so this
///         // function is a no-op.
///     }
/// }
/// ```
pub trait GlobalCache {
    /// Get the cost of moving a rigid molecule in the system.
    ///
    /// This function is passed the current `configuration`, the index of the
    /// molecule in the configuration; and the `new_positions` of the
    /// particles. The previous positions of the particles are still in the
    /// system.
    fn move_molecule_cost(
        &self,
        configuration: &Configuration,
        molecule_id: usize,
        new_positions: &[Vector3D],
    ) -> f64;

    /// Update the cache as needed after a call to `move_molecule_cost`.
    ///
    /// If the Monte Carlo move is accepted, this function will be called and
    /// should update any cached quantity so that further call to
    /// `GlobalPotential::energy` gives the right value.
    fn update(&self);
}

/// Electrostatic potential solver.
///
/// This trait is a marker trait for [global potentials][GlobalPotential] that
/// are actually coulombic potential solvers.
///
/// [GlobalPotential]: trait.GlobalPotential.html
pub trait CoulombicPotential: GlobalPotential + BoxCloneCoulombic {
    /// Set the pair restriction scheme to use to the given `restriction`. All
    /// future call to `GlobalPotential::energy`, `GlobalPotential::force` or
    /// `GlobalPotential::virial` should use this restriction.
    fn set_restriction(&mut self, restriction: PairRestriction);
}

impl_box_clone!(CoulombicPotential, BoxCloneCoulombic, box_clone_coulombic);

mod wolf;
pub use self::wolf::Wolf;

mod ewald;
pub use self::ewald::{Ewald, SharedEwald};
