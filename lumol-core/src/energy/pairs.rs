// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

use crate::{PairPotential, PairRestriction};
use crate::{Matrix3, Vector3D};

/// The different way to compute non-bonded pair interactions
#[derive(Clone, Copy, Debug)]
enum PairComputation {
    /// Using only a cutoff distance
    Cutoff,
    /// Using a cutoff distance and a shift
    Shifted(f64),
}

/// A non-bonded interaction between two particle.
///
/// This is a thin wrapper around a [`Box<PairPotential>`][PairPotential]
/// associated with a [pair restriction][PairRestriction]. It ensure that the
/// potential is computed up to a cutoff distance. An additional shifting of the
/// potential can be used in molecular dynamics, to ensure that the energy is
/// continuous at the cutoff distance.
///
/// [PairPotential]: trait.PairPotential.html
/// [PairRestriction]: enum.PairRestriction.html
#[derive(Clone)]
pub struct PairInteraction {
    /// The potential of this interaction
    potential: Box<dyn PairPotential>,
    /// The associated pair restrictions
    restriction: PairRestriction,
    /// The cutoff distance
    cutoff: f64,
    /// Should we use tail corrections
    tail: bool,
    /// The computation mode
    computation: PairComputation,
}

impl PairInteraction {
    /// Create a new `PairInteraction` for the given `potential` and using the
    /// given `cutoff` distance.
    ///
    /// # Examples
    ///
    /// ```
    /// use lumol_core::energy::PairInteraction;
    /// use lumol_core::energy::Harmonic;
    ///
    /// let potential = Box::new(Harmonic{x0: 0.5, k: 4.2});
    /// let interaction = PairInteraction::new(potential, 2.0);
    ///
    /// assert_eq!(interaction.energy(1.0), 0.525);
    /// // energy at and after the cutoff is zero
    /// assert_eq!(interaction.energy(2.0), 0.0);
    /// ```
    pub fn new(potential: Box<dyn PairPotential>, cutoff: f64) -> PairInteraction {
        PairInteraction {
            potential: potential,
            cutoff: cutoff,
            restriction: PairRestriction::None,
            computation: PairComputation::Cutoff,
            tail: false,
        }
    }

    /// Create a new `PairInteraction` with the given `cutoff`, using shifted
    /// computation of the energy.
    ///
    /// # Examples
    ///
    /// ```
    /// use lumol_core::energy::PairInteraction;
    /// use lumol_core::energy::Harmonic;
    ///
    /// let potential = Box::new(Harmonic{x0: 0.5, k: 4.2});
    /// let interaction = PairInteraction::shifted(potential, 2.0);
    ///
    /// assert_eq!(interaction.energy(1.0), -4.2);
    ///
    /// // energy goes smoothly to zero at the cutoff
    /// assert!(interaction.energy(1.999).abs() < 0.01);
    /// // energy after the cutoff is zero
    /// assert_eq!(interaction.energy(2.0), 0.0);
    /// ```
    pub fn shifted(potential: Box<dyn PairPotential>, cutoff: f64) -> PairInteraction {
        let shift = potential.energy(cutoff);
        PairInteraction {
            potential: potential,
            cutoff: cutoff,
            restriction: PairRestriction::None,
            computation: PairComputation::Shifted(shift),
            tail: false,
        }
    }

    /// Enable the use of tail corrections for energy and virial contribution
    /// of this pair interaction.
    ///
    /// # Examples
    ///
    /// ```
    /// use lumol_core::energy::PairInteraction;
    /// use lumol_core::energy::LennardJones;
    ///
    /// let potential = Box::new(LennardJones{sigma: 0.5, epsilon: 4.2});
    /// let mut interaction = PairInteraction::new(potential, 2.0);
    ///
    /// assert_eq!(interaction.tail_energy(), 0.0);
    ///
    /// interaction.enable_tail_corrections();
    /// assert!(interaction.tail_energy() < 0.0);
    /// ```
    pub fn enable_tail_corrections(&mut self) {
        self.tail = true;
    }

    /// Get the associated pair restriction. The default is to have no pair
    /// restriction.
    ///
    /// # Examples
    ///
    /// ```
    /// use lumol_core::energy::PairInteraction;
    /// use lumol_core::energy::{NullPotential, PairRestriction};
    /// let interaction = PairInteraction::new(Box::new(NullPotential), 2.0);
    ///
    /// assert_eq!(interaction.restriction(), PairRestriction::None);
    /// ```
    pub fn restriction(&self) -> PairRestriction {
        self.restriction
    }

    /// Set the pair restriction associated with this interaction.
    ///
    /// # Examples
    ///
    /// ```
    /// use lumol_core::energy::PairInteraction;
    /// use lumol_core::energy::{NullPotential, PairRestriction};
    /// let mut interaction = PairInteraction::new(Box::new(NullPotential), 2.0);
    ///
    /// assert_eq!(interaction.restriction(), PairRestriction::None);
    /// interaction.set_restriction(PairRestriction::Exclude12);
    /// assert_eq!(interaction.restriction(), PairRestriction::Exclude12);
    /// ```
    pub fn set_restriction(&mut self, restriction: PairRestriction) {
        self.restriction = restriction;
    }

    /// Return the cutoff radius
    ///
    /// # Examples
    ///
    /// ```
    /// use lumol_core::energy::PairInteraction;
    /// use lumol_core::energy::LennardJones;
    ///
    /// let ar = LennardJones{sigma: 3.405, epsilon: 1.0};
    /// let interaction = PairInteraction::new(Box::new(ar), 9.1935);
    /// let rc = interaction.cutoff();
    /// assert_eq!(rc, 2.7f64 * 3.405);
    /// ```
    pub fn cutoff(&self) -> f64 {
        self.cutoff
    }
}

impl PairInteraction {
    /// Get the energy for this pair interaction at the distance `r`.
    ///
    /// # Examples
    ///
    /// ```
    /// use lumol_core::energy::PairInteraction;
    /// use lumol_core::energy::Harmonic;
    ///
    /// let potential = Box::new(Harmonic{x0: 0.5, k: 4.2});
    /// let interaction = PairInteraction::new(potential, 2.0);
    ///
    /// assert_eq!(interaction.energy(1.0), 0.525);
    /// // energy at and after the cutoff is zero
    /// assert_eq!(interaction.energy(2.0), 0.0);
    /// ```
    pub fn energy(&self, r: f64) -> f64 {
        if r >= self.cutoff {
            0.0
        } else {
            let energy = self.potential.energy(r);
            match self.computation {
                PairComputation::Cutoff => energy,
                PairComputation::Shifted(shift) => energy - shift,
            }
        }
    }

    /// Get the norm of the force for this pair interaction at the distance `r`.
    ///
    /// # Examples
    ///
    /// ```
    /// use lumol_core::energy::PairInteraction;
    /// use lumol_core::energy::Harmonic;
    ///
    /// let potential = Box::new(Harmonic{x0: 0.5, k: 4.2});
    /// let interaction = PairInteraction::new(potential, 2.0);
    ///
    /// assert_eq!(interaction.force(1.0), -2.1);
    /// // force at and after the cutoff is zero
    /// assert_eq!(interaction.force(2.0), 0.0);
    /// ```
    pub fn force(&self, r: f64) -> f64 {
        if r >= self.cutoff {
            0.0
        } else {
            self.potential.force(r)
        }
    }

    /// Get the virial contribution for this pair interaction at the distance
    /// `r`.
    ///
    /// # Examples
    ///
    /// ```
    /// use lumol_core::energy::PairInteraction;
    /// use lumol_core::energy::Harmonic;
    /// # use lumol_core::types::Vector3D;
    ///
    /// let potential = Box::new(Harmonic{x0: 0.5, k: 4.2});
    /// let interaction = PairInteraction::shifted(potential, 2.0);
    ///
    /// let r = Vector3D::new(1.0, 0.0, 0.3);
    /// let force = interaction.force(r.norm()) * r / r.norm();
    /// assert_eq!(interaction.virial(&r), r.tensorial(&force));
    /// ```
    pub fn virial(&self, r: &Vector3D) -> Matrix3 {
        if r.norm() >= self.cutoff {
            Matrix3::zero()
        } else {
            self.potential.virial(r)
        }
    }

    /// Get the tail correction to the energy for this pair interaction
    ///
    /// # Examples
    ///
    /// ```
    /// use lumol_core::energy::PairInteraction;
    /// use lumol_core::energy::LennardJones;
    ///
    /// let potential = Box::new(LennardJones{sigma: 0.5, epsilon: 4.2});
    /// let mut interaction = PairInteraction::new(potential, 2.0);
    /// interaction.enable_tail_corrections();
    ///
    /// assert_eq!(interaction.tail_energy(), -0.010936609903971353);
    /// ```
    pub fn tail_energy(&self) -> f64 {
        if self.tail {
            self.potential.tail_energy(self.cutoff)
        } else {
            0.0
        }
    }

    /// Get the tail correction to the virial for this pair interaction
    ///
    /// # Examples
    ///
    /// ```
    /// use lumol_core::energy::PairInteraction;
    /// use lumol_core::energy::LennardJones;
    /// # use lumol_core::types::Matrix3;
    ///
    /// let potential = Box::new(LennardJones{sigma: 0.5, epsilon: 4.2});
    /// let mut interaction = PairInteraction::new(potential, 2.0);
    /// interaction.enable_tail_corrections();
    ///
    /// let w = -0.02187143961588542;
    /// let virial = Matrix3::new([
    ///     [w, 0.0, 0.0],
    ///     [0.0, w, 0.0],
    ///     [0.0, 0.0, w],
    /// ]);
    ///
    /// assert_eq!(interaction.tail_virial(), virial);
    /// ```
    pub fn tail_virial(&self) -> Matrix3 {
        if self.tail {
            let tensor = Matrix3::one() / 3.0;
            return self.potential.tail_virial(self.cutoff) * tensor;
        }

        return Matrix3::zero();
    }
}

#[cfg(test)]
#[allow(clippy::unreadable_literal)]
mod tests {
    use super::*;
    use crate::{LennardJones, NullPotential, PairRestriction};
    use crate::Potential;

    use approx::assert_ulps_eq;

    #[test]
    fn restriction() {
        let mut pairs = PairInteraction::new(Box::new(NullPotential), 10.0);
        assert_eq!(pairs.restriction(), PairRestriction::None);

        pairs.set_restriction(PairRestriction::Exclude13);
        assert_eq!(pairs.restriction(), PairRestriction::Exclude13);
    }

    #[test]
    fn test_cutoff() {
        let lj = LennardJones {
            sigma: 1.0,
            epsilon: 2.0,
        };
        let pairs = PairInteraction::new(Box::new(lj), 4.0);

        assert_eq!(pairs.force(2.5), lj.force(2.5));
        assert_eq!(pairs.energy(2.5), lj.energy(2.5));

        assert_eq!(pairs.force(4.1), 0.0);
        assert_eq!(pairs.energy(4.1), 0.0);

        assert_eq!(pairs.cutoff(), 4.0);
    }

    #[test]
    fn shifted() {
        let lj = LennardJones {
            sigma: 1.0,
            epsilon: 2.0,
        };
        let pairs = PairInteraction::shifted(Box::new(lj), 4.0);

        assert_ulps_eq!(pairs.force(2.5), lj.force(2.5));
        assert_ulps_eq!(pairs.energy(2.5), -0.030681134109158216);

        assert_eq!(pairs.force(4.1), 0.0);
        assert_eq!(pairs.energy(4.1), 0.0);
    }

    #[test]
    fn tail_corrections() {
        let lj = LennardJones {
            sigma: 1.0,
            epsilon: 2.0,
        };
        let mut pairs = PairInteraction::shifted(Box::new(lj), 4.0);
        pairs.enable_tail_corrections();

        assert_eq!(pairs.tail_energy(), -0.041663275824652776);
        assert_ulps_eq!(pairs.tail_virial().trace(), -0.24995930989583334);
    }
}
