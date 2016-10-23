// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

use energy::{Potential, PairPotential, PairRestriction};

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
    potential: Box<PairPotential>,
    /// The associated pair restrictions
    restriction: PairRestriction,
    /// The cutoff distance
    cutoff: f64,
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
    /// # use lumol::energy::{Potential, PairInteraction};
    /// use lumol::energy::Harmonic;
    ///
    /// let potential = Box::new(Harmonic{x0: 0.5, k: 4.2});
    /// let interaction = PairInteraction::new(potential, 2.0);
    ///
    /// assert_eq!(interaction.energy(1.0), 0.525);
    /// // energy at and after the cutoff is zero
    /// assert_eq!(interaction.energy(2.0), 0.0);
    /// ```
    pub fn new(potential: Box<PairPotential>, cutoff: f64) -> PairInteraction {
        PairInteraction {
            potential: potential,
            cutoff: cutoff,
            restriction: PairRestriction::None,
            computation: PairComputation::Cutoff,
        }
    }

    /// Create a new `PairInteraction` with the given `cutoff`, using shifted
    /// computation of the energy.
    ///
    /// # Examples
    ///
    /// ```
    /// # use lumol::energy::{Potential, PairInteraction};
    /// use lumol::energy::Harmonic;
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
    pub fn shifted(potential: Box<PairPotential>, cutoff: f64) -> PairInteraction {
        let shift = potential.energy(cutoff);
        PairInteraction {
            potential: potential,
            cutoff: cutoff,
            restriction: PairRestriction::None,
            computation: PairComputation::Shifted(shift),
        }
    }

    /// Get the associated pair restriction. The default is to have no pair
    /// restriction.
    ///
    /// # Examples
    ///
    /// ```
    /// # use lumol::energy::PairInteraction;
    /// use lumol::energy::{NullPotential, PairRestriction};
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
    /// # use lumol::energy::PairInteraction;
    /// use lumol::energy::{NullPotential, PairRestriction};
    /// let mut interaction = PairInteraction::new(Box::new(NullPotential), 2.0);
    ///
    /// assert_eq!(interaction.restriction(), PairRestriction::None);
    /// interaction.set_restriction(PairRestriction::Exclude12);
    /// assert_eq!(interaction.restriction(), PairRestriction::Exclude12);
    /// ```
    pub fn set_restriction(&mut self, restriction: PairRestriction) {
        self.restriction = restriction;
    }
}

impl Potential for PairInteraction {
    fn energy(&self, r: f64) -> f64 {
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

    fn force(&self, r: f64) -> f64 {
        if r >= self.cutoff {
            0.0
        } else {
            self.potential.force(r)
        }
    }
}

impl PairPotential for PairInteraction {}

#[cfg(test)]
mod tests {
    use super::*;
    use energy::{NullPotential, LennardJones, PairRestriction};
    use energy::Potential;

    #[test]
    fn restriction() {
        let mut pairs = PairInteraction::new(Box::new(NullPotential), 10.0);
        assert_eq!(pairs.restriction(), PairRestriction::None);

        pairs.set_restriction(PairRestriction::Exclude13);
        assert_eq!(pairs.restriction(), PairRestriction::Exclude13);
    }

    #[test]
    fn cutoff() {
        let lj = LennardJones{sigma: 1.0, epsilon: 2.0};
        let pairs = PairInteraction::new(Box::new(lj), 4.0);

        assert_eq!(pairs.force(2.5), lj.force(2.5));
        assert_eq!(pairs.energy(2.5), lj.energy(2.5));

        assert_eq!(pairs.force(4.1), 0.0);
        assert_eq!(pairs.energy(4.1), 0.0);
    }

    #[test]
    fn shifted() {
        let lj = LennardJones{sigma: 1.0, epsilon: 2.0};
        let pairs = PairInteraction::shifted(Box::new(lj), 4.0);

        assert_eq!(pairs.force(2.5), lj.force(2.5));
        assert_approx_eq!(pairs.energy(2.5), -0.030681134109158216);

        assert_eq!(pairs.force(4.1), 0.0);
        assert_eq!(pairs.energy(4.1), 0.0);
    }
}
