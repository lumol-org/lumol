// Cymbalum, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

use potentials::{Potential, PairPotential, PairRestriction};

/// The different way to compute non-bonded pair interactions
#[derive(Clone, Copy, Debug)]
enum PairComputation {
    /// Using only a cutoff distance
    Cutoff,
    /// Using a cutoff distance and a shift
    Shifted(f64),
}

/// A non-bonded interaction between two particle. The potential `T` is
/// associated with a [pair restriction][PairRestriction]. For all non-bonded
/// pairs, the potential is computed up to a cutoff distance. An additional
/// shifting of the potential can be used in molecular dynamics, to ensure that
/// the energy is continuous at the cutoff distance.
///
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
    /// Create a new `PairInteraction` with the given cutoff
    pub fn new(potential: Box<PairPotential>, cutoff: f64) -> PairInteraction {
        PairInteraction {
            potential: potential,
            cutoff: cutoff,
            restriction: PairRestriction::None,
            computation: PairComputation::Cutoff,
        }
    }

    /// Create a new `PairInteraction` with the given cutoff, using shifted
    /// computation of the energy.
    pub fn shifted(potential: Box<PairPotential>, cutoff: f64) -> PairInteraction {
        let shift = potential.energy(cutoff);
        PairInteraction {
            potential: potential,
            cutoff: cutoff,
            restriction: PairRestriction::None,
            computation: PairComputation::Shifted(shift),
        }
    }

    /// Get the associated pair restriction
    pub fn restriction(&self) -> PairRestriction {
        self.restriction
    }

    /// Set the restriction associated with this potential
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

#[cfg(test)]
mod tests {
    use super::*;
    use potentials::{NullPotential, LennardJones, PairRestriction};
    use potentials::Potential;

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
