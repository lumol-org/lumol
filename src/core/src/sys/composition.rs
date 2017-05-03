// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 Lumol's contributors — BSD license

use std::ops::{Index, IndexMut};
use sys::ParticleKind;

/// The system composition contains the number of particles of each kind
/// in the system.
///
/// # Examples
/// ```
/// # use lumol::sys::{Composition, ParticleKind};
/// let mut composition = Composition::new();
/// composition.resize(10);
///
/// composition[ParticleKind(2)] = 56;
/// composition[ParticleKind(8)] = 2;
/// composition[ParticleKind(5)] = 42;
///
/// assert_eq!(composition[ParticleKind(2)], 56);
/// assert_eq!(composition[ParticleKind(8)], 2);
/// assert_eq!(composition[ParticleKind(5)], 42);
/// ```
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Composition(Vec<usize>);

impl Composition {
    /// Create a new empty composition
    ///
    /// # Examples
    /// ```
    /// # use lumol::sys::Composition;
    /// let composition = Composition::new();
    /// assert_eq!(composition.len(), 0);
    /// ```
    pub fn new() -> Composition {
        Composition(Vec::new())
    }

    /// Get the size of the composition, *i.e.* the number of different
    /// particle kinds in the composition.
    ///
    /// # Examples
    /// ```
    /// # use lumol::sys::Composition;
    /// let mut composition = Composition::new();
    /// assert_eq!(composition.len(), 0);
    ///
    /// composition.resize(10);
    /// assert_eq!(composition.len(), 10);
    /// ```
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// Resize the composition to hold `size` items. The new particles kinds
    /// start with no associated particles.
    ///
    /// # Examples
    /// ```
    /// # use lumol::sys::{Composition, ParticleKind};
    /// let mut composition = Composition::new();
    /// assert_eq!(composition.len(), 0);
    ///
    /// composition.resize(10);
    /// assert_eq!(composition.len(), 10);
    /// assert_eq!(composition[ParticleKind(8)], 0);
    /// ```
    pub fn resize(&mut self, size: usize) {
        self.0.resize(size, 0)
    }
}


impl Index<ParticleKind> for Composition {
    type Output = usize;

    #[inline]
    fn index(&self, i: ParticleKind) -> &usize {
        &self.0[i.0 as usize]
    }
}

impl IndexMut<ParticleKind> for Composition {
    #[inline]
    fn index_mut(&mut self, i: ParticleKind) -> &mut usize {
        &mut self.0[i.0 as usize]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn len() {
        let mut composition = Composition::new();
        assert_eq!(composition.len(), 0);

        composition.resize(10);
        assert_eq!(composition.len(), 10);
    }

    #[test]
    fn index() {
        let mut composition = Composition::new();
        composition.resize(10);

        composition[ParticleKind(2)] = 56;
        composition[ParticleKind(8)] = 2;
        composition[ParticleKind(5)] = 42;

        assert_eq!(composition[ParticleKind(2)], 56);
        assert_eq!(composition[ParticleKind(8)], 2);
        assert_eq!(composition[ParticleKind(5)], 42);
    }
}
