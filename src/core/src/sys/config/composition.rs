// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

use std::ops::{Index, IndexMut};
use std::iter::{Map, Enumerate};
use std::vec;
use std::slice;

use sys::ParticleKind;

/// The system composition contains the number of particles of each kind
/// in the system.
///
/// # Examples
/// ```
/// # use lumol_core::sys::{Composition, ParticleKind};
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
///
/// for (kind, number) in &composition {
///     println!("We have {} particles of kind {}", number, kind);
/// }
/// ```
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Composition(Vec<usize>);

impl Composition {
    /// Create a new empty composition
    ///
    /// # Examples
    /// ```
    /// # use lumol_core::sys::Composition;
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
    /// # use lumol_core::sys::Composition;
    /// let mut composition = Composition::new();
    /// assert_eq!(composition.len(), 0);
    ///
    /// composition.resize(10);
    /// assert_eq!(composition.len(), 10);
    /// ```
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// Check if the composition is empty, *i.e.* if it contains no particle
    /// kind. A composition with only zero entries (`0 => 0, 2 => 0, 3 => 0`)
    /// is not empty.
    ///
    /// # Examples
    /// ```
    /// # use lumol_core::sys::Composition;
    /// let mut composition = Composition::new();
    /// assert_eq!(composition.len(), 0);
    /// assert!(composition.is_empty());
    /// ```
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Resize the composition to hold `size` items. The new particles kinds
    /// start with no associated particles.
    ///
    /// # Examples
    /// ```
    /// # use lumol_core::sys::{Composition, ParticleKind};
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

    /// Get an iterator over `(kind, number of particles)`
    ///
    /// # Examples
    /// ```
    /// # use lumol_core::sys::{Composition, ParticleKind};
    /// let mut composition = Composition::new();
    /// composition.resize(2);
    /// composition[ParticleKind(0)] = 10;
    /// composition[ParticleKind(1)] = 55;
    ///
    /// for (kind, number) in composition.iter() {
    ///     println!("{}: {}", kind, number);
    /// }
    ///
    /// let mut iter = composition.iter();
    /// assert_eq!(iter.next(), Some((ParticleKind(0), &10)));
    /// assert_eq!(iter.next(), Some((ParticleKind(1), &55)));
    /// assert_eq!(iter.next(), None);
    /// ```
    pub fn iter(&self) -> Iter {
        self.into_iter()
    }

    /// Get an iterator over `(kind, mutable number of particles)`
    ///
    /// # Examples
    /// ```
    /// # use lumol_core::sys::{Composition, ParticleKind};
    /// let mut composition = Composition::new();
    /// composition.resize(2);
    /// composition[ParticleKind(0)] = 10;
    /// composition[ParticleKind(1)] = 55;
    ///
    /// for (kind, number) in composition.iter_mut() {
    ///     if kind == ParticleKind(0) {
    ///         *number += 10;
    ///     }
    /// }
    ///
    /// let mut iter = composition.iter_mut();
    /// assert_eq!(iter.next(), Some((ParticleKind(0), &mut 20)));
    /// assert_eq!(iter.next(), Some((ParticleKind(1), &mut 55)));
    /// assert_eq!(iter.next(), None);
    /// ```
    pub fn iter_mut(&mut self) -> IterMut {
        self.into_iter()
    }
}

pub struct IntoIter {
    iter: Map<Enumerate<vec::IntoIter<usize>>, fn((usize, usize)) -> (ParticleKind, usize)>
}

impl Iterator for IntoIter {
    type Item = (ParticleKind, usize);
    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next()
    }
}

impl IntoIterator for Composition {
    type Item = (ParticleKind, usize);
    type IntoIter = IntoIter;
    fn into_iter(self) -> Self::IntoIter {
        fn convert((i, ni): (usize, usize)) -> (ParticleKind, usize) {
            (ParticleKind(i as u32), ni)
        }

        IntoIter {
            iter: self.0.into_iter().enumerate().map(convert)
        }
    }
}

pub struct Iter<'a> {
    iter: Map<Enumerate<slice::Iter<'a, usize>>, fn((usize, &'a usize)) -> (ParticleKind, &'a usize)>
}

impl<'a> Iterator for Iter<'a> {
    type Item = (ParticleKind, &'a usize);
    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next()
    }
}

impl<'a> IntoIterator for &'a Composition {
    type Item = (ParticleKind, &'a usize);
    type IntoIter = Iter<'a>;
    fn into_iter(self) -> Self::IntoIter {
        fn convert((i, ni): (usize, &usize)) -> (ParticleKind, &usize) {
            (ParticleKind(i as u32), ni)
        }

        Iter {
            iter: self.0.iter().enumerate().map(convert)
        }
    }
}

pub struct IterMut<'a> {
    iter: Map<Enumerate<slice::IterMut<'a, usize>>, fn((usize, &'a mut usize)) -> (ParticleKind, &'a mut usize)>
}

impl<'a> Iterator for IterMut<'a> {
    type Item = (ParticleKind, &'a mut usize);
    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next()
    }
}

impl<'a> IntoIterator for &'a mut Composition {
    type Item = (ParticleKind, &'a mut usize);
    type IntoIter = IterMut<'a>;
    fn into_iter(self) -> Self::IntoIter {
        fn convert((i, ni): (usize, &mut usize)) -> (ParticleKind, &mut usize) {
            (ParticleKind(i as u32), ni)
        }

        IterMut {
            iter: self.0.iter_mut().enumerate().map(convert)
        }
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
        assert!(composition.is_empty());

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

    #[test]
    fn iter() {
        let mut composition = Composition::new();
        composition.resize(4);
        composition[ParticleKind(0)] = 2;
        composition[ParticleKind(1)] = 4;
        composition[ParticleKind(2)] = 6;
        composition[ParticleKind(3)] = 8;

        let mut expected = 2;
        for (_, &count) in &composition {
            assert_eq!(count, expected);
            expected += 2;
        }

        for (_, count) in &mut composition {
            *count *= 3;
        }

        for (kind, count) in composition {
            match kind {
                ParticleKind(0) => assert_eq!(count, 6),
                ParticleKind(1) => assert_eq!(count, 12),
                ParticleKind(2) => assert_eq!(count, 18),
                ParticleKind(3) => assert_eq!(count, 24),
                _ => panic!(),
            }
        }
    }
}
