/* Cymbalum, Molecular Simulation in Rust - Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
 */

//! This module provide type for counting mutable borrow of a value. The `Bc<T>`
//! type is a small wrapper on top of a value of type `T` which count the number
//! of time the value has been mutably borrowed since it's creation.

use std::ops::{Deref, DerefMut};

/// The borrow counter struct for type `T`.
pub struct Bc<T> {
    counter: usize,
    val: T
}

impl<T> Bc<T> {
    /// Create a new `Bc<T>` containing the value `val`.
    pub fn new(val: T) -> Bc<T> {
        Bc {
            val: val,
            counter: 0,
        }
    }

    /// Reset the borrow counter
    pub fn reset(&mut self) {
        self.counter = 0;
    }

    /// Get the number of time this structure has been mutably borrowed.
    pub fn count(&self) -> usize {
        self.counter
    }
}

impl<T> Deref for Bc<T> {
    type Target = T;
    fn deref(&self) -> &T {
        &self.val
    }
}

impl<T> DerefMut for Bc<T> {
     fn deref_mut(&mut self) -> &mut T {
         self.counter = self.counter.wrapping_add(1);
         &mut self.val
     }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::usize;

    fn do_nothing(_: &mut f64) {}

    #[test]
    fn count() {
        let mut a = Bc::new(5.0);
        assert_eq!(a.count(), 0);

        *a = 89.0;
        assert_eq!(a.count(), 1);

        do_nothing(&mut a);
        assert_eq!(a.count(), 2);
    }

    #[test]
    fn reset() {
        let mut a = Bc::new(3);

        assert_eq!(a.count(), 0);

        *a = 18;
        *a = 42;
        assert_eq!(a.count(), 2);

        a.reset();
        assert_eq!(a.count(), 0);
    }

    #[test]
    fn overflow() {
        let mut a = Bc::new(3);
        a.counter = usize::MAX - 1;

        *a = 18;
        *a = 18;
        assert_eq!(a.count(), 0);
    }

    #[test]
    fn non_mutable() {
        fn observe(_: &f64) {/* Do nothing */}

        let a = Bc::new(3.0);
        assert_eq!(a.count(), 0);

        observe(&a);
        observe(&a);
        observe(&a);
        assert_eq!(a.count(), 0);
    }
}
