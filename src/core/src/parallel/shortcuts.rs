// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

use std::ops::Range;

use rayon::prelude::{ParallelIterator, IntoParallelIterator};
use rayon::iter::{Map, MapFn};
use ndarray_parallel::NdarrayIntoParallelIterator;
use ndarray::Zip;


/// Utility trait that adds shortcuts for `IntoParallelIterator` structs.
///
/// # Example
///
/// ```
/// use lumol::parallel::prelude::*;
///
/// let s = (0..100_i32).par_map(|i| -i).sum();
/// assert_eq!(-4950, s);
/// ```
pub trait ParallelShortcuts: Sized {
    /// The iterator type
    type Iter: ParallelIterator;

    /// Dummy function to make it work
    fn _into_par_iter(self) -> Self::Iter;

    /// Shortcut for `into_par_iter().map()`
    fn par_map<F, R>(self, map_op: F) -> Map<Self::Iter, MapFn<F>>
        where F: Fn(<Self::Iter as ParallelIterator>::Item) -> R + Sync, R: Send {
        self._into_par_iter().map(map_op)
    }
}

impl<T> ParallelShortcuts for Range<T>
    where Range<T>: IntoParallelIterator {
    type Iter = <Range<T> as IntoParallelIterator>::Iter;

    fn _into_par_iter(self) -> Self::Iter {
        self.into_par_iter()
    }
}

impl<Parts, D> ParallelShortcuts for Zip<Parts, D>
    where Zip<Parts, D>: NdarrayIntoParallelIterator {
    type Iter = <Zip<Parts, D> as NdarrayIntoParallelIterator>::Iter;

    fn _into_par_iter(self) -> Self::Iter {
        self.into_par_iter()
    }
}
