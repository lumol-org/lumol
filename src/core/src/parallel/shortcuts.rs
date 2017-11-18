// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

use std::ops::Range;

use ndarray::Zip;
use ndarray_parallel::NdarrayIntoParallelIterator;
use rayon::iter::Map;
use rayon::prelude::{IntoParallelIterator, ParallelIterator};


/// Utility trait that adds shortcuts for `IntoParallelIterator` structs.
///
/// # Example
///
/// ```
/// use lumol_core::parallel::prelude::*;
///
/// let s = (0..100_i32).par_map(|i| -i).sum();
/// assert_eq!(-4950, s);
/// ```
pub trait ParallelShortcuts: Sized {
    /// The iterator type
    type Iter: ParallelIterator;

    /// Shortcut for `into_par_iter().map()`
    fn par_map<F, R>(self, map_op: F) -> Map<Self::Iter, F>
    where
        F: Fn(<Self::Iter as ParallelIterator>::Item) -> R + Sync + Send,
        R: Send;
}

impl<T> ParallelShortcuts for Range<T>
where
    Range<T>: IntoParallelIterator,
{
    type Iter = <Range<T> as IntoParallelIterator>::Iter;

    fn par_map<F, R>(self, map_op: F) -> Map<Self::Iter, F>
    where
        F: Fn(<Self::Iter as ParallelIterator>::Item) -> R + Sync + Send,
        R: Send,
    {
        self.into_par_iter().map(map_op)
    }
}

impl<Parts, D> ParallelShortcuts for Zip<Parts, D>
where
    Zip<Parts, D>: NdarrayIntoParallelIterator,
{
    type Iter = <Zip<Parts, D> as NdarrayIntoParallelIterator>::Iter;

    fn par_map<F, R>(self, map_op: F) -> Map<Self::Iter, F>
    where
        F: Fn(<Self::Iter as ParallelIterator>::Item) -> R + Sync + Send,
        R: Send,
    {
        self.into_par_iter().map(map_op)
    }
}
