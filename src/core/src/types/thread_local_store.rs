use std::cell::RefCell;
use std::ops::AddAssign;

use thread_local::CachedThreadLocal;

/// A simple struct that delivers thread local variables.
///
/// Some computations need thread local copies of the memory
/// to avoid data races. This is a simple wrapper over
/// the `thread_local` crate that exposes a simple API
/// for doing so.
///
/// # Example
///
/// Computing in parallel how many times a remainder of
/// the division by 3 (so either 0, 1 or 2) appears in the
/// range `(0..100)`.
///
/// ```
/// extern crate rayon;
/// extern crate lumol;
///
/// use rayon::prelude::*;
/// use lumol::types::ThreadLocalStore;
///
/// fn main() {
///     let store = ThreadLocalStore::new(|| vec![0, 0, 0]);
///
///     (0..100_usize).into_par_iter().for_each(|i| {
///          let mut thread_local_vec = store.get().borrow_mut();
///          thread_local_vec[i % 3] += 1;
///     });
///
///     let mut result = vec![0, 0, 0];
///     store.add_local_values(&mut result);
///
///     assert_eq!(result, vec![34, 33, 33]);
/// }
///
/// ```
///
pub struct ThreadLocalStore<T, F>
    where T: Sized + Send, F: Fn() -> T {
    inner: CachedThreadLocal<RefCell<T>>,
    init: F
}

pub struct IntoIter<T: Sized + Send> {
    inner: <CachedThreadLocal<RefCell<T>> as IntoIterator>::IntoIter
}

impl<T: Sized + Send> Iterator for IntoIter<T> {
    type Item = T;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        self.inner.next().map(|rc| {
            rc.into_inner()
        })
    }
}

impl<T, F> IntoIterator for ThreadLocalStore<T, F>
    where T: Sized + Send, F: Fn() -> T {
    type Item = T;
    type IntoIter = IntoIter<T>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        IntoIter { inner: self.inner.into_iter() }
    }
}

impl<T, F> ThreadLocalStore<T, F>
    where T: Sized + Send, F: Fn() -> T {
    /// Create a new `ThreadLocalStore`.
    ///
    /// Takes a closure that returns the initial value
    /// of the thread local data.
    pub fn new(init: F) -> Self {
        ThreadLocalStore { inner: CachedThreadLocal::new(), init: init }
    }

    /// Get a `RefCell` containing the thread local data.
    pub fn get(&self) -> &RefCell<T> {
        self.inner.get_or(|| Box::new(RefCell::new((self.init)())))
    }

    /// Shortcut to sum the local values if the local data
    /// can be iterated into `AddAssign` items.
    pub fn add_local_values<I>(self, output: &mut [I])
        where T: IntoIterator<Item=I>, I: AddAssign<I> {
        for local_values in self {
            for (o, v) in output.iter_mut().zip(local_values) {
                *o += v;
            }
        }
    }
}



