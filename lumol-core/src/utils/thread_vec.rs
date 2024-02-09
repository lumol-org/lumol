// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

use std::cell::{RefCell, RefMut};
use std::ops::AddAssign;

use thread_local::ThreadLocal;

/// A collection of vectors, one by thread using this struct. All the vectors
/// are wrapped in `RefCell` to ensure mutability of the values. The thread
/// creating the `ThreadLocalVec` will get faster access to the underlying data.
pub struct ThreadLocalVec<T> where T: Send {
    inner: ThreadLocal<RefCell<Vec<T>>>,
    size: usize
}

impl<T: Send + Default + Clone> ThreadLocalVec<T> {
    /// Create a new `ThreadLocalVec` with the given size, initializing the
    /// values with `T::default`.
    pub fn with_size(size: usize) -> Self {
        let inner = ThreadLocal::new();
        // Set the current thread as owner of the data
        let _ = inner.get_or(|| RefCell::new(vec![T::default(); size]));
        ThreadLocalVec {
            inner: inner,
            size: size,
        }
    }

    /// Mutably borrow the thread local vector if it already exists, or create
    /// it and then borrow it.
    pub fn borrow_mut(&self) -> RefMut<'_, Vec<T>> {
        self.inner
            .get_or(|| RefCell::new(vec![T::default(); self.size]))
            .borrow_mut()
    }
}

impl<T: Send> ThreadLocalVec<T> {
    /// Get an iterator over all the vectors created by the different threads
    pub fn into_iter(self) -> impl Iterator<Item = Vec<T>> {
        self.inner
            .into_iter()
            .map(|cell| cell.into_inner())
    }

    /// Sum the values from all the vectors created by the different threads in
    /// the `output` buffer
    pub fn sum_into(self, output: &mut [T]) where T: AddAssign {
        for local in self.into_iter() {
            for (a, b) in zip!(&mut *output, local) {
                *a += b;
            }
        }
    }
}
