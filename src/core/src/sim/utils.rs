// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

//! Module for small utility structs

/// Helper struct that can wrap an algorithm to make
/// it run only a fraction of the times it is called.
///
/// # Example
///
/// ```
/// use lumol_core::sim::Alternator;
///
/// struct HelloWorld;
///
/// trait SaySomething {
///     fn say_something(&mut self) -> &'static str;
/// }
///
/// impl SaySomething for HelloWorld {
///     fn say_something(&mut self) -> &'static str {
///         "Hello world!"
///     }
/// }
///
/// impl<T> SaySomething for Alternator<T> where T: SaySomething {
///     fn say_something(&mut self) -> &'static str {
///         if self.can_run() {
///             self.as_mut().say_something()
///         } else {
///             ""
///         }
///     }
/// }
///
/// let mut alternator = Alternator::new(3, HelloWorld);
///
/// assert_eq!(alternator.say_something(), "");
/// assert_eq!(alternator.say_something(), "");
/// assert_eq!(alternator.say_something(), "Hello world!");
/// assert_eq!(alternator.say_something(), "");
/// assert_eq!(alternator.say_something(), "");
/// assert_eq!(alternator.say_something(), "Hello world!");
///
/// ```
#[derive(Debug)]
pub struct Alternator<T> {
    every: u64,
    count: u64,
    base: T
}


impl<T> Alternator<T> {
    /// Wrap the algorithm `base` to call it only
    /// every `every` time.
    pub fn new(every: u64, base: T) -> Alternator<T> {
        Alternator { every: every, base: base, count: 0 }
    }

    /// Check if this is the appropriate time to run the algorithm.
    pub fn can_run(&mut self) -> bool {
        self.count += 1;
        self.count % self.every == 0
    }
}

impl<T> AsRef<T> for Alternator<T> {
    /// Access the base algorithm as a reference.
    fn as_ref(&self) -> &T {
        &self.base
    }
}

impl<T> AsMut<T> for Alternator<T> {
    /// Access the base algorithm as a mutable reference.
    fn as_mut(&mut self) -> &mut T {
        &mut self.base
    }
}
