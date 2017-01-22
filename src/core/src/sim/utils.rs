/// Helper struct that can wrap an algorithm to make
/// it run only a fraction of the times it is called.
///
/// # Example
///
/// ```
/// use lumol::sim::Alternator;
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
///             "Hello world!"
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

    /// Access the base algorithm.
    pub fn base(&self) -> &T {
        &self.base
    }

    /// Access the base algorithm.
    pub fn base_mut(&mut self) -> &mut T {
        &mut self.base
    }
}
