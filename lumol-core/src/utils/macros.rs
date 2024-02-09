// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

macro_rules! zip {
    (@map $pattern:pat => $tuple:expr) => {
        |$pattern| $tuple
    };
    (@map $pattern:pat => ( $($tuple:tt)* ) , $_iter:expr $(, $tail:expr )*) => {
        zip!(@map ($pattern, b) => ( $($tuple)*, b ) $( , $tail )*)
    };
    ($first:expr $( , $rest:expr )* $(,)*) => {
        ::std::iter::IntoIterator::into_iter($first)
            $(.zip($rest))*
            .map(
                zip!(@map a => (a) $( , $rest )*)
            )
    };
}

/// A simple macro to implement Clone for Box<Trait>, without requiring that
/// Trait is Clone. This works by creating a new trait (`BoxCloneTrait`) and
/// making the first Trait inherit the `BoxCloneTrait`. `BoxCloneTrait` is
/// automatically implemented for all `T: Clone + Trait`.
///
/// Usage:
///
/// ```ignore
/// trait Foo: BoxCloneFoo {}
///
/// impl_box_clone!(Foo, BoxCloneFoo, box_clone_foo);
/// ```
macro_rules! impl_box_clone {
    ($Trait: ident, $BoxClone: ident, $box_clone: ident) => (
        #[doc(hidden)]
        /// This is an internal implementation detail for cloning `Box<Trait>`
        pub trait $BoxClone {
            /// Get a cloned of self as a boxed trait.
            fn $box_clone(&self) -> Box<dyn $Trait>;
        }

        impl<T: Clone + $Trait + 'static> $BoxClone for T {
            fn $box_clone(&self) -> Box<dyn $Trait> {
                Box::new(self.clone())
            }
        }

        impl Clone for Box<dyn $Trait> {
            fn clone(&self) -> Box<dyn $Trait> {
                self.$box_clone()
            }
        }
    );
}

/// Implement $Lhs -- $Rhs arithmetic operations for all variation of by
/// value, by reference and by mutable reference of $Rhs and $Lhs.
macro_rules! impl_arithmetic {
    ($Lhs:ty, $Rhs:ty, $Op:ident, $op:ident, $Output:ty, $sel:ident, $other:ident, $res:expr) => (
        impl $Op<$Rhs> for $Lhs {
            type Output = $Output;
            #[inline] fn $op($sel, $other: $Rhs) -> $Output {
                $res
            }
        }

        impl<'a> $Op<$Rhs> for &'a $Lhs {
            type Output = $Output;
            #[inline] fn $op($sel, $other: $Rhs) -> $Output {
                $res
            }
        }

        impl<'a> $Op<&'a $Rhs> for $Lhs {
            type Output = $Output;
            #[inline] fn $op($sel, $other: &'a $Rhs) -> $Output {
                $res
            }
        }

        impl<'a, 'b> $Op<&'a $Rhs> for &'b $Lhs {
            type Output = $Output;
            #[inline] fn $op($sel, $other: &'a $Rhs) -> $Output {
                $res
            }
        }

        impl<'a, 'b> $Op<&'a mut $Rhs> for &'b mut $Lhs {
            type Output = $Output;
            #[inline] fn $op($sel, $other: &'a mut $Rhs) -> $Output {
                $res
            }
        }

        impl<'a, 'b> $Op<&'a mut $Rhs> for &'b $Lhs {
            type Output = $Output;
            #[inline] fn $op($sel, $other: &'a mut $Rhs) -> $Output {
                $res
            }
        }

        impl<'a, 'b> $Op<&'a $Rhs> for &'b mut $Lhs {
            type Output = $Output;
            #[inline] fn $op($sel, $other: &'a $Rhs) -> $Output {
                $res
            }
        }

        impl<'a> $Op<&'a mut $Rhs> for $Lhs {
            type Output = $Output;
            #[inline] fn $op($sel, $other: &'a mut $Rhs) -> $Output {
                $res
            }
        }

        impl<'a> $Op<$Rhs> for &'a mut $Lhs {
            type Output = $Output;
            #[inline] fn $op($sel, $other: $Rhs) -> $Output {
                $res
            }
        }
    );
}

/// Implement operators `@=` for all variations of references for the right-hand
/// side.
macro_rules! impl_in_place_arithmetic {
    ($Lhs:ty, $Rhs:ty, $Op:ident, $op:ident, $sel:ident, $other:ident, $res:expr) => (
        impl<'a> $Op<$Rhs> for $Lhs {
            #[inline] fn $op(&mut $sel, $other: $Rhs) {
                $res
            }
        }

        impl<'a> $Op<&'a $Rhs> for $Lhs {
            #[inline] fn $op(&mut $sel, $other: &'a $Rhs) {
                $res
            }
        }

        impl<'a> $Op<&'a mut $Rhs> for $Lhs {
            #[inline] fn $op(&mut $sel, $other: &'a mut $Rhs) {
                $res
            }
        }
    )
}

/// Implement $Lhs -- scalar arithmetic operations for all variation of by
/// value, by reference and by mutable reference $Lhs.
macro_rules! lsh_scalar_arithmetic {
    ($Lhs: ty, $Op:ident, $op:ident, $Output:ty, $sel:ident, $other:ident, $res:expr) => (
        impl $Op<f64> for $Lhs {
            type Output = $Output;
            #[inline] fn $op($sel, $other: f64) -> $Output {
                $res
            }
        }

        impl<'a> $Op<f64> for &'a $Lhs {
            type Output = $Output;
            #[inline] fn $op($sel, $other: f64) -> $Output {
                $res
            }
        }

        impl<'a> $Op<f64> for &'a mut $Lhs {
            type Output = $Output;
            #[inline] fn $op($sel, $other: f64) -> $Output {
                $res
            }
        }
    );
}

/// Implement scalar -- $Rhs arithmetic operations for all variation of by
/// value, by reference and by mutable reference of $Rhs.
macro_rules! rhs_scalar_arithmetic {
    ($Rhs:ty, $Op:ident, $op:ident, $Output:ty, $sel:ident, $other:ident, $res:expr) => (
        impl $Op<$Rhs> for f64 {
            type Output = $Output;
            #[inline] fn $op($sel, $other: $Rhs) -> $Output {
                $res
            }
        }

        impl<'a> $Op<&'a $Rhs> for f64 {
            type Output = $Output;
            #[inline] fn $op($sel, $other: &'a $Rhs) -> $Output {
                $res
            }
        }

        impl<'a> $Op<&'a mut $Rhs> for f64 {
            type Output = $Output;
            #[inline] fn $op($sel, $other: &'a mut $Rhs) -> $Output {
                $res
            }
        }
    );
}
