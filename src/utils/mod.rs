/* Cymbalum, Molecular Simulation in Rust - Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
 */
//! Various internal utilities, which do not have there own module

mod bcount;
pub use self::bcount::Bc;

#[cfg(test)]
mod xyz;
#[cfg(test)]
pub use self::xyz::system_from_xyz;

/// A simple macro to implement Clone for Box<Trait>, without requiring that
/// Trait is Clone. This works by creating a new trait (`BoxCloneTrait`) and
/// making the first Trait inherit the `BoxCloneTrait`. `BoxCloneTrait` is
/// automatically implemented for all `T: Clone + Trait`.
///
/// Usage:
///
/// ```rust
/// trait Foo: BoxCloneFoo {}
///
/// impl_box_clone!(Foo, BoxCloneFoo, box_clone_foo);
/// ```
macro_rules! impl_box_clone {
    ($Trait: ident, $BoxClone: ident, $box_clone: ident) => (
        /// This is an internal implementation detail for cloning Box<Trait>,
        /// where Trait is used as a trait object.
        pub trait $BoxClone {
            /// Get a cloned of self as a boxed trait.
            fn $box_clone(&self) -> Box<$Trait>;
        }

        impl<T: Clone + $Trait + 'static> $BoxClone for T {
            fn $box_clone(&self) -> Box<$Trait> {
                Box::new(self.clone())
            }
        }

        impl Clone for Box<$Trait> {
            fn clone(&self) -> Box<$Trait> {
                self.$box_clone()
            }
        }
    );
}
