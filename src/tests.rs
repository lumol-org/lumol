/*
 * Cymbalum, Molecular Simulation in Rust
 * Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
*/

macro_rules! assert_approx_eq {
    ($left: expr, $right: expr, $tol: expr) => (
        {
            match (&($left), &($right), ($tol)) {
                (left_val , right_val, tol_val) => {
                    if !(f64::abs(*left_val - *right_val) < tol_val) {
                        panic!(
                            "assertion failed: `(left == right)` \
                             (left: `{:?}`, right: `{:?}`)"
                        , *left_val , *right_val)
                    }
                }
            }
        }
    )
}
