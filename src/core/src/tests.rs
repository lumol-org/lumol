// Lumol, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

/// Macro for testing approximated equality of floats.
macro_rules! assert_approx_eq {
    ($left: expr, $right: expr, $tol: expr) => (
        {
            match (&($left), &($right), ($tol)) {
                (left_val , right_val, tol_val) => {
                    let delta: f64 = (*left_val - *right_val) as f64;
                    let delta = delta.abs();
                    if !(delta < tol_val) {
                        panic!(
                            "assertion failed: `(left ≈ right)` \
                             (left: `{:?}`, right: `{:?}`) \
                             with ∆={:1.1e} (allowed ∆={:e})\
                             "
                        , *left_val , *right_val, delta, tol_val)
                    }
                }
            }
        }
    );
    ($left: expr, $right: expr) => (assert_approx_eq!(($left), ($right), 1e-15))
}

#[cfg(test)]
mod tests {
    //use std::num::Float;

    #[test]
    fn approx_eq() {
        assert_approx_eq!(1.0, 1.0, 1e-90);
        assert_approx_eq!(1.1, 1.1000003, 1e-3);
    }

    #[test]
    fn default_vaue() {
        assert_approx_eq!(1.0, 1.0);
        assert_approx_eq!(1.1, 1.1 + 1e-16);
    }

    #[test]
    #[should_panic]
    fn not_eq_default_vaue() {
        assert_approx_eq!(1.1, 1.1 + 1e-15);
    }

    #[test]
    #[should_panic]
    fn not_approx_eq() {
        assert_approx_eq!(1.1, 1.1000003, 1e-8);
    }

    #[test]
    fn all_float_types() {
        assert_approx_eq!(1.0f32, 1.0f32, 1e-90);
        assert_approx_eq!(1.0f64, 1.0f64, 1e-90);
    }
}
