use lumol::energy::{Potential, PairPotential};

#[derive(Clone, Copy)]
pub struct Mie {
    /// Distance constant
    sigma: f64,
    /// Exponent of repulsive contribution
    n: f64,
    /// Exponent of attractive contribution
    m: f64,
    /// Energetic prefactor computed from the exponents and epsilon
    prefactor: f64,
}

impl Mie {
    pub fn new(sigma: f64, epsilon: f64, n: f64, m: f64) -> Mie {
        if m >= n {
            panic!("The repulsive exponent n has to be larger than the attractive exponent m")
        };
        let prefactor = n / (n - m) * (n / m).powf(m / (n - m)) * epsilon;
        Mie {
            sigma,
            n,
            m,
            prefactor,
        }
    }
}

impl Potential for Mie {
    fn energy(&self, r: f64) -> f64 {
        let sigma_r = self.sigma / r;
        let repulsive = f64::powf(sigma_r, self.n);
        let attractive = f64::powf(sigma_r, self.m);

        self.prefactor * (repulsive - attractive)
    }

    fn force(&self, r: f64) -> f64 {
        let sigma_r = self.sigma / r;
        let repulsive = f64::powf(sigma_r, self.n);
        let attractive = f64::powf(sigma_r, self.m);

        -self.prefactor * (self.n * repulsive - self.m * attractive) / r
    }
}

impl PairPotential for Mie {
    fn tail_energy(&self, cutoff: f64) -> f64 {
        if self.m < 3.0 {
            return 0.0;
        };
        let sigma_rc = self.sigma / cutoff;
        let n_3 = self.n - 3.0;
        let m_3 = self.m - 3.0;
        let repulsive = f64::powf(sigma_rc, n_3);
        let attractive = f64::powf(sigma_rc, m_3);

        -self.prefactor * self.sigma.powi(3) * (repulsive / n_3 - attractive / m_3)
    }

    fn tail_virial(&self, cutoff: f64) -> f64 {
        if self.m < 3.0 {
            return 0.0;
        };
        let sigma_rc = self.sigma / cutoff;
        let n_3 = self.n - 3.0;
        let m_3 = self.m - 3.0;
        let repulsive = f64::powf(sigma_rc, n_3);
        let attractive = f64::powf(sigma_rc, m_3);

        -self.prefactor * self.sigma.powi(3) * (repulsive * self.n / n_3 - attractive * self.m / m_3)
    }
}
