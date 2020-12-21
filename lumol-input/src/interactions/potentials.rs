// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

#![allow(clippy::wildcard_imports)]

use toml::value::Table;
use log::warn;

use lumol_core::units;
use lumol_core::energy::*;
use lumol_core::Configuration;

use crate::{Error, FromToml, FromTomlWithData, FromTomlWithRefData};
use crate::extract;

impl FromToml for NullPotential {
    fn from_toml(_: &Table) -> Result<NullPotential, Error> {
        Ok(NullPotential)
    }
}

impl FromToml for Harmonic {
    fn from_toml(table: &Table) -> Result<Harmonic, Error> {
        let k = extract::str("k", table, "harmonic potential")?;
        let x0 = extract::str("x0", table, "harmonic potential")?;
        Ok(Harmonic {
            k: units::from_str(k)?,
            x0: units::from_str(x0)?,
        })
    }
}

impl FromToml for LennardJones {
    fn from_toml(table: &Table) -> Result<LennardJones, Error> {
        let sigma = extract::str("sigma", table, "Lennard-Jones potential")?;
        let epsilon = extract::str("epsilon", table, "Lennard-Jones potential")?;
        Ok(LennardJones {
            sigma: units::from_str(sigma)?,
            epsilon: units::from_str(epsilon)?,
        })
    }
}

impl FromToml for Mie {
    fn from_toml(table: &Table) -> Result<Mie, Error> {
        let sigma = extract::str("sigma", table, "Mie potential")?;
        let epsilon = extract::str("epsilon", table, "Mie potential")?;
        let m = extract::number("m", table, "Mie potential")?;
        let n = extract::number("n", table, "Mie potential")?;

        if m < 3.0 {
            warn!("'m' is smaller than 3. Tail corrections for Mie potential are set to zero.");
        };

        Ok(Mie::new(
                units::from_str(sigma)?,
                units::from_str(epsilon)?,
                n as f64,
                m as f64)
        )
    }
}

impl FromToml for CosineHarmonic {
    fn from_toml(table: &Table) -> Result<CosineHarmonic, Error> {
        let k = extract::str("k", table, "cosine harmonic potential")?;
        let x0 = extract::str("x0", table, "cosine harmonic potential")?;
        Ok(CosineHarmonic::new(units::from_str(k)?, units::from_str(x0)?))
    }
}

impl FromToml for Torsion {
    fn from_toml(table: &Table) -> Result<Torsion, Error> {
        let n = extract::uint("n", table, "torsion potential")?;
        let k = extract::str("k", table, "torsion potential")?;
        let delta = extract::str("delta", table, "torsion potential")?;
        Ok(Torsion {
            n: n as usize,
            k: units::from_str(k)?,
            delta: units::from_str(delta)?,
        })
    }
}

impl FromToml for Buckingham {
    fn from_toml(table: &Table) -> Result<Buckingham, Error> {
        let a = extract::str("A", table, "Buckingham potential")?;
        let c = extract::str("C", table, "Buckingham potential")?;
        let rho = extract::str("rho", table, "Buckingham potential")?;

        Ok(Buckingham {
            a: units::from_str(a)?,
            c: units::from_str(c)?,
            rho: units::from_str(rho)?,
        })
    }
}

impl FromToml for BornMayerHuggins {
    fn from_toml(table: &Table) -> Result<BornMayerHuggins, Error> {
        let a = extract::str("A", table, "Born-Mayer-Huggins potential")?;
        let c = extract::str("C", table, "Born-Mayer-Huggins potential")?;
        let d = extract::str("D", table, "Born-Mayer-Huggins potential")?;
        let rho = extract::str("rho", table, "Born-Mayer-Huggins potential")?;
        let sigma = extract::str("sigma", table, "Born-Mayer-Huggins potential")?;

        Ok(BornMayerHuggins {
            a: units::from_str(a)?,
            c: units::from_str(c)?,
            d: units::from_str(d)?,
            sigma: units::from_str(sigma)?,
            rho: units::from_str(rho)?,
        })
    }
}

impl FromToml for Morse {
    fn from_toml(table: &Table) -> Result<Morse, Error> {
        let a = extract::str("A", table, "Morse potential")?;
        let depth = extract::str("depth", table, "Morse potential")?;
        let x0 = extract::str("x0", table, "Morse potential")?;
        Ok(Morse {
            a: units::from_str(a)?,
            depth: units::from_str(depth)?,
            x0: units::from_str(x0)?,
        })
    }
}

impl FromToml for Gaussian {
    fn from_toml(table: &Table) -> Result<Gaussian, Error> {
        let a = units::from_str(extract::str("A", table, "Gaussian potential")?)?;
        let b = units::from_str(extract::str("B", table, "Gaussian potential")?)?;

        if b <= 0.0 {
            Err(Error::from("'B' parameter has to be positive in Gaussian potential"))
        } else {
            Ok(Gaussian::new(a, b))
        }
    }
}

impl FromTomlWithData for TableComputation {
    type Data = Box<dyn PairPotential>;

    fn from_toml(table: &Table, potential: Box<dyn PairPotential>) -> Result<TableComputation, Error> {
        let table = table["table"].as_table().ok_or(
            Error::from("'table' key in computation must be a TOML table")
        )?;

        let n = extract::uint("n", table, "table computation")?;
        let max = extract::str("max", table, "table computation")?;
        Ok(TableComputation::new(potential, n as usize, units::from_str(max)?))
    }
}

impl FromToml for Wolf {
    fn from_toml(table: &Table) -> Result<Wolf, Error> {
        let cutoff = extract::str("cutoff", table, "Wolf coulombic potential")?;
        Ok(Wolf::new(units::from_str(cutoff)?))
    }
}

impl FromTomlWithRefData for Ewald {
    type Data = Configuration;

    fn from_toml(table: &Table, configuration: &Configuration) -> Result<Ewald, Error> {
        let cutoff = extract::str("cutoff", table, "Ewald coulombic potential")?;
        let cutoff = units::from_str(cutoff)?;

        // Check first for the accuracy key
        if table.contains_key("accuracy") {
            if table.contains_key("kmax") || table.contains_key("alpha") {
                return Err(Error::from(
                    "can not have both accuracy and kmax/alpha in Ewald coulombic potential"
                ));
            }
            let accuracy = extract::number("accuracy", table, "Ewald coulombic potential")?;
            return Ok(Ewald::with_accuracy(cutoff, accuracy, configuration));
        }

        // Else use directly specified parameters
        let kmax = extract::uint("kmax", table, "Ewald coulombic potential")?;
        let alpha = if table.contains_key("alpha") {
            let alpha = extract::str("alpha", table, "Ewald coulombic potential")?;
            Some(units::from_str(alpha)?)
        } else {
            None
        };
        Ok(Ewald::new(cutoff, kmax as usize, alpha))
    }
}
