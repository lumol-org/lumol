// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

//! Convert TOML values to Lumol types.


use FromToml;
use FromTomlWithData;
use error::{Error, Result};
use extract;

use lumol::energy::{BornMayerHuggins, Buckingham, Gaussian, MorsePotential, Torsion};
use lumol::energy::{CosineHarmonic, Harmonic, LennardJones, NullPotential, Mie};
use lumol::energy::{Ewald, Wolf};
use lumol::energy::{PairPotential, TableComputation};
use lumol::units;
use toml::value::Table;

impl FromToml for NullPotential {
    fn from_toml(_: &Table) -> Result<NullPotential> {
        Ok(NullPotential)
    }
}

impl FromToml for Harmonic {
    fn from_toml(table: &Table) -> Result<Harmonic> {
        let k = extract::str("k", table, "harmonic potential")?;
        let x0 = extract::str("x0", table, "harmonic potential")?;
        Ok(Harmonic {
            k: units::from_str(k)?,
            x0: units::from_str(x0)?,
        })
    }
}

impl FromToml for LennardJones {
    fn from_toml(table: &Table) -> Result<LennardJones> {
        let sigma = extract::str("sigma", table, "Lennard-Jones potential")?;
        let epsilon = extract::str("epsilon", table, "Lennard-Jones potential")?;
        Ok(LennardJones {
            sigma: units::from_str(sigma)?,
            epsilon: units::from_str(epsilon)?,
        })
    }
}

impl FromToml for Mie {
    fn from_toml(table: &Table) -> Result<Mie> {
        let sigma = extract::str("sigma", table, "Mie potential")?;
        let epsilon = extract::str("epsilon", table, "Mie potential")?;
        let m = extract::number("m", table, "Mie potential")?;
        let n = extract::number("n", table, "Mie potential")?;
        Ok(Mie::new(units::from_str(sigma)?, units::from_str(epsilon)?, n as f64, m as f64))
    }
}

impl FromToml for CosineHarmonic {
    fn from_toml(table: &Table) -> Result<CosineHarmonic> {
        let k = extract::str("k", table, "cosine harmonic potential")?;
        let x0 = extract::str("x0", table, "cosine harmonic potential")?;
        Ok(CosineHarmonic::new(units::from_str(k)?, units::from_str(x0)?))
    }
}

impl FromToml for Torsion {
    fn from_toml(table: &Table) -> Result<Torsion> {
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
    fn from_toml(table: &Table) -> Result<Buckingham> {
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
    fn from_toml(table: &Table) -> Result<BornMayerHuggins> {
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

impl FromToml for MorsePotential {
    fn from_toml(table: &Table) -> Result<MorsePotential> {
        let a = extract::str("A", table, "Morse potential")?;
        let depth = extract::str("depth", table, "Morse potential")?;
        let x0 = extract::str("x0", table, "Morse potential")?;
        Ok(MorsePotential {
            a: units::from_str(a)?,
            depth: units::from_str(depth)?,
            x0: units::from_str(x0)?,
        })
    }
}

impl FromToml for Gaussian {
    fn from_toml(table: &Table) -> Result<Gaussian> {
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
    type Data = Box<PairPotential>;

    fn from_toml(table: &Table, potential: Box<PairPotential>) -> Result<TableComputation> {
        let table = try!(table["table"]
            .as_table()
            .ok_or(Error::from("'table' key in computation must be a TOML table")));
        let n = extract::uint("n", table, "table computation")?;
        let max = extract::str("max", table, "table computation")?;
        Ok(TableComputation::new(potential, n as usize, units::from_str(max)?))
    }
}

impl FromToml for Wolf {
    fn from_toml(table: &Table) -> Result<Wolf> {
        let cutoff = extract::str("cutoff", table, "Wolf coulombic potential")?;
        Ok(Wolf::new(units::from_str(cutoff)?))
    }
}

impl FromToml for Ewald {
    fn from_toml(table: &Table) -> Result<Ewald> {
        let cutoff = extract::str("cutoff", table, "Ewald coulombic potential")?;
        let kmax = extract::uint("kmax", table, "Ewald coulombic potential")?;
        Ok(Ewald::new(units::from_str(cutoff)?, kmax as usize))
    }
}
