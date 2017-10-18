// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

//! Convert TOML values to Lumol types.
use toml::value::Table;

use error::{Error, Result};
use FromToml;
use FromTomlWithData;

use lumol::energy::{Harmonic, LennardJones, NullPotential, CosineHarmonic};
use lumol::energy::{Torsion, Buckingham, BornMayerHuggins, MorsePotential, Gaussian};
use lumol::energy::{Wolf, Ewald};
use lumol::energy::{PairPotential, TableComputation};

macro_rules! try_extract_parameter {
    ($table: expr, $key: expr, $context: expr) => (
        match $table.get($key) {
            Some(value) => value,
            None => return Err(Error::from(
                String::from("Missing '") + $key + "' in " + $context
            ))
        }
    )
}

impl FromToml for NullPotential {
    fn from_toml(_: &Table) -> Result<NullPotential> {
        Ok(NullPotential)
    }
}

impl FromToml for Harmonic {
    fn from_toml(table: &Table) -> Result<Harmonic> {
        let k = try_extract_parameter!(table, "k", "harmonic potential");
        let x0 = try_extract_parameter!(table, "x0", "harmonic potential");

        if let (Some(k), Some(x0)) = (k.as_str(), x0.as_str()) {
            let k = try!(::lumol::units::from_str(k));
            let x0 = try!(::lumol::units::from_str(x0));
            Ok(Harmonic { k: k, x0: x0 })
        } else {
            Err(Error::from("'k' and 'x0' must be strings in harmonic potential"))
        }
    }
}

impl FromToml for LennardJones {
    fn from_toml(table: &Table) -> Result<LennardJones> {
        let sigma = try_extract_parameter!(table, "sigma", "Lennard-Jones potential");
        let epsilon = try_extract_parameter!(table, "epsilon", "Lennard-Jones potential");

        if let (Some(sigma), Some(epsilon)) = (sigma.as_str(), epsilon.as_str()) {
            let sigma = try!(::lumol::units::from_str(sigma));
            let epsilon = try!(::lumol::units::from_str(epsilon));
            Ok(LennardJones {
                sigma: sigma,
                epsilon: epsilon,
            })
        } else {
            Err(Error::from("'epsilon' and 'sigma' must be strings in Lennard-Jones potential"))
        }
    }
}

impl FromToml for CosineHarmonic {
    fn from_toml(table: &Table) -> Result<CosineHarmonic> {
        let k = try_extract_parameter!(table, "k", "cosine harmonic potential");
        let x0 = try_extract_parameter!(table, "x0", "cosine harmonic potential");

        if let (Some(k), Some(x0)) = (k.as_str(), x0.as_str()) {
            let k = try!(::lumol::units::from_str(k));
            let x0 = try!(::lumol::units::from_str(x0));
            Ok(CosineHarmonic::new(k, x0))
        } else {
            Err(Error::from("'k' and 'x0' must be strings in cosine harmonic potential"))
        }
    }
}

impl FromToml for Torsion {
    fn from_toml(table: &Table) -> Result<Torsion> {
        let k = try_extract_parameter!(table, "k", "torsion potential");
        let n = try_extract_parameter!(table, "n", "torsion potential");
        let delta = try_extract_parameter!(table, "delta", "torsion potential");

        if let (Some(n), Some(k), Some(delta)) = (n.as_integer(), k.as_str(), delta.as_str()) {
            let k = try!(::lumol::units::from_str(k));
            let delta = try!(::lumol::units::from_str(delta));
            Ok(Torsion {
                n: n as usize,
                k: k,
                delta: delta,
            })
        } else {
            Err(Error::from("'k' and 'delta' must be strings in torsion potential, and 'n' must \
                             be an integer"))
        }
    }
}

impl FromToml for Buckingham {
    fn from_toml(table: &Table) -> Result<Buckingham> {
        let a = try_extract_parameter!(table, "A", "Buckingham potential");
        let c = try_extract_parameter!(table, "C", "Buckingham potential");
        let rho = try_extract_parameter!(table, "rho", "Buckingham potential");

        if let (Some(a), Some(c), Some(rho)) = (a.as_str(), c.as_str(), rho.as_str()) {
            let a = try!(::lumol::units::from_str(a));
            let c = try!(::lumol::units::from_str(c));
            let rho = try!(::lumol::units::from_str(rho));
            Ok(Buckingham {
                a: a,
                c: c,
                rho: rho,
            })
        } else {
            Err(Error::from("'A', 'C' and 'rho' must be strings in Buckingham potential"))
        }
    }
}

impl FromToml for BornMayerHuggins {
    fn from_toml(table: &Table) -> Result<BornMayerHuggins> {
        let a = try_extract_parameter!(table, "A", "Born-Mayer-Huggins potential");
        let c = try_extract_parameter!(table, "C", "Born-Mayer-Huggins potential");
        let d = try_extract_parameter!(table, "D", "Born-Mayer-Huggins potential");
        let sigma = try_extract_parameter!(table, "sigma", "Born-Mayer-Huggins potential");
        let rho = try_extract_parameter!(table, "rho", "Born-Mayer-Huggins potential");

        if let (Some(a), Some(c), Some(d), Some(sigma), Some(rho)) =
            (a.as_str(), c.as_str(), d.as_str(), sigma.as_str(), rho.as_str()) {
            let a = try!(::lumol::units::from_str(a));
            let c = try!(::lumol::units::from_str(c));
            let d = try!(::lumol::units::from_str(d));
            let sigma = try!(::lumol::units::from_str(sigma));
            let rho = try!(::lumol::units::from_str(rho));
            Ok(BornMayerHuggins {
                a: a,
                c: c,
                d: d,
                sigma: sigma,
                rho: rho,
            })
        } else {
            Err(Error::from("'A', 'C', 'D', 'sigma' and 'rho' must be strings in \
                             Born-Mayer-Huggins potential"))
        }
    }
}

impl FromToml for MorsePotential {
    fn from_toml(table: &Table) -> Result<MorsePotential> {
        let a = try_extract_parameter!(table, "a", "Morse potential");
        let depth = try_extract_parameter!(table, "depth", "Morse potential");
        let x0 = try_extract_parameter!(table, "x0", "Morse potential");

        if let (Some(a), Some(depth), Some(x0)) = (a.as_str(), depth.as_str(), x0.as_str()) {
            let a = try!(::lumol::units::from_str(a));
            let depth = try!(::lumol::units::from_str(depth));
            let x0 = try!(::lumol::units::from_str(x0));
            Ok(MorsePotential {
                a: a,
                depth: depth,
                x0: x0,
            })
        } else {
            Err(Error::from("'depth', 'a' and 'x0' must be strings in Morse potential"))
        }
    }
}

impl FromToml for Gaussian {
    fn from_toml(table: &Table) -> Result<Gaussian> {
        let a = try_extract_parameter!(table, "a", "Gaussian");
        let b = try_extract_parameter!(table, "b", "Gaussian");

        if let (Some(a), Some(b)) = (a.as_str(), b.as_str()) {
            let a = try!(::lumol::units::from_str(a));
            let b = try!(::lumol::units::from_str(b));
            if b <= 0.0 {
                Err(Error::from("'b' has to be positive"))
            } else {
                Ok(Gaussian::new(a, b))
            }
        } else {
            Err(Error::from("'a' and 'b' must be strings in Gaussian"))
        }
    }
}
/// ***************************************************************************

impl FromTomlWithData for TableComputation {
    type Data = Box<PairPotential>;

    fn from_toml(table: &Table, potential: Box<PairPotential>) -> Result<TableComputation> {
        let table = try!(table["table"]
            .as_table()
            .ok_or(Error::from("'table' key in computation must be a TOML table")));
        let n = try_extract_parameter!(table, "n", "table computation");
        let max = try_extract_parameter!(table, "max", "table computation");
        if let (Some(n), Some(max)) = (n.as_integer(), max.as_str()) {
            let max = try!(::lumol::units::from_str(max));
            Ok(TableComputation::new(potential, n as usize, max))
        } else {
            Err(Error::from("'max' must be a string and 'n' and integere in table computation"))
        }
    }
}

/// ***************************************************************************

impl FromToml for Wolf {
    fn from_toml(table: &Table) -> Result<Wolf> {
        let cutoff = try_extract_parameter!(table, "cutoff", "wolf potential");
        if let Some(cutoff) = cutoff.as_str() {
            let cutoff = try!(::lumol::units::from_str(cutoff));
            Ok(Wolf::new(cutoff))
        } else {
            Err(Error::from("'cutoff' parameter must be a string in Wolf potential"))
        }
    }
}

impl FromToml for Ewald {
    fn from_toml(table: &Table) -> Result<Ewald> {
        let cutoff = try_extract_parameter!(table, "cutoff", "ewald potential");
        let kmax = try_extract_parameter!(table, "kmax", "ewald potential");

        if let (Some(cutoff), Some(kmax)) = (cutoff.as_str(), kmax.as_integer()) {
            let cutoff = try!(::lumol::units::from_str(cutoff));
            if kmax < 0 {
                Err(Error::from("'kmax' can not be negative in Ewald potential"))
            } else {
                Ok(Ewald::new(cutoff, kmax as usize))
            }
        } else {
            Err(Error::from("'cutoff' must be a string and 'kmax' an integer in Ewald potential"))
        }
    }
}
