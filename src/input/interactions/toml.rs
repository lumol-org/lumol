// Cymbalum, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

//! Convert TOML values to Cymbalum types.
use toml::Table;

use super::{Error, Result, FromToml, FromTomlWithPairs};

use potentials::{Harmonic, LennardJones, NullPotential, CosineHarmonic, Torsion};
use potentials::{PairPotential, TableComputation, CutoffComputation};

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
            let k = try!(::units::from_str(k));
            let x0 = try!(::units::from_str(x0));
            Ok(Harmonic{k: k, x0: x0})
        } else {
            Err(
                Error::from("'k' and 'x0' must be strings in harmonic potential")
            )
        }
    }
}

impl FromToml for LennardJones {
    fn from_toml(table: &Table) -> Result<LennardJones> {
        let sigma = try_extract_parameter!(table, "sigma", "Lennard-Jones potential");
        let epsilon = try_extract_parameter!(table, "epsilon", "Lennard-Jones potential");

        if let (Some(sigma), Some(epsilon)) = (sigma.as_str(), epsilon.as_str()) {
            let sigma = try!(::units::from_str(sigma));
            let epsilon = try!(::units::from_str(epsilon));
            Ok(LennardJones{sigma: sigma, epsilon: epsilon})
        } else {
            Err(
                Error::from("'epsilon' and 'sigma' must be strings in Lennard-Jones potential")
            )
        }
    }
}

impl FromToml for CosineHarmonic {
    fn from_toml(table: &Table) -> Result<CosineHarmonic> {
        let k = try_extract_parameter!(table, "k", "cosine harmonic potential");
        let x0 = try_extract_parameter!(table, "x0", "cosine harmonic potential");

        if let (Some(k), Some(x0)) = (k.as_str(), x0.as_str()) {
            let k = try!(::units::from_str(k));
            let x0 = try!(::units::from_str(x0));
            Ok(CosineHarmonic::new(k, x0))
        } else {
            Err(
                Error::from("'k' and 'x0' must be strings in cosine harmonic potential")
            )
        }
    }
}

impl FromToml for Torsion {
    fn from_toml(table: &Table) -> Result<Torsion> {
        let k = try_extract_parameter!(table, "k", "torsion potential");
        let n = try_extract_parameter!(table, "n", "torsion potential");
        let delta = try_extract_parameter!(table, "delta", "torsion potential");

        if let (Some(n), Some(k), Some(delta)) = (n.as_integer(), k.as_str(), delta.as_str()) {
            let k = try!(::units::from_str(k));
            let delta = try!(::units::from_str(delta));
            Ok(Torsion{n: n as usize, k: k, delta: delta})
        } else {
            Err(
                Error::from("'k' and 'delta' must be strings in torsion potential, \
                and 'n' must be an integer")
            )
        }
    }
}

impl FromTomlWithPairs for CutoffComputation {
    fn from_toml(table: &Table, potential: Box<PairPotential>) -> Result<CutoffComputation> {
        let cutoff = try_extract_parameter!(table, "cutoff", "cutoff computation");
        if let Some(cutoff) = cutoff.as_str() {
            let cutoff = try!(::units::from_str(cutoff));
            Ok(CutoffComputation::new(potential, cutoff))
        } else {
            Err(
                Error::from("'cutoff' must be a string in cutoff computation")
            )
        }
    }
}

impl FromTomlWithPairs for TableComputation {
    fn from_toml(table: &Table, potential: Box<PairPotential>) -> Result<TableComputation> {
        let table = table["table"].as_table().unwrap();
        let n = try_extract_parameter!(table, "n", "table computation");
        let max = try_extract_parameter!(table, "max", "table computation");
        if let (Some(n), Some(max)) = (n.as_integer(), max.as_str()) {
            let max = try!(::units::from_str(max));
            Ok(TableComputation::new(potential, n as usize, max))
        } else {
            Err(Error::from(
                "'max' must be a string and 'n' and integere in table computation"
            ))
        }
    }
}

#[cfg(test)]
mod tests {
    use system::System;
    use input::read_interactions;
    use input::interactions::testing::bad_interactions;

    #[test]
    fn bad_potentials() {
        for path in bad_interactions("toml") {
            println!("{:?}", path.display());
            let mut system = System::new();
            assert!(read_interactions(&mut system, path).is_err());
        }
    }
}
