/* Cymbalum, Molecular Simulation in Rust - Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
 */
extern crate yaml_rust as yaml;
use self::yaml::{ScanError, Yaml, YamlLoader};

use std::io::prelude::*;
use std::io;
use std::result;
use std::fs::File;
use std::path::Path;

use universe::Universe;
use units::UnitParsingError;
use potentials::*;

#[derive(Debug)]
/// Possible causes of error when reading potential files
pub enum Error {
    /// Error in the YAML input file
    YAMLError(ScanError),
    /// IO error
    FileError(io::Error),
    /// File content error: missing sections, bad data types
    ConfigError{
        /// Error message
        msg: String,
    },
    /// Unit parsing error
    UnitError(UnitParsingError),
}

impl From<ScanError> for Error {
    fn from(err: ScanError) -> Error {Error::YAMLError(err)}
}

impl From<io::Error> for Error {
    fn from(err: io::Error) -> Error {Error::FileError(err)}
}

impl<'a> From<&'a str> for Error {
    fn from(err: &'a str) -> Error {
        Error::ConfigError{msg: String::from(err)}
    }
}

impl From<String> for Error {
    fn from(err: String) -> Error {
        Error::ConfigError{msg: err}
    }
}

impl From<UnitParsingError> for Error {
    fn from(err: UnitParsingError) -> Error {Error::UnitError(err)}
}

/// Custom Result for input files
pub type Result<T> = result::Result<T, Error>;

/******************************************************************************/
trait FromYaml where Self: Sized {
    fn from_yaml(node: &Yaml) -> Result<Self>;
}

/// Read a interactions from a Yaml file at `path`, and add these interactions
/// to the `universe`.
///
/// The following format is accepted:
///
/// ```YAML
/// pairs:  # Non bonded atoms pairs
///   - atoms: [He, He]
///     type: LennardJones
///     sigma: 3.4 A
///     epsilon: 0.45 kJ/mol
///     # computations specification is optional
///     computation:
///       type: table
///       numpoints: 5000
///       max: 20.0 A
///   - atoms: [He, Ar]
///     type: LennardJones
///     sigma: 3.8 A
///     epsilon: 0.67 kJ/mol
///     # restriction section is also optional
///     restriction:
///       type: IntraMolecular
///   - atoms: [He, Ar]
///     type: NullPotential
/// bond: # Bonded atoms pairs
///   - atoms: [C, C]
///     type: Harmonic
///     k: 67 kJ/mol/A^2
///     x0: 1.20 A
/// angles: # Molecular angles
///   - atoms: [O, C, C]
///     type: CosineHarmonic
///     k: 67 kJ/mol/deg^2
///     x0: 120 deg
///   - atoms: [C, C, C]
///     type: harmonic
///     k: 300 kJ/mol/deg^2
///     x0: 120 deg
/// dihedrals: # Dihedral angles
///   - atoms: [O, C, C, O]
///     type: harmonic
///     k: 42 kJ/mol/deg^2
///     x0: 180 deg
///   - atoms: [C, C, C, C]
///     type: torsion
///     k: 40 kJ/mol
///     delta: 120 deg
///     n: 4
/// coulomb:
///   - type: wolf
///     cutoff: 10 A
///     charges:
///         O: -1.8
///         Na: 0.9
///     # restriction section is optional here too
///     restriction:
///       type: Scale14
///       scaling: 0.8
/// ```
///
/// The main items are the `"pairs"`, `"bonds"`, `"angles"` and `"dihedrals"`;
/// which contains the data for pair potentials, bonds, angles and dihedral
/// angles potential. This data is orgnised as a vector of hash maps. Each maps
/// has at least a vector of `"atoms"` defining which particles will get this
/// interaction applied to, and a `"type"` parameter defining the type of the
/// potential. Others keys may be supplied depending on the potential type.
///
/// The `"coulomb"` section specify how to compute the coulombic interactions.
/// It should contains the `type` and `charge` key, with data about charges to
/// assign to the system, and which method to use to compute these interactions.
/// Additional keys can exist depending on the actual coulombic solver used.
pub fn read_interactions<P: AsRef<Path>>(universe: &mut Universe, path: P) -> Result<()> {
    let mut file = try!(File::open(path));
    let mut buffer = String::new();
    try!(file.read_to_string(&mut buffer));
    return read_interactions_string(universe, &buffer);
}

/// This is the same as `read_interactions`, but directly read a YAML formated
/// string.
pub fn read_interactions_string(universe: &mut Universe, string: &str) -> Result<()> {
    let doc = &try!(YamlLoader::load_from_str(string))[0];
    if let Some(config) = doc["pairs"].as_vec() {
        try!(read_pairs(universe, config, true));
    }

    if let Some(config) = doc["bonds"].as_vec() {
        try!(read_pairs(universe, config, false));
    }

    if let Some(config) = doc["angles"].as_vec() {
        try!(read_angles(universe, config));
    }

    if let Some(config) = doc["dihedrals"].as_vec() {
        try!(read_dihedrals(universe, config));
    }

    if doc["coulomb"].as_hash().is_some() {
        try!(read_coulomb(universe, &doc["coulomb"]));
    }

    Ok(())
}

/// Read the "pairs" or the "bonds" section in the file. If `pair_potentials`
/// is `true`, then the interactions are added to the pair interactions. Else,
/// the interaction are added to the bond interactions.
fn read_pairs(universe: &mut Universe, pairs: &[Yaml], pair_potentials: bool) -> Result<()> {
    for node in pairs {
        let (a, b) = match node["atoms"].as_vec() {
            Some(vec) => {
                if vec.len() != 2 {
                    return Err(Error::from(
                        format!("Wrong size for 'atoms' section in pair potentials. Should be 2, is {}", vec.len())
                    ));
                }
                if let (Some(a), Some(b)) = (vec[0].as_str(), vec[1].as_str()) {
                    (a, b)
                } else {
                    return Err(Error::from("The atoms names are not strings in pair potential"));
                }
            },
            None => {
                return Err(Error::from("Missing 'atoms' section in pair potential"));
            }
        };

        let restriction = if node["restriction"].as_hash().is_some() {
            Some(try!(read_restriction(&node["restriction"])))
        } else {
            None
        };

        let potential = try!(read_pair_potential(node));
        let potential = if node["computation"].as_hash().is_some() {
            try!(read_pair_computation(&node["computation"], potential))
        } else {
            potential
        };

        if pair_potentials {
            if let Some(restriction) = restriction {
                universe.add_pair_interaction_with_restriction(a, b, potential, restriction);
            } else {
                universe.add_pair_interaction(a, b, potential);
            }
        } else {
            universe.add_bond_interaction(a, b, potential);
        }
    }
    Ok(())
}

fn read_pair_potential(node: &Yaml) -> Result<Box<PairPotential>> {
    match node["type"].as_str() {
        Some(val) => {
            let val: &str = &val.to_lowercase();
            match val {
                "harmonic" => Ok(Box::new(try!(Harmonic::from_yaml(node)))),
                "lennard-jones" | "lennardjones" => Ok(Box::new(try!(LennardJones::from_yaml(node)))),
                "null" | "nullpotential" => Ok(Box::new(try!(NullPotential::from_yaml(node)))),
                val => Err(
                    Error::from(format!("Unknown potential type '{}'", val))
                ),
            }
        },
        None => {
            Err(Error::from("Missing 'type' parameter in pair potential"))
        }
    }
}
/******************************************************************************/
fn read_restriction(node: &Yaml) -> Result<PairRestriction> {
    match node["type"].as_str() {
        Some(val) => {
            let val: &str = &val.to_lowercase();
            match val {
                "none" => Ok(PairRestriction::None),
                "intramolecular" => Ok(PairRestriction::IntraMolecular),
                "intermolecular" => Ok(PairRestriction::InterMolecular),
                "exclude12" => Ok(PairRestriction::Exclude12),
                "exclude13" => Ok(PairRestriction::Exclude13),
                "exclude14" => Ok(PairRestriction::Exclude14),
                "scale14" => {
                    if let Some(scaling) = node["scaling"].as_f64() {
                        if 0.0 <= scaling && scaling <= 1.0 {
                            Ok(PairRestriction::Scale14{scaling: scaling})
                        } else {
                            Err(Error::from(format!("Scaling parameter for Scale14 restriction must be between 0 and 1")))
                        }
                    } else {
                        Err(Error::from(format!("Missing 'scaling' parameter in Scale14 restriction")))
                    }
                },
                val => Err(
                    Error::from(format!("Unknown potential type '{}'", val))
                ),
            }
        },
        None => {
            Err(Error::from("Missing 'type' parameter in restriction section"))
        }
    }
}

/******************************************************************************/
fn read_angles(universe: &mut Universe, angles: &[Yaml]) -> Result<()> {
    for potential in angles {
        let (a, b, c) = match potential["atoms"].as_vec() {
            Some(vec) => {
                if vec.len() != 3 {
                    return Err(Error::from(
                        format!("Wrong size for 'atoms' section in angle potential. Should be 3, is {}", vec.len())
                    ));
                }
                if let (Some(a), Some(b), Some(c)) = (vec[0].as_str(), vec[1].as_str(), vec[2].as_str()) {
                    (a, b, c)
                } else {
                    return Err(Error::from("The atoms names are not strings in angle potential"));
                }
            },
            None => {
                return Err(Error::from("Missing 'atoms' section in angles potential"));
            }
        };

        let potential = try!(read_angle_potential(potential));
        universe.add_angle_interaction(a, b, c, potential);
    }
    Ok(())
}

fn read_angle_potential(node: &Yaml) -> Result<Box<AnglePotential>> {
    match node["type"].as_str() {
        Some(val) => {
            let val: &str = &val.to_lowercase();
            match val {
                "harmonic" => Ok(Box::new(try!(Harmonic::from_yaml(node)))),
                "cosine-harmonic" | "cosineharmonic" => Ok(Box::new(try!(CosineHarmonic::from_yaml(node)))),
                "null" => Ok(Box::new(try!(NullPotential::from_yaml(node)))),
                val => Err(
                    Error::from(format!("Unknown potential type '{}'", val))
                ),
            }
        },
        None => {
            Err(Error::from(format!("Missing 'type' parameter in angle potential")))
        }
    }
}

/******************************************************************************/
fn read_dihedrals(universe: &mut Universe, dihedrals: &[Yaml]) -> Result<()> {
    for potential in dihedrals {
        let (a, b, c, d) = match potential["atoms"].as_vec() {
            Some(vec) => {
                if vec.len() != 4 {
                    return Err(Error::from(
                        format!("Wrong size for 'atoms' section in dihedral potential. Should be 4, is {}", vec.len())
                    ));
                }
                if let (Some(a), Some(b), Some(c), Some(d)) = (vec[0].as_str(), vec[1].as_str(), vec[2].as_str(), vec[3].as_str()) {
                    (a, b, c, d)
                } else {
                    return Err(Error::from("The atoms names are not strings in dihedral potential"));
                }
            },
            None => {
                return Err(Error::from("Missing 'atoms' section in dihedral potential"));
            }
        };

        let potential = try!(read_dihedral_potential(potential));
        universe.add_dihedral_interaction(a, b, c, d, potential);
    }
    Ok(())
}

fn read_dihedral_potential(node: &Yaml) -> Result<Box<DihedralPotential>> {
    match node["type"].as_str() {
        Some(val) => {
            let val: &str = &val.to_lowercase();
            match val {
                "harmonic" => Ok(Box::new(try!(Harmonic::from_yaml(node)))),
                "cosine-harmonic" | "cosineharmonic" => Ok(Box::new(try!(CosineHarmonic::from_yaml(node)))),
                "torsion" => Ok(Box::new(try!(Torsion::from_yaml(node)))),
                "null" => Ok(Box::new(try!(NullPotential::from_yaml(node)))),
                val => Err(
                    Error::from(format!("Unknown potential type '{}'", val))
                ),
            }
        },
        None => {
            Err(Error::from(format!("Missing 'type' parameter in dihedral potential")))
        }
    }
}

/******************************************************************************/
impl FromYaml for NullPotential {
    fn from_yaml(_: &Yaml) -> Result<NullPotential> {
        Ok(NullPotential)
    }
}

impl FromYaml for Harmonic {
    fn from_yaml(node: &Yaml) -> Result<Harmonic> {
        if let (Some(k), Some(x0)) = (node["k"].as_str(), node["x0"].as_str()) {
            let k = try!(::units::from_str(k));
            let x0 = try!(::units::from_str(x0));
            Ok(Harmonic{k: k, x0: x0})
        } else {
            Err(
                Error::from("Missing 'k' or 'x0' in harmonic potential")
            )
        }
    }
}

impl FromYaml for LennardJones {
    fn from_yaml(node: &Yaml) -> Result<LennardJones> {
        if let (Some(sigma), Some(epsilon)) = (node["sigma"].as_str(), node["epsilon"].as_str()) {
            let sigma = try!(::units::from_str(sigma));
            let epsilon = try!(::units::from_str(epsilon));
            Ok(LennardJones{sigma: sigma, epsilon: epsilon})
        } else {
            Err(
                Error::from("Missing 'sigma' or 'espilon' in Lennard-Jones potential")
            )
        }
    }
}

impl FromYaml for CosineHarmonic {
    fn from_yaml(node: &Yaml) -> Result<CosineHarmonic> {
        if let (Some(k), Some(x0)) = (node["k"].as_str(), node["x0"].as_str()) {
            let k = try!(::units::from_str(k));
            let x0 = try!(::units::from_str(x0));
            Ok(CosineHarmonic::new(k, x0))
        } else {
            Err(
                Error::from("Missing 'k' or 'x0' in cosine harmonic potential")
            )
        }
    }
}

impl FromYaml for Torsion {
    fn from_yaml(node: &Yaml) -> Result<Torsion> {
        if let (Some(n), Some(k), Some(delta)) = (node["n"].as_i64(), node["k"].as_str(), node["delta"].as_str()) {
            let k = try!(::units::from_str(k));
            let delta = try!(::units::from_str(delta));
            Ok(Torsion{n: n as usize, k: k, delta: delta})
        } else {
            Err(
                Error::from("Missing 'n', 'k' or 'delta' in torsion potential")
            )
        }
    }
}

/******************************************************************************/
fn read_pair_computation(node: &Yaml, potential: Box<PairPotential>) -> Result<Box<PairPotential>> {
    match node["type"].as_str() {
        Some(val) => {
            let val: &str = &val.to_lowercase();
            match val {
                "cutoff" => Ok(Box::new(try!(CutoffComputation::from_yaml(node, potential)))),
                "table" => Ok(Box::new(try!(TableComputation::from_yaml(node, potential)))),
                val => Err(
                    Error::from(format!("Unknown computation type '{}'", val))
                ),
            }
        },
        None => {
            Err(Error::from(format!("Missing 'type' parameter for potential computation")))
        }
    }
}

trait FromYamlWithPairPotential where Self: Sized {
    fn from_yaml(node: &Yaml, potential: Box<PairPotential>) -> Result<Self>;
}

impl FromYamlWithPairPotential for CutoffComputation {
    fn from_yaml(node: &Yaml, potential: Box<PairPotential>) -> Result<CutoffComputation> {
        if let Some(cutoff) = node["cutoff"].as_str() {
            let cutoff = try!(::units::from_str(cutoff));
            Ok(CutoffComputation::new(potential, cutoff))
        } else {
            Err(
                Error::from("Missing 'cutoff' parameter in cutoff computation")
            )
        }
    }
}

impl FromYamlWithPairPotential for TableComputation {
    fn from_yaml(node: &Yaml, potential: Box<PairPotential>) -> Result<TableComputation> {
        if let (Some(n), Some(max)) = (node["n"].as_i64(), node["max"].as_str()) {
            let max = try!(::units::from_str(max));
            Ok(TableComputation::new(potential, n as usize, max))
        } else {
            Err(
                Error::from("Missing 'max' or 'n' parameter in cutoff computation")
            )
        }
    }
}

/******************************************************************************/
fn read_coulomb(universe: &mut Universe, config: &Yaml) -> Result<()> {
    let mut potential = try!(read_coulomb_potential(config));

    if config["restriction"].as_hash().is_some() {
        let restriction = try!(read_restriction(&config["restriction"]));
        potential.set_restriction(restriction);
    }

    universe.set_coulomb_interaction(potential);

    if config["charges"].as_hash().is_some() {
        try!(assign_charges(universe, &config["charges"]));
    }
    Ok(())
}

fn read_coulomb_potential(node: &Yaml) -> Result<Box<CoulombicPotential>> {
    match node["type"].as_str() {
        Some(val) => {
            let val: &str = &val.to_lowercase();
            match val {
                "wolf" => Ok(Box::new(try!(Wolf::from_yaml(node)))),
                "ewald" => Ok(Box::new(try!(Ewald::from_yaml(node)))),
                val => Err(Error::from(format!("Unknown coulomb solver type '{}'", val))),
            }
        },
        None => {
            Err(Error::from(format!("Missing 'type' parameter for coulomb section")))
        }
    }
}

impl FromYaml for Wolf {
    fn from_yaml(node: &Yaml) -> Result<Wolf> {
        if let Some(cutoff) = node["cutoff"].as_str() {
            let cutoff = try!(::units::from_str(cutoff));
            Ok(Wolf::new(cutoff))
        } else {
            Err(Error::from("Missing 'cutoff' parameter in Wolf potential"))
        }
    }
}

impl FromYaml for Ewald {
    fn from_yaml(node: &Yaml) -> Result<Ewald> {
        if let (Some(cutoff), Some(kmax)) = (node["cutoff"].as_str(), node["kmax"].as_i64()) {
            let cutoff = try!(::units::from_str(cutoff));
            if kmax < 0 {
                Err(Error::from("'kmax' can not be negative in Ewald potential"))
            } else {
                Ok(Ewald::new(cutoff, kmax as usize))
            }
        } else {
            Err(Error::from("Missing 'cutoff' or 'kmax' parameter in Ewald potential"))
        }
    }
}

fn assign_charges(universe: &mut Universe, config: &Yaml) -> Result<()> {
    let charges = config.as_hash().unwrap();
    for (name, charge) in charges {
        if let (Some(name), Some(charge)) = (name.as_str(), charge.as_f64()) {
            let mut n_changed = 0;
            for particle in universe.iter_mut() {
                if particle.name() == name {
                    particle.charge = charge;
                    n_changed += 1;
                }
            }
            if n_changed == 0 {
                return Err(Error::from(format!("No particle with the name {} was found", name)));
            } else {
                info!("Charge was set to {} for {} {} particles", charge, n_changed, name);
            }
        } else {
            return Err(
                Error::from(format!("Bad Yaml format in charges section: {:?}, {:?}", name, charge))
            );
        }
    }
    Ok(())
}

/******************************************************************************/
#[cfg(test)]
mod tests {
    use super::*;
    use universe::Universe;
    use std::path::{Path, PathBuf};
    use std::fs;

    fn bad_files(motif: &str) -> Vec<PathBuf> {
        let DATA = Path::new(file!()).parent().unwrap().join("data").join("bad");
        let paths = fs::read_dir(DATA).unwrap();

        // Convert the list of DirEntry to a list of PathBuf, and return only
        // the one whose filename starts with the given motif
        return paths.filter_map(|entry| entry.ok())
                    .map(|entry| entry.path())
                    .filter(|path| {
                        path.file_name().unwrap().to_str().unwrap().starts_with(motif)
                    }).collect();
    }

    #[test]
    fn pairs() {
        let DATA = Path::new(file!()).parent().unwrap().join("data");
        let mut universe = Universe::new();
        read_interactions(&mut universe, DATA.join("pairs.yml")).unwrap();
    }

    #[test]
    fn bad_pairs() {
        for path in bad_files("pairs") {
            let mut universe = Universe::new();
            assert!(read_interactions(&mut universe, path).is_err());
        }
    }

    #[test]
    fn bonds() {
        let DATA = Path::new(file!()).parent().unwrap().join("data");
        let mut universe = Universe::new();
        read_interactions(&mut universe, DATA.join("bonds.yml")).unwrap();
    }

    #[test]
    fn bad_bonds() {
        for path in bad_files("bonds") {
            let mut universe = Universe::new();
            assert!(read_interactions(&mut universe, path).is_err());
        }
    }

    #[test]
    fn angles() {
        let DATA = Path::new(file!()).parent().unwrap().join("data");
        let mut universe = Universe::new();
        read_interactions(&mut universe, DATA.join("angles.yml")).unwrap();
    }

    #[test]
    fn bad_angles() {
        for path in bad_files("angles") {
            let mut universe = Universe::new();
            assert!(read_interactions(&mut universe, path).is_err());
        }
    }

    #[test]
    fn dihedrals() {
        let DATA = Path::new(file!()).parent().unwrap().join("data");
        let mut universe = Universe::new();
        read_interactions(&mut universe, DATA.join("dihedrals.yml")).unwrap();
    }

    #[test]
    fn bad_dihedrals() {
        for path in bad_files("dihedrals") {
            let mut universe = Universe::new();
            assert!(read_interactions(&mut universe, path).is_err());
        }
    }

    #[test]
    fn coulomb() {
        let DATA = Path::new(file!()).parent().unwrap().join("data");
        let mut universe = Universe::new();
        read_interactions(&mut universe, DATA.join("wolf.yml")).unwrap();

        read_interactions(&mut universe, DATA.join("ewald.yml")).unwrap();
    }

    #[test]
    fn bad_coulomb() {
        for path in bad_files("coulomb") {
            let mut universe = Universe::new();
            assert!(read_interactions(&mut universe, path).is_err());
        }
    }
}
