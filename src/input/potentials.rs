/* Cymbalum, Molecular Simulation in Rust - Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
 */
#![allow(unused_variables)]

use std::io::prelude::*;
use std::fs::File;
use std::path::Path;

use ::universe::Universe;
use ::potentials::*;

use super::yaml::{YamlLoader, Yaml};
use super::{Error, Result};

trait FromYaml where Self: Sized {
    fn from_yaml(node: &Yaml) -> Result<Self>;
}

// TODO: Add support for * in atoms sections.
// TODO: Add support for computation section

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
/// dihedrals: # Molecular angles
///   - atoms: [O, C, C, O]
///     type: harmonic
///     k: 42 kJ/mol/deg^2
///     x0: 180 deg
/// ```
///
/// The main items are the `"pairs"`, `"bonds"`, `"angles"` and `"dihedrals"`;
/// which contains the data for pair potentials, bonds, angles and dihedral
/// angles potential. This data is orgnised as a vector of hash maps. Each maps
/// has at least a vector of `"atoms"` defining which particles will get this
/// interaction applied to, and a `"type"` parameter defining the type of the
/// potential. Others keys may be supplied depending on the potential type.
pub fn read_potentials<P: AsRef<Path>>(universe: &mut Universe, path: P) -> Result<()> {
    let mut file = try!(File::open(path));
    let mut buff = String::new();
    try!(file.read_to_string(&mut buff));

    let doc = &try!(YamlLoader::load_from_str(&buff))[0];
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

    Ok(())
}

/// Read the "pairs" or the "bonds" section in the file. If `pair_potentials`
/// is `true`, then the interactions are added to the pair interactions. Else,
/// the interaction are added to the bond interactions.
fn read_pairs(universe: &mut Universe, pairs: &Vec<Yaml>, pair_potentials: bool) -> Result<()> {
    for potential in pairs {
        let (a, b) = match potential["atoms"].as_vec() {
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

        let potential = try!(read_pair_potential(potential));
        if pair_potentials {
            universe.add_pair_interaction(a, b, potential);
        } else {
            universe.add_bond_interaction(a, b, potential);
        }
    }
    Ok(())
}

fn read_pair_potential(node: &Yaml) -> Result<Box<PairPotential>> {
    match node["type"].as_str() {
        None => {
            return Err(Error::from("Missing 'type' section in pair potential"));
        }
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
        }
    }
}

/******************************************************************************/

fn read_angles(universe: &mut Universe, angles: &Vec<Yaml>) -> Result<()> {
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
        None => {
            return Err(Error::from(
                format!("Missing 'type' section in: {:?}", node)
            ));
        }
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
        }
    }
}

/******************************************************************************/

fn read_dihedrals(universe: &mut Universe, dihedrals: &Vec<Yaml>) -> Result<()> {
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
        None => {
            return Err(Error::from(
                format!("Missing 'type' section in: {:?}", node)
            ));
        }
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
        }
    }
}

/******************************************************************************/

impl FromYaml for NullPotential {
    fn from_yaml(node: &Yaml) -> Result<NullPotential> {
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
