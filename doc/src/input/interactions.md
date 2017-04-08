# Interactions

Interactions describe the energies between atoms - or more general - between
interaction sites. These energies arise due to covalent bonding of atoms to
form molecules, short-range van der Waals or long-range Coulombic energies. To
describe interactions we use *potentials* which are functions  that take a set
of parameters and geometrical coordinates (for example distances  or angles) as
input and yield energies and forces. A set of potential functions  and
parameters to describe energies and forces of a molecule is also often called
*force field*. We will use the terms interactions and force field
interchangeably.

To be more specific, we distinguish between the following contributions:
  - `pairs` are van der Waals interactions between pairs of atoms;
  - `bonds` describe the energy between bonded atoms;
  - `angles` and `dihedrals` describe energy contributions due to bending and
  twisting of bonded atoms;
  - `coulomb` and `charges` describe long-range contributions due to
  electrostatic interactions;
  - the `global` section describes additional parameter that apply to all the
  energy contributions.

Information about interactions for `pairs`, `bonds`, `angles` and `dihedrals`
are organized as arrays of TOML tables. The `coulomb` section contains
information about the treatment of long-range electrostatic interactions and the
`charges` section defines the partial charges of the atoms.

An example of an input file for the f-SPC model of water is given bellow:

```toml
[input]
version = 1

[global]
cutoff = "10 A"

[[pairs]]
atoms = ["O", "O"]
lj = {sigma = "3.16 A", epsilon = "0.155 kcal/mol"}

[[pairs]]
atoms = ["O", "H"]
null = {}

[[pairs]]
atoms = ["H", "H"]
harmonic = {k = "79.8 kcal/mol/A^2", x0 = "1.633 A"}
restriction = "intra-molecular"

[[bonds]]
atoms = ["O", "H"]
harmonic = {k = "1054.2 kcal/mol/A^2", x0 = "1.0 A"}

[[angles]]
atoms = ["H", "O", "H"]
harmonic = {k = "75.9 kcal/mol/rad^2", x0 = "109.5 deg"}

[coulomb]
ewald = {cutoff = "10 A", kmax = 7}

[charges]
O = -0.82
H = 0.41
```

## van der Waals and covalent interactions

The `pairs`, `bonds`, `angles` and `dihedrals` sections are arrays, in which
every entry must contain at least two keys: the `atoms` key, and a
[potential](input/potentials.html) key. With the `atoms` key you can specify the
two atom types to which the interaction should be applied. The number of atoms
depends on the type of interaction: You have to provide two atoms for `pairs`
and `bonds`, three atoms for `angles` and four atoms for `dihedrals`.

For example you can use a `harmonic` bond potential for all your `C-H` bonds:

```toml
[[bonds]]
atoms = ["C", "H"]
harmonic = {x0 = "3.405 A", k = "2385 kcal/mol/A^2"}
```

The `[[pairs]]` entries can be customized further with a specific cutoff, or
pair restriction. See the [corresponding](input/pairs.html) documentation.

### Restrictions

In some force fields, neighboring atoms in molecules may interact solely via
covalent potentials (and not van der Waals or electrostatic). In this case, we
have to *exclude* interactions between neighbors in molecules. You can do that
(and more) in the `pairs` section by specifying a
[restriction](input/pairs.html#pairs-restrictions) using the `restriction` key.

```toml
[[pairs]]
atoms = ["Ar", "Ar"]
lj = {sigma = "3.405 A", epsilon = "0.2385 kcal/mol"}
cutoff = "10 A"
restriction  = "IntraMolecular"
```

## Coulombic interactions

The method for treatment of electrostatic interactions is specified in the
`coulomb` section. There are multiple available solvers for [electrostatic
interactions](input/electrostatic.html). Optionally, an additional
[restriction](input/pairs.html#pairs-restrictions) can be specified with the
`restriction` key.

```toml
[coulomb]
ewald = {cutoff = "10 A", kmax = 7}
restriction  = "exclude13"

[charges]
Na = 1
Cl = -1
```
