# Interactions

The interactions in a system can be specified using [TOML] syntax. An
interaction is the association of atomic types, a potential function and
optionally restrictions and mode of computation.

[TOML]: https://github.com/toml-lang/toml

The `pairs`, `bonds`, `angles` and `dihedrals` sections are arrays of
interactions, listing possible atomic combinations and associated potentials.
The `coulomb` section contains information about the treatment of long-range
electrostatic interactions, and charges section defines the partial charges of
the atoms.

All the potentials input file must have the `input.potentials.version` key
defined to the value 1. This value will change if the input files formating
changes.

An example of an input file for the f-SPC model of water is given bellow:

```toml
[input.potentials]
version = 1

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

## Pairs and molecular interactions

The `pairs`, `bonds`, `angles` and `dihedrals` sections are arrays, in which
every entry must contain at least two keys: the `atoms` key, and a
[potential](input/potentials.html#Available%20potentials) key. The `atoms` key
specifies the atom types to which the interaction should be applied and should
contain two atoms for `pairs` and `bonds`, three atoms for `angles` and four
atoms for `dihedrals`.

For example, to use a harmonic bond potential for all `C-H` bonds:

```toml
[[bonds]]
atoms = ["C", "H"]
harmonic = {x0 = "3.405 A", k = "2385 kcal/mol/A^2"}
```

It is also possible in the `pairs` section to specify an additional
[restriction](input/potentials.html#Restrictions) using the `restriction` key;
or a [computation method](input/potentials.html#Potential%20computations) using
the `computation` key.

```toml
[[pairs]]
atoms = ["Ar", "Ar"]
lj = {sigma = "3.405 A", epsilon = "0.2385 kcal/mol"}
computation = {table = {n = 2000, max = "20.0 A"}}
restriction  = "IntraMolecular"
```

## Coulombic interactions

The method for treatment of electrostatic interactions is specified in the
`coulomb` section. There are multiple
[available](input/potentials.html#Electrostatic%20interactions) solvers for
electrostatic interactions. Optionally, an additional
[restriction](input/potentials.html#Restrictions) can be specified with the
`restriction` key.

```toml
[coulomb]
ewald = {cutoff = "10 A", kmax = 7}
restriction  = "exclude13"

[charges]
Na = 1
Cl = -1
```
