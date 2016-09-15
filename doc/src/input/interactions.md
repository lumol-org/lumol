# Interactions

An interaction is the association of atomic types, a potential function and
optionally restrictions and mode of computation. Input files for interaction
contains all this information, ordered in multiple sections. The `pairs`,
`bonds`, `angles` and `dihedrals` sections are arrays of interactions, listing
possible atomic combinations and associated potentials. The `coulomb` section
contains information about the treatment of long-range electrostatic
interactions, and the `charges` section defines the partial charges of the
atoms.

An example of an input file for the f-SPC model of water is given bellow:

```toml
[input]
version = 1

[global]
cutoff = "12.5 A"

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

Within the `pairs` section, additional information about how to compute the
pair potential can be given:

- A [restriction](input/potentials.html#Restrictions) restrict the potential
  to a subset of the pairs in the system with the `restriction` keyword;
- A different [computation method](   
  input/potentials.html#Potential%20computations) can be specified with the
  `computation` keyword;
- A different [cutoff treatment](
  input/interactions.html#Cutoff%20treatment%20for%20pair%20interactions) is
  specified by the `cutoff` keyword.

```toml
[[pairs]]
atoms = ["Ar", "Ar"]
lj = {sigma = "3.405 A", epsilon = "0.2385 kcal/mol"}
computation = {table = {n = 2000, max = "20.0 A"}}
restriction  = "IntraMolecular"
```

### Cutoff treatment for pair interactions

Two different cutoff treatments are available for pair interactions: a simple
cutoff distance, or a shifted potential with a cutoff. With a simple cutoff the
potential is only set to zero after the cutoff distance. With a shifted
potential, the potential energy is shifted so that it is zero at the cutoff
distance and after. This mean that the energy is continuous for the shifted
potential, which is a desirable property for molecular dynamics stability.

In the input files, a simple cutoff can be used with `cutoff = "8 A"`, and a
shifted potential can be used with `cutoff = {shifted = "8 A"}`. The cutoff treatment can be set either globally in the `[global]` table:

```toml
[input]
version = 1

[global]
cutoff = "8 A"

[[pairs]]
# ...
# All pairs interaction will use a 8 A cutoff, unless specified otherwise.
```

or specifically for each pair interaction:

```toml
[input]
version = 1

[[pairs]]
cutoff = "8 A"
# ...

[[pairs]]
cutoff = {shifted: "7 A"}
# ...

[[pairs]]
cutoff = "12 A"
# ...
```

If a global cutoff is defined, defining another cutoff in a `pairs` section
will override the global cutoff.

## Coulombic interactions

The method for treatment of electrostatic interactions is specified in the
`coulomb` section. There are multiple available solvers for [electrostatic
interactions](input/potentials.html#Electrostatic%20interactions). Optionally,
an additional [restriction](input/potentials.html#Restrictions) can be
specified with the `restriction` key.

```toml
[coulomb]
ewald = {cutoff = "10 A", kmax = 7}
restriction  = "exclude13"

[charges]
Na = 1
Cl = -1
```
