# Systems

The system to use in a simulation must be specified in the `[[systems]]` section
of the input file. The minimal input file contains only the `file` key, giving
the path to a file to read to get the system.

```toml
[[systems]]
file = "data/CO2.pdb"
```

The program read the file, and use information from this file to build the
system. If the file is a trajectory containing multiple steps, the first frame
is used. File reading is provided by [chemfiles](http://chemfiles.github.io/),
and supported formats are documented [here][formats].

[formats]: http://chemfiles.readthedocs.io/en/latest/formats.html

Additional keys allow to change the system after it is read from the file.

## Setting the topology

When using a simple `.xyz` file as the system definition, no bonding information
are present in the system. There are two ways to add bonding informations:

- Use the `topology` key to give a file containing topological informations
  ```toml
  [[systems]]
  file = "water.xyz"
  topology = "topology.pdb"
  ```
- Use the `guess_bonds` key, to use a distance-based algorithm to guess the
  bonds in the system. The algorithm used is the same as the one in [VMD][VMD].
  You should always check if the bonding informations guessed by the algorithm
  is coherent with your system.
  ```toml
  [[systems]]
  file = "water.xyz"
  guess_bonds = true
  ```

[VMD]: http://www.ks.uiuc.edu/Research/vmd/

## Setting the unit cell

The `cell` key can be used to set the unit cell of the system in three different
ways:
- `cell = <length>` will set the unit cell to a cubic cell with the given side
  length. `<length>` should be a numeric value in Angstrom.
  ```toml
  [[systems]]
  file = "water.xyz"
  cell = 40
  ```
- `cell = [<a>, <b>, <c>]` will set the unit cell to a orthorhombic cell with   
  the given side lengths. `<a>`, `<b>`, and `<c>` should be numeric values in
  Angstrom.
  ```toml
  [[systems]]
  file = "water.xyz"
  cell = [24, 24, 76]
  ```
- `cell = [<a>, <b>, <c>, <alpha>, <beta>, <gamma>]` will set the unit cell to a
  triclinic cell with the given side lengths and angles. `<a>`, `<b>`, and `<c>`
  should be numeric values in Angstrom; and `<alpha>`, `<beta>`, and `<gamma>`
  numeric values in degree.
  ```toml
  [[systems]]
  file = "water.xyz"
  cell = [24, 24, 22, 90, 82.33, 110.4]
  ```

## Initializing velocities

For molecular dynamics simulations, the `velocities` key can be used to
initialize the velocities in the system. The syntax is the following:

```toml
[[systems]]
file = "data/CO2.pdb"
velocities = {init = "300 K"}
```

Where `init` key refers to a temperature string. The velocities will the be
initialized from a Boltzmann distribution at this temperature.

## Specifying interactions

Interactions between atoms can be specified in two ways: either inside the main
input file, or as a separated file.

If the `potentials` key contains a string, this string is interpreted as the
path to another input file containing only interactions definitions, as
documented [here][interactions].

```toml
[[systems]]
file = "polymer.pdb"
potentials = "potentials.toml"
```

If the `potentials` key is actually a sub-table of `systems`, then the different
sections for interactions input (`pairs`, `charges`, `angles`, ...) are looked
inside this table. They use the same syntax as in [standalone
input][interactions].

```toml
[[systems]]
file = "polymer.pdb"

[[systems.potentials.pairs]]
atoms = ["C", "O"]
lj = {sigma = "3 A", epsilon = "0.5 kJ/mol"}

[[systems.potentials.angles]]
atoms = ["O", "C", "O"]
harmonic = {x0 = "109 deg", k = "5 kJ/mol/deg"}
```

[interactions]: input/interactions.html
