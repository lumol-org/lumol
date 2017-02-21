# Systems

Let's talk about how you can set up your system. The system contains information
about:
- the configuration, i.e. the positions (and velocities) of your atoms;
- which atoms are connected (bonded) to form molecules;
- how atoms will interact with each other;
- and the simulation cell (i.e. volume).

All these details are listed after the `[[systems]]` keyword. The double
brackets indicate arrays of tables in TOML. Don't get confused too much, we will
talk about these in more detail while we go through the different parts of the
input file.


## Setting the initial configuration

A convenient way to get initial atom positions is by reading them from a
file using the `file` key:

```toml
[[systems]]
file = "data/water.pdb"
```

Lumol will read the file to build the system accordingly. If the file is a
trajectory containing multiple steps, the first frame is used. Under the hood,
we utilize [chemfiles](http://chemfiles.github.io/) to parse the data. You can
read about available file formats in the [chemfiles documentation][formats].

[formats]: http://chemfiles.readthedocs.io/en/latest/formats.html

## Setting the topology

You might want to use simple file formats such as `*.xyz` that don't specify
bonding information. You can do that by adding the bonding information
separately from your configuration file.

- Use the `topology` key to provide a file containing topological details.
  ```toml
  [[systems]]
  # the *.xyz format has no bonding information
  file = "water.xyz"
  # the *.pdb format specifies connections between atoms
  topology = "topology.pdb"
  ```

  where `topology.pdb` (here for water) may look like:

  ```
  HEADER    water
  COMPND
  SOURCE
  HETATM    1  H    ...
  HETATM    2  H    ...
  HETATM    3  O    ...
  CONECT    1    3
  CONECT    2    3
  CONECT    3    1    2
  END
  ```
  with `CONECT` entries detailing for bonding information.

  Often, you can get ready-to-use topology files from databases or create your
  own topologies. If you want to know more about the Protein Data Bank (PDB)
  format have a look at the [PDB website][PDB];

- The second option is to use the `guess_bonds` key to utilize a distance-based
  algorithm that will guess the bonds in the system. The algorithm is the same
  as the one in [VMD][VMD]. You should always check and confirm that the bonds
  guessed by the algorithm are coherent with your system.
  ```toml
  [[systems]]
  file = "water.xyz"
  guess_bonds = true
  ```

  Note that `guess_bonds` takes a boolean value as argument: there are no
  quotation marks around `true`. Also, TOML is case sensitive, i.e. writing
  `guess_bonds = True` will throw an error.

[PDB]: http://wwpdb.org/
[VMD]: http://www.ks.uiuc.edu/Research/vmd/

## Setting the unit cell

To set up the (initial) simulation cell you use the `cell` key.
We offer three different ways to set the cell:
- `cell = <length>` creates a cubic unit cell with the given side
  length. `<length>` should be a numeric value (no quotation marks) in Angstrom.
  ```toml
  [[systems]]
  file = "water.xyz"
  topology = "topology.pdb"
  cell = 40
  ```
- `cell = []` create an infinite unit cell, without boundaries. This can be
  used when periodic boundary conditions are undesirables, for example to
  simulate aggregates in the void;
  ```toml
  [[systems]]
  file = "water.xyz"
  topology = "topology.pdb"
  cell = []
  ```
- `cell = [<a>, <b>, <c>]` creates an orthorhombic unit cell.
  You should provide the lengths of the cell, `<a>`, `<b>`, and `<c>` as numeric
  values in Angstrom.
  ```toml
  [[systems]]
  file = "water.xyz"
  topology = "topology.pdb"
  cell = [24, 24, 76]
  ```
- `cell = [<a>, <b>, <c>, <alpha>, <beta>, <gamma>]` creates a triclinic unit
  cell with the given side lengths and angles. `<a>`, `<b>`, and `<c>`
  should be numeric values in Angstrom and `<alpha>`, `<beta>`, and `<gamma>`
  numeric values in degree.
  ```toml
  [[systems]]
  file = "water.xyz"
  topology = "topology.pdb"
  cell = [24., 24., 22., 90., 82.33, 110.4]
  ```

  Note that in an TOML array, all values have to have the same type.
  `cell = [24, 24, 76]` will work since we use all integer values, while
  `cell = [24., 24., 76]` will throw an error.

## Initializing velocities

For molecular dynamics (MD) simulations you need initial positions and initial
velocities of all atoms in your system. Use the `velocities` key to initialize
the velocities in the following way:

```toml
[[systems]]
file = "data/water.xyz"
topology = "topology.pdb"
velocities = {init = "300 K"}
```

where the `init` key will take the temperature as *string*. The velocities will
be initialized from a Boltzmann distribution at the given temperature.
Monte-Carlo simulations will not make any use of velocities since transition
probabilities (i.e. how the system evolves) are based on the positions (and the
underlying interactions) only.


## Specifying interactions

Interactions between atoms are formulated via potentials; functions that give us
information about energies and forces acting between atoms. We distinguish
between *inter*- and *intramolecular* potentials. For an overview of available
potential functions, have a look [here][potentials].

[potentials]: input/potentials.html

You can specify interactions between atoms in two ways: either inside the main
input file or in a separate file. As example we will use a flexible SPC water
model and put it directly into our main input file:

```toml
# system configuration: initial positions, topology and cell
[[systems]]
file = "data/water.xyz"
topology = "topology.pdb"
cell = 40

# intermolecular potentials
[[systems.potentials.pairs]]
atoms = ["O", "O"]
lj = {sigma = "3.165 A", epsilon = "0.155 kcal/mol"}

[[systems.potentials.pairs]]
atoms = ["H", "H"]
null = {}

[[systems.potentials.pairs]]
atoms = ["O", "H"]
null = {}

# intramolecular potentials
[[systems.potentials.bonds]]
atoms = ["O", "H"]
harmonic = {k = "1059.162 kcal/mol/A^2", x0 = "1.012 A"}

[[systems.potentials.angles]]
atoms = ["H", "O", "H"]
harmonic = {k = "75.90 kcal/mol/deg^2", x0 = "113.24 deg"}

# additional interactions omitted
```

As you can see, there is a lot of bracket notation going on here. First, in
`[[systems.potentials.xxx]]`, the `potentials` key is actually a nested table of
`systems` indicated by the dot notation. Accordingly, `pairs`, `bonds`,
`angles`, etc. are nested tables of `potentials`. Second, `harmonic = {k =
"75.90 kcal/mol/deg", x0 = "113.24 deg"}` is the notation for an inline table.

Input files can get very big and hard to read when you simulate complex systems
with a large number of different atoms. In these scenarios it may be better to
define a separate input file for your interactions like so:

```toml
[[systems]]
file = "data/water.xyz"
topology = "topology.pdb"
cell = 40
potentials = "water.toml"
```

Here, the `potentials` key contains a string that is interpreted as the path to
another input file containing only definitions of interactions. We will give
examples of interaction input files [here][interactions] along with available
keywords.

[interactions]: input/interactions.html
