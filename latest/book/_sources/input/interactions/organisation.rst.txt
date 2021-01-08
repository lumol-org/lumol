Organizing interactions
=======================

You can specify interactions between atoms in two ways: either inside the main
input file or in a separate file. For example take a look at the input file for
the flexible SPC water model where we put all interactions directly into our
main input file:

.. code::

    # system configuration: initial positions, topology and cell
    [[systems]]
    file = "data/water.xyz"
    topology = "topology.pdb"
    cell = 40

    # intermolecular potentials
    [systems.potentials.pairs]
    O-O = {type = "lj", sigma = "3.165 A", epsilon = "0.155 kcal/mol"}
    H-H = {type = "null"}
    O-H = {type = "null"}

    # intramolecular potentials
    [systems.potentials.bonds]
    O-H = {type = "harmonic", k = "1059.162 kcal/mol/A^2", x0 = "1.012 A"}

    [systems.potentials.angles]
    H-O-H = {type = "harmonic", k = "75.90 kcal/mol/deg^2", x0 = "113.24 deg"}

    # ... additional interactions omitted

As you can see, there is a lot of bracket notation going on here. First, in
``[systems.potentials.xxx]``, the ``potentials`` key is actually a nested table
of ``systems`` indicated by the dot notation. Accordingly, ``pairs``, ``bonds``,
``angles``, etc. are nested tables of ``potentials``. Second, ``{type =
"harmonic", k = "75.90 kcal/mol/deg", x0 = "113.24 deg"}`` is the notation for
an inline table.

Input files can get very big and hard to read when you simulate complex systems
with a large number of different atoms. In these scenarios it may be better to
define a separate input file for your interactions like so:

.. code::

    [[systems]]
    file = "data/water.xyz"
    topology = "topology.pdb"
    cell = 40
    potentials = "water.toml"

Here, the ``potentials`` key contains a string that is interpreted as the path
to another input file containing only definitions of interactions. This way, you
can build your own library of force field files.
