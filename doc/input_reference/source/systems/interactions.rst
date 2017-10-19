Specifying interactions
-----------------------

Interactions between atoms are formulated via potentials; functions that
give us information about energies and forces acting between atoms. We
distinguish between *inter*- and *intramolecular* potentials. For an
overview of available potential functions, have a look
`here <../interactions/potentials.html>`__.

You can specify interactions between atoms in two ways: either inside
the main input file or in a separate file. As example we will use a
flexible SPC water model and put it directly into our main input file:

.. code::

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

As you can see, there is a lot of bracket notation going on here. First,
in ``[[systems.potentials.xxx]]``, the ``potentials`` key is actually a
nested table of ``systems`` indicated by the dot notation. Accordingly,
``pairs``, ``bonds``, ``angles``, etc. are nested tables of
``potentials``. Second,
``harmonic = {k = "75.90 kcal/mol/deg", x0 = "113.24 deg"}`` is the
notation for an inline table.

Input files can get very big and hard to read when you simulate complex
systems with a large number of different atoms. In these scenarios it
may be better to define a separate input file for your interactions like
so:

.. code::

    [[systems]]
    file = "data/water.xyz"
    topology = "topology.pdb"
    cell = 40
    potentials = "water.toml"

Here, the ``potentials`` key contains a string that is interpreted as
the path to another input file containing only definitions of
interactions. We will give examples of interaction input files
`here <input/interactions.html>`__ along with available keywords.
