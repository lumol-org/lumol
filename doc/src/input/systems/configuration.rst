Setting the initial configuration
---------------------------------

A convenient way to get initial atom positions is by reading them from a
file using the ``file`` key:

.. code::

    [[systems]]
    file = "data/water.pdb"

Lumol will read the file to build the system accordingly. If the file is a
trajectory containing multiple steps, only the first frame is used. Under the
hood, we utilize `chemfiles`_ to parse the data and thus we can read in plenty
different file formats. All possible formats are listed in the `chemfiles`_
documentation.

From the file, we will read in the unit cell, the atomic positions, the atomic
masses, and use the atomic types as particles names. We will also read the list
of bonds from the topology.

.. _chemfiles: http://chemfiles.org/


Initializing velocities
-----------------------

For molecular dynamics (MD) simulations you need initial positions and initial
velocities of all atoms in your system. If no velocities are present within the
read in configuration you can use the ``velocities`` key to initialize the
velocities in the following way:

.. code::

    [[systems]]
    file = "data/water.xyz"
    topology = "topology.pdb"
    velocities = {init = "300 K"}

where the ``init`` key will take the temperature as *string*. The velocities
will be initialized from a Boltzmann distribution at the given temperature.
Monte Carlo simulations will not make any use of velocities since transition
probabilities (i.e. how the system evolves) are based on the positions (and the
underlying interactions) only.

Setting the simulation cell
---------------------------

To set up the (initial) simulation cell you can use the `cell` key. This key is
only needed if the configuration file does not contain this information (for
example XYZ file), or if you want to override the cell from the file.

We offer three different ways to set the cell:

-  ``cell = <length>`` creates a cubic unit cell with the given side length.
   ``<length>`` should be a numeric value (no quotation marks) in Angstrom.

   .. code:

    [[systems]]
    file = "water.xyz"
    topology = "topology.pdb"
    cell = 40

- ``cell = []`` creates an infinite unit cell, without boundaries. This can be
  used when periodic boundary conditions are undesirable, for example to
  simulate aggregates in the void;

  .. code:

    [[systems]]
    file = "water.xyz"
    topology = "topology.pdb"
    cell = []

-  ``cell = [<a>, <b>, <c>]`` creates an orthorhombic unit cell.
   You should provide the lengths of the cell, ``<a>``, ``<b>``, and ``<c>`` as
   numeric values in Angstrom.

   .. code:

    [[systems]]
    file = "water.xyz"
    topology = "topology.pdb"
    cell = [24, 24, 76]

-  ``cell = [<a>, <b>, <c>, <alpha>, <beta>, <gamma>]`` creates a triclinic unit
   cell with the given side lengths and angles. ``<a>``, ``<b>``, and ``<c>``
   should be numeric values in Angstrom and ``<alpha>``, ``<beta>``, and
   ``<gamma>`` numeric values in degree.

   .. code:

    [[systems]]
    file = "water.xyz"
    topology = "topology.pdb"
    cell = [24., 24., 22., 90., 82.33, 110.4]

.. note::
    In an TOML array, all values have to have the same type.  ``cell = [24, 24,
    76]`` will work since we use all integer values, while ``cell = [24., 24.,
    76]`` will throw an error.
