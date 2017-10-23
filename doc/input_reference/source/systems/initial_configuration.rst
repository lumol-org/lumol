Setting the initial configuration
---------------------------------

A convenient way to get initial atom positions is by reading them from a
file using the ``file`` key:

.. code::

    [[systems]]
    file = "data/water.pdb"

Lumol will read the file to build the system accordingly. If the file is
a trajectory containing multiple steps, only the first frame is used. Under
the hood, we utilize `chemfiles <http://chemfiles.github.io/>`__ to
parse the data and thus we can read in plenty different file formats.
All possible formats are listed in the `chemfiles documentation <http://chemfiles.readthedocs.io/en/latest/formats.html>`__.


Initializing velocities
-----------------------

For molecular dynamics (MD) simulations you need initial positions and
initial velocities of all atoms in your system. If no velocities are present
within the read in configuration you can use the ``velocities``
key to initialize the velocities in the following way:

.. code::

    [[systems]]
    file = "data/water.xyz"
    topology = "topology.pdb"
    velocities = {init = "300 K"}

where the ``init`` key will take the temperature as *string*. The
velocities will be initialized from a Boltzmann distribution at the
given temperature. Monte Carlo simulations will not make any use of
velocities since transition probabilities (i.e. how the system evolves)
are based on the positions (and the underlying interactions) only.
