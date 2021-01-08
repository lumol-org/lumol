Hello Sodium Chloride
=====================

In this example, we will simulate a Sodium Chloride crystal using molecular
dynamics to propagate the position throughout time. Sodium Chloride adds a
challenge in simulations because each atom carries a charge.  These charges
interact with a Coulomb potential which goes to zero as :math:`1 / r`. The
problem is that the cutoff scheme used for pair potentials in most molecular
simulations can not be used for anything that goes to zero slower than :math:`1
/ r^3`. So we need to use alternative methods to compute the potential for the
interactions between pairs of charges.

For this simulation, you will need the following files:

* the initial configuration ``nacl.xyz```
* the input file ``nacl.toml``

.. only:: html

    You can download both files :download:`here <../data/nacl.zip>`.

.. only:: latex

    You can download both files at the following URL:
    `<https://lumol.org/lumol/latest/book/_downloads/nacl.zip>`_.

Again, you can run the simulation which should complete in a minute with ``lumol
nacl.toml``. This will perform a molecular dynamics simulation of a NaCl crystal
using electrostatic interactions between atomic charges.

The input file commented
------------------------

We start with the input version again:

.. literalinclude:: ../data/nacl.toml
    :lines: 1-2

Then we load the system from the ``nacl.xyz`` file and define the unit cell.

.. literalinclude:: ../data/nacl.toml
    :lines: 4-6

Next, we define our potential. This section is way bigger than the one
for our previous Lennard-Jones example:

.. literalinclude:: ../data/nacl.toml
    :lines: 8-21

Let's break it down. First, we define some global values for the interactions:
setting ``systems.potentials.global.cutoff`` will use the given cutoff for all
pair interactions in the system.
We set the charges of atoms in the ``systems.potentials.charges`` section.

.. literalinclude:: ../data/nacl.toml
    :lines: 8-13

Then, we need to define the pair interactions for all possible pair combinations
in the system, *i.e.* (Na, Na), (Cl, Cl), and (Na, Cl).

.. literalinclude:: ../data/nacl.toml
    :lines: 15-18

Because our system contains charges, we need to use an electrostatic potential
solver. Here we are going for the ``Wolf`` solver, with a cutoff of 8 A.
Note that if we'd chose a cutoff different from the global one defined above,
we would overwrite the global one *for the coulombic interactions*.

.. literalinclude:: ../data/nacl.toml
    :lines: 20-21

We can now define the simulation and the outputs for this simulation.
We are using a molecular dynamics simulation of the NaCl crystal with a
timestep of 1 fs for integration (this will produce a NVE ensemble).

.. literalinclude:: ../data/nacl.toml
    :lines: 23-31

As before, run the simulation via

.. code-block:: bash

    lumol nacl.toml

Until now, the force field we used for the system was defined in the same input
file (in the ``system.potential`` section) as rest of the simulation settings.
It can sometimes be interesting to separate the force field from the rest of the
input, in particular when using the same force-field for multiple simulations.
In the :doc:`next example <water>`, we will do exactly this.
