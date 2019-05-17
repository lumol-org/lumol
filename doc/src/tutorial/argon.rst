Monte Carlo simulation of Argon
===============================

So let's run out first simulation with Lumol. The easiest system to simulate is
a Lennard-Jones fluid, which is a good model for noble gas fluids. Here we
will simulate super-critical argon using the Metropolis Monte Carlo algorithm.

For this simulation, you will need the following files:

* the initial configuration ``argon.xyz``
* the input file ``argon.toml``

.. only:: html

    You can download both files :download:`here <../data/argon.zip>`.

.. only:: latex

    You can download both files at the following URL:
    `<https://lumol.org/lumol/latest/book/_downloads/argon.zip>`_.

After extracting the archive, you can run the simulation with ``lumol
argon.toml``.  The simulation should complete in a few seconds and produce two
files: ``energy.dat`` and ``trajectory.xyz``.

Input file anatomy
------------------

The input file is written using the TOML syntax, you can learn more about this
`syntax here <https://github.com/toml-lang/toml>`__. The file starts with a
header declaring the version of the input file syntax used, here the version 1:

.. literalinclude:: ../data/argon.toml
    :lines: 1-2

Then, we declare which system that we want to simulate in the ``systems`` array.
We define this system using an XYZ file for which we also have to provide the
unit cell size.

.. literalinclude:: ../data/argon.toml
    :lines: 4-6

We also need to define the interaction potential between the atoms in the
system, which we'll do in the ``system.potential.pairs`` section.
Here, we are using a Lennard-Jones potential with a cutoff distance of 10
Angstrom for all Argon (Ar) pairs.

.. literalinclude:: ../data/argon.toml
    :lines: 8-9

Next, we define how we want to simulate our system. For this brief tutorial
simulating ``100000`` steps should be enough. We also specify some output:
the energy is written to ``energy.dat`` every 100 steps, and the trajectory
is written to ``trajectory.xyz`` also every 100 steps.

.. literalinclude:: ../data/argon.toml
    :lines: 11-16

Finally, we define how to propagate the system from one step to another.
We are using a Monte Carlo simulation in this tutorial. We need to specify the
temperature (set to 500 K) and the set of moves that we want to perform to
change the systems. The only Monte Carlo move in this example is a translation
for which we set the maximum amplitude (the furthest a particle can be
translated in a single move) to 1 Angstrom.

.. literalinclude:: ../data/argon.toml
    :lines: 18-23

Wrapping all this together, here is the complete input file:

.. literalinclude:: ../data/argon.toml

As mentioned above, you can start the simulation using

.. code-block:: bash

    lumol argon.toml

So now we know how to run a simulation of a Lennard-Jones fluid. How about we
add electrostatic interactions in the :doc:`next example <nacl>`?
