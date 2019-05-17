Molecular dynamics simulation of water
======================================

In this example, we will simulate a molecular liquid of high scientific
interest: water. At the same time we will learn how we can reuse a potential
definition between multiple simulations.

You will need three input files for this simulation:

- the initial configuration ``water.xyz``;
- the simulation input ``water.toml``;
- the force field (potential definitions) ``water-fSCP.toml``.

.. only:: html

    You can download them :download:`here <../data/water.zip>`.

.. only:: latex

    You can download these files at the following URL:
    `<https://lumol.org/lumol/latest/book/_downloads/water.zip>`_.

As usual, you can run the simulation with ``lumol water.toml``. The simulation
should finish in a few minutes.

The input files commented
-------------------------

``water.toml`` file
^^^^^^^^^^^^^^^^^^^

The main input file is pretty similar to the previous examples, with two
novelties:

-  The ``guess_bonds = true`` entry tells lumol to try to guess the bonds from
   the distances in the XYZ file. This is needed because we want to simulate a
   molecule but there is no bonding information inside the XYZ format. If we
   were to use a PDB file with connectivity instead, this would not be needed;
-  The ``potentials = "water-fSCP.toml"`` tells lumol to read the potentials
   from the ``water-fSCP.toml`` file. Using this type of input for the
   potentials we can reuse the same potential for multiple simulations and also
   easily change the potential used in the simulation.

.. literalinclude:: ../data/water.toml

``water-fSCP.toml`` file
^^^^^^^^^^^^^^^^^^^^^^^^

``water-fSCP.toml`` is a standalone potential input file. It contains the same
data as a potential definition inside a ``[[systems]]`` section, but without the
``system.potential`` prefix on all section names.

In this input, we start by defining which version of the input format we are
using.

.. literalinclude:: ../data/water-fSCP.toml
    :lines: 1-2

Then, we can define some global input data: the pair potential cutoff and the
atomic charges.

.. literalinclude:: ../data/water-fSCP.toml
    :lines: 4-9

The pair potential section contains the usual declarations for pairs, with a few
additional options.

.. literalinclude:: ../data/water-fSCP.toml
    :lines: 11-12

We can add a ``restriction`` to restrict a specific pair interaction to some
kind of pairs. Here, we will only account for H-H pairs inside the same molecule
(``"intra-molecular"`` interactions).

.. literalinclude:: ../data/water-fSCP.toml
    :lines: 13

We can also define a non-interacting pair interaction! This is useful when some
atoms do not interact in the model we use.

.. literalinclude:: ../data/water-fSCP.toml
    :lines: 14

Next, the bonds and angles are defined. These are interactions used
between bonded particles (or angles formed by two bonds and dihedral angles
formed by three bonds). This section follows the same pattern as the
``[pairs]`` section.

.. literalinclude:: ../data/water-fSCP.toml
    :lines: 16-20

Finally we specify how to compute the electrostatic interaction, this time using
the Ewald summation method. We can restrict the coulombic interactions to only
apply between particles not in the same molecule by using
``restriction = "inter-molecular"``.

.. literalinclude:: ../data/water-fSCP.toml
    :lines: 22-24

The simulation is run via

.. code::

    lumol water.toml
