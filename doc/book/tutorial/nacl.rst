Hello Sodium Chloride
=====================

In this example, we will simulate a Sodium Chloride crystal, using molecular
dynamics to propagate the position throughout time. Sodium Chloride add a
challenge in simulation because each atom carry a charge.  These charges
interact with a Coulomb potential which goes to zero as :math:`1 / r`. The
problem is that the cutoff scheme used for pair potential in most molecular
simulations can not be used for anything that goes to zero slower than :math:`1
/ r^3`. So we need to use alternate methods to compute the potential for the
charges-charges interactions.

For this simulation, you will need the following files:

* the initial configuration ``nacl.xyz```
* the input file ``nacl.toml``

You can download both files :download:`here <../data/nacl.zip>`. Again you can
run the simulation which should complete in a minute with ``lumol nacl.toml``.

This will perform a molecular dynamics simulation of a NaCl crystal, using
electrostatic interactions between the atomic charges.

The input file commented
------------------------

We start with the input version again:

.. literalinclude:: ../data/nacl.toml
    :lines: 1-2

Then we load the system from the ``nacl.xyz`` file and define the unit cell.

.. literalinclude:: ../data/nacl.toml
    :lines: 4-6

Then comes the potential definition section, which is way bigger than the one
for our previous Lennard-Jones example:

.. literalinclude:: ../data/nacl.toml
    :lines: 8-28

Let's break it down. First we define some global values for the interactions:
setting ``systems.potentials.global.cutoff`` will use the given cutoff for all
the pair interactions. The ``systems.potentials.charges`` section defined the
atomic charges in the system.

.. literalinclude:: ../data/nacl.toml
    :lines: 8-13

Then, we need to define the pair interactions for all the pair combinations in the
system, *i.e.* (Na, Na); (Cl, Cl); and (Na, Cl).

.. literalinclude:: ../data/nacl.toml
    :lines: 15-25

Because our system have charges, we need to use an electrostatic potential
solver. Here we are going for the ``Wolf`` solver, with a cutoff of 8 A.

.. literalinclude:: ../data/nacl.toml
    :lines: 27-28

We can now define the simulation and the outputs for this simulation.

.. literalinclude:: ../data/nacl.toml
    :lines: 30-34

We are using here a molecular dynamics simulation of the NaCl crystal, and a
timestep of 1 fs for integration.

.. literalinclude:: ../data/nacl.toml
    :lines: 36-38
