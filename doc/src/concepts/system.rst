System
******

A ``System`` contains all the data about the physical system we are simulating.
It contains four types of data:

-  A list of **Particles**, which are physical objects with a position, a
   velocity, a mass and a name;
-  A list of **Molecules** containing information about how the particles are
   bounded together;
-  An **UnitCell**, *i.e.* the bounding box of the simulation.
-  **Interactions**, sometimes called a force field.

Particles
---------

Particles are the basic building blocks of a system. They can be atoms or more
complex: coarse-grained sites, dummy atoms, anisotropic particles ...

They have a name, a mass, a position, a velocity, and most importantly a
particle kind, defined by the name of the particle. All particles with the same
name share the same kind: all ``H`` are the same, as well as all ``C``, *etc.*

Molecules
---------

When particles are bonded together, they form molecules. A ``Molecule`` contains
the list of its bonds; molecules make for the *molecular* part in
*molecular simulation*.

Unit cells
----------

Lumol knows about three types of unit cells:

-  Infinite cells do not have any boundaries and can be used to simulate
   systems in vacuum;
-  Orthorhombic cells have up to three independent lengths whereas all angles
   of the cell are set to 90°;
-  Triclinic cells have 6 independent parameters: 3 lengths and 3 angles.

Orthorhombic and Triclinic cells are used in combination with `periodic boundary
conditions <pbc_>`_ to simulate infinite systems.

.. _pbc: https://en.wikipedia.org/wiki/Periodic_boundary_conditions

Interactions
------------

Interactions are potentials acting on or between particles.
Lumol provides functions for various potential types:

-  Non-bonding pair potentials;
-  Bonds potentials in molecules;
-  Angles potentials in molecules;
-  Dihedral angles potentials in molecules;
-  Long-ranges coulombic potentials (Ewald and Wolf methods);
-  Arbitrary external potential applying on the whole system at once.
