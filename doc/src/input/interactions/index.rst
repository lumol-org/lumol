Interactions input
##################

Interactions describe the energies between atoms - or more general - between
interaction sites. These energies arise due to covalent bonding of atoms to form
molecules, short-range van der Waals or long-range Coulombic energies. To
describe interactions we use *potentials* which are functions that take a set of
parameters and geometrical coordinates (for example distances or angles) as
input and yield energies and forces. A set of potential functions and parameters
to describe energies and forces of a molecule is also often called *force
field*. We will use the terms interactions and force field interchangeably.

To be more specific, we distinguish between the following contributions:

- ``pairs`` are van der Waals interactions between pairs of atoms;
- ``bonds`` describe the energy between bonded atoms;
- ``angles`` and ``dihedrals`` describe energy contributions due to bending and
  twisting of bonded atoms;
- ``coulomb`` and ``charges`` describe long-range contributions due to
  electrostatic interactions;
- the ``global`` section describes additional parameter that apply to all the
  energy contributions.

Information about interactions for ``pairs``, ``bonds``, ``angles`` and
``dihedrals`` are organized as TOML tables. The ``coulomb`` section contains
information about the treatment of long-range electrostatic interactions and the
``charges`` section defines the partial charges of the atoms.

.. toctree::
   :maxdepth: 2

   organisation
   non_bonded
   electrostatic
   potentials
   restrictions
