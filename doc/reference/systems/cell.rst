Setting the simulation cell
---------------------------

To set up the (initial) simulation cell you use the `cell` key.
We offer three different ways to set the cell:

-  ``cell = <length>`` creates a cubic unit cell with the given side
   length. ``<length>`` should be a numeric value (no quotation marks) in Angstrom.

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
   You should provide the lengths of the cell, ``<a>``, ``<b>``, and ``<c>`` as numeric
   values in Angstrom.

   .. code:

     [[systems]]
     file = "water.xyz"
     topology = "topology.pdb"
     cell = [24, 24, 76]
-  ``cell = [<a>, <b>, <c>, <alpha>, <beta>, <gamma>]`` creates a triclinic unit
   cell with the given side lengths and angles. ``<a>``, ``<b>``, and ``<c>``
   should be numeric values in Angstrom and ``<alpha>``, ``<beta>``, and ``<gamma>``
   numeric values in degree.

   .. code:

     [[systems]]
     file = "water.xyz"
     topology = "topology.pdb"
     cell = [24., 24., 22., 90., 82.33, 110.4]

.. note::
    In an TOML array, all values have to have the same type.
    ``cell = [24, 24, 76]`` will work since we use all integer values, while
    ``cell = [24., 24., 76]`` will throw an error.
