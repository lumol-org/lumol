Electrostatic interactions
==========================

When some particles in a system are charged, they interact with a Coulomb
potential:

.. math::

    V(x) = \frac{Z_i Z_j}{4 \pi \epsilon r_{ij}},

where :math:`Z_{i,j}` are the net charges of the particles, :math:`r_{ij}` the
distance between them and :math:`\epsilon` the dielectric permittivity of the
current medium (usually the one of void).  Because this potential goes to zero
at infinity slower than :math:`1/r^3`, it can not be computed in periodic
simulations using a cutoff distance. This section present the available solvers
for electrostatic interactions.

In the input files, electrostatic interactions are specified with two sections:
the ``[charges]`` section sets the values of the charges of the atoms in the
system, and the ``[coulomb]`` section sets the solver to use for the
interaction.

Charge section
--------------

Charges for the particles in the system are set in a ``[charges]`` section in
the potential input file. This section should contain multiple ``name =
<charge>`` entries, one for each charged particle in the system.

.. code::

    # Some salt here
    [charges]
    Na = 1
    Cl = -1

Ewald solver
------------

Ewald's idea to compute electrostatic interactions is to split the interaction
into a short-range term which can be handled with a cutoff scheme; and a long
range term that can be computed using a Fourier transform. For more information about
the Ewald summation and its variants, see `[Frenkel2002]`_.

.. _[Frenkel2002]: http://dx.doi.org/10.1063/1.881812

The ``[coulomb]`` section for using an Ewald solver looks like this in the input
file:

.. code::

    [coulomb]
    ewald = {cutoff = "9 A", accuracy = 1e-5}

The ``cutoff`` parameter specifies the cutoff distance for the short-range and
long-range interactions splitting. The ``accuracy`` parameter is used to request
a relative relative error in forces, and should be smaller than 1. The
``accuracy`` is used to set the other Ewald parameters (alpha and kmax).

It is also possible to manually set the Ewald parameters:

.. code::

    [coulomb]
    ewald = {cutoff = "9 A", kmax = 7, alpha = "0.33451 A^-1"}

The ``kmax`` parameter gives the number of points to use in the reciprocal space
(the long-range part of interactions). Usually 7-8 is a good value for pure
water, for a very periodic charges distribution (like a crystal) a lower value,
such as 5 is sufficient, and for more heterogeneous system, higher values of
``kmax`` are needed. The ``alpha`` parameter specifies the width of the charges
spreading used to smooth the distribution in reciprocal space. A good value of
``alpha`` is one that satisfies :math:`\exp \left(-\alpha \frac L 2 \right) <<
1`. If only ``kmax`` is provided in the input file, the default value of
:math:`\pi / \text{cutoff}` is used for ``alpha``.

Wolf solver
-----------

The Wolf summation method is another method for computing electrostatic
interactions presented in `[Wolf1999]`_.  This method replaces the expensive
computation in reciprocal space from Ewald by a corrective term, and can be
expressed as a converging sum over the charged pairs in the system.

.. _[Wolf1999]: http://dx.doi.org/10.1063/1.478738

It is accessible using the ``wolf`` keyword in the input files:

.. code::

    [coulomb]
    wolf = {cutoff = "11 A"}

The only parameter is a ``cutoff``, which - as a rule of thumb - should be
larger than the corresponding cutoff from Ewald summation. For example, ``cutoff
= "11 A"`` should be suitable for pure water.

--------------

[Frenkel2002] Frenkel, D. & Smith, B. *Understanding molecular simulation.*
(Academic press, 2002).

[Wolf1999] Wolf, D., Keblinski, P., Phillpot, S. R. & Eggebrecht, J.  *Exact
method for the simulation of Coulombic systems by spherically truncated,
pairwise 1/r summation.* The Journal of Chemical Physics **110**, 8254 (1999).
