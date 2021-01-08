Non-bonded interactions
=======================

Cutoff treatment
----------------

When computing the energy and forces for non-bonded pair interactions, Lumol
uses a cutoff radius :math:`r_c`. This means that the force and energy
associated with any pair at a distance bigger than :math:`r_c` will be zero. We
can use two different cutoff schemes, presented in the following section.

In the potentials input file, the cutoff should be specified for all the pairs.
It can be specified once for all the pairs in the ``global`` section, and then
overridden for specific interactions:

.. code::

    # Use a simple cutoff with a radius of 10 A for all pair interactions
    [global]
    cutoff = "10 A"

    [pairs]
    A-A = {type = "lj", sigma = "3 A", epsilon = "123 kJ/mol"}
    B-B = {type = "lj", sigma = "3 A", epsilon = "123 kJ/mol"}

    # Except for this one, use a shifted cutoff with a radius of 8 A
    [pairs.A-B]
    type = "lj"
    sigma = "3 A"
    epsilon = "123 kJ/mol"
    cutoff = {shifted = "8 A"}

Simple cutoff (potential truncation)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This scheme just sets the potential :math:`U(r)` to zero for any distance bigger
than :math:`rc`:

.. math::

    V(r) = \begin{cases}
       U(r) & r <= rc \\\\
       0 & r > rc
   \end{cases}

To use this potential truncation, we specify a string containing the cutoff
distance as the ``cutoff`` value.

.. code::

    [global]
    cutoff = "10 A"

    [pairs]
    O-O = {type = "lj", x0 = "3 A", k = "5.9 kJ/mol/A^2", cutoff = "8 A"}

Truncation with energy shift
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The potential :math:`U` can be additionally shifted to make sure it is
continuous at :math:`r = rc`. This is important for molecular dynamics, where a
discontinuity means an infinite force in the integration.

.. math::

    V(r) = \begin{cases}
       U(r) - U(rc) & r <= rc \\\\
       0 & r > rc
   \end{cases}

In the input, this uses a table containing the shifted value, which must be a
string containing the cutoff radius.

.. code::

    [global]
    cutoff = {shifted = "8 A"}

    [pairs]
    O-O = {type = "lj", x0 = "3 A", k = "5.9 kJ/mol/A^2", cutoff = {shifted = "10 A"}}

Tail correction
~~~~~~~~~~~~~~~

Tail corrections (also called long range corrections) are a way to account for
the error we introduce by cutting off the potential. If the simulated system is
homogeneous and isotropic beyond the the cutoff (if the pair distribution
function :math:`g(r)` is 1 after the cutoff) we have an expression for the
corrections we can evaluate. For a potential :math:`V(r)` with associated forces
:math:`\vec f(r)`, the missing energy and virial (which is used to compute
instantaneous pressure) expressions are at density :math:`\rho`:

.. math::  U_\text{tail} = 2 \pi \rho \int_{rc}^{+\infty} r^2 V(r) \ dr,

.. math::  P_\text{tail} = 2 \pi \rho^2 \int_{rc}^{+\infty} r^2 \vec r \cdot \vec f(r) \ dr.

In the input, these additional energetic and pressure terms are controlled by
the ``tail_correction`` keyword, which can be placed either in the ``[global]``
section, or in any specific ``[pairs]`` section.

.. code::

    # Use tail corrections for every pair interaction
    [global]
    tail_correction = true

    # Except for this one.
    [pairs]
    O-O = {type = "lj", x0 = "3 A", k = "5.9 kJ/mol/A^2", tail_correction = false}

Potentials computation
----------------------

The same potential function (Lennard-Jones, Harmonic, *etc.*) can be computed
with different methods: directly, by shifting at the cutoff distance, using a
table interpolation, *etc.* This is the purpose of computation. The default way
is to use the mathematical function corresponding to a potential to compute it.
To use a different type of computation, the ``computation`` keyword can be used
in the ``[pairs]`` section.

Table interpolation
-------------------

Another way to compute a potential is to compute it on a regularly spaced grid,
and then to interpolate values for points in the grid. In some cases, this can
be faster than recomputing the function every time.

This can be done with the ``table`` computation, which does a linear
interpolation in regularly spaced values in the ``[0, max)`` segment.  You need
to provide the ``max`` value, and the number of points ``n`` for the
interpolation:

.. code::

    [pairs.O-O]
    type = "lj"
    sigma = "3 A"
    epsilon = "123 kJ/mol"
    computation = {table = {max = "8 A", n = 5000}}
