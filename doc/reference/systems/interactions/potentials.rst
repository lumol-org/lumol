Available potentials
====================

This section is a list of all the available potentials in Lumol, with the
associated parameters. All potentials have to provide additional parameters in
there definition, as a TOML table. Using inline tables is the easiest way to do
so:

.. code::

    [[pairs]]
    atoms = ["A", "B"]
    # Additional parameters here are 'sigma' and 'epsilon'.
    lj = {sigma = "3 A", epsilon = "123 kJ/mol"}

The same potential can be used for either pairs (at distance :math:`r`); or for
angles (at angle :math:`\phi`). In all the formulas, the :math:`x` parameter
represents either a distance or an angle.

Null potential
--------------

This potential is always 0, for all values of :math:`x`. It should be used to
remove interactions between atoms in a pair/bond/angle/dihedral that are present
in the system but should not be interacting.

This potential can be used by specifying the ``null`` key with an empty table
``{}`` as value.

.. code::

    [[pairs]]
    atoms = ["O", "O"]
    null = {}

Lennard-Jones potential
-----------------------

The Lennard-Jones potential is a classical potential for pair interactions
expressed as:

.. math::

    V(x) = 4 \epsilon \left[\left(\frac{\sigma}{x}\right)^{12} -
   \left(\frac{\sigma}{x}\right)^6\right].

The Lennard-Jones potential is defined using the ``lj`` key. The parameters are
``sigma`` (:math:`\sigma`) and ``epsilon`` (:math:`\epsilon`), which should be
provided as strings.

.. code::

    [[pairs]]
    atoms = ["O", "O"]
    lj = {sigma = "3.16 A", epsilon = "0.155 kcal/mol"}

Buckingham potential
--------------------

The Buckingham potential is a potential for pair interactions expressed as:

.. math::


   V(x) = A \exp(-r / \rho) - \frac{C}{r^6}.

The potential type keyword is ``buckingham``, and the parameters ``A``, ``rho``
(:math:`\rho`) and ``C`` should be provided as strings.

.. code::

    [[pairs]]
    atoms = ["A", "B"]
    buckingham = {A = "40 kJ/mol", C = "120e-6 kJ/mol/A^6", rho = "3.0 A"}

Born-Mayer-Huggins potential
----------------------------

The Born-Mayer-Huggins potential is a potential for pair interactions, used in
particular for halide alkali. Its expression is:

.. math::

    V(x) = A
   \exp\left(\frac{\sigma -r}{\rho}\right) - \frac{C}{r^6} + \frac{D}{r^8}.

The potential type keyword is ``born``, and the parameters ``A``, ``C``, ``D``,
``sigma`` (:math:`\sigma`) and ``rho`` (:math:`\rho`) should be provided as
strings.

.. code::

    [[pairs]]
    atoms = ["A", "B"]
    [pairs.born]
    A = "40 kJ/mol"
    C = "120e-6 kJ/mol/A^6"
    D = "23e-6 kJ/mol/A^8"
    rho = "3.0 A"
    sigma = "2.2 A"

Harmonic potential
------------------

The Harmonic potential is usually used for intramolecular interactions such as
bonds, angles or dihedrals. It is expressed as:

.. math::  V(x) = \frac 12 k \ (x - x_0)^2

The potential type keyword is ``harmonic``, and the parameters are ``k`` and
``x0``, provided as strings.

.. code::

    [[bonds]]
    atoms = ["O", "H"]
    harmonic = {k = "1054.2 kcal/mol/A^2", x0 = "1.0 A"}

    [[angles]]
    atoms = ["H", "O", "H"]
    harmonic = {k = "75.9 kcal/mol/rad^2", x0 = "109.5 deg"}

Cosine-Harmonic potential
-------------------------

This potential is usually used for angles and dihedral angles interactions,
because it presents a :math:`2\pi` periodicity. It is expressed as:

.. math::

    V(x) = \frac 12
   k \ (\cos x - \cos x_0)^2

The potential type keyword is ``cosine-harmonic``, and the parameters ``k`` and
``x0`` should be provided as strings.

.. code::

    [[angles]]
    atoms = ["H", "C", "H"]
    cosine-harmonic = {k = "67 kJ/mol", x0 = "120 deg"}

Torsion potential
-----------------

This potential is usually used for dihedral interactions. It is expressed as:

.. math::  V(x) = k \ (1 + \cos(n x - \delta))

The potential type keyword is ``torsion``, and the parameters ``k`` and
``delta`` (:math:`\delta`) should be provided as strings, and ``n`` should be
provided as an integer.

.. code::

    [[dihedrals]]
    atoms = ["C", "C", "C", "C"]
    torsion = {k = "40 kJ/mol", delta = "120 deg", n: 4}

Morse potential
---------------

This potential is usually used for intramolecular interaction such as bonds,
angles or dihedrals. It is a better approximation for the vibrational structure
of the molecule than the Harmonic potential. It is expressed as:

.. math::

    V(r) = depth
   (1 - \exp(- a (r - x_0))^2

The potential type keyword is ``morse``, and the parameters ``a``, ``x0`` and
``depth`` should be provided as strings.

.. code::

    [[pairs]]
    atoms = ["A", "B"]
    morse = {depth = "40 kJ/mol", a = "2.0 A^-1", x0 = "1.3 A"}

For angles and dihedral angles, ``x0`` and ``a`` should be provided in angle
units:

.. code::

    [[pairs]]
    atoms = ["A", "B"]
    morse = {depth = "40 kJ/mol", a = "2.0 rad^-1", x0 = "109.7 deg"}


Gaussian potential
------------------

This potential is usually used to describe energy wells and is expressed as:

.. math::

    V(r) = -a \exp(-b r^2)

The potential type keyword is ``gaussian``, and the parameters ``a`` (well depth)
and ``b`` (well width) should be provided as strings.

.. code::

    [[pairs]]
    atoms = ["A", "B"]
    gaussian = {a = "8.0 kJ/mol", b = "0.2 A^-2"}

.. caution::
    ``b`` has to be positive
