Restrictions
============

Some force fields define additional restrictions concerning which particles
should interact together and which ones should not. For example, sometimes
bonded particles should not interact through electrostatic potential, or some
interactions should only be taken in account for particles not in the same
molecule. The way to specify this is to use restrictions. Restrictions can be
used in two places: in the ``[pairs]`` section, and in the ``[coulomb]``
section. In both cases, they are specified with the ``restriction`` keyword, and
one of the possible values.

.. code::

    [pairs]
    O-O = {type = "lj", sigma = "3 A", epsilon = "123 kJ/mol", restriction = {scale14 = 0.5}}

    [coulomb]
    ewald = {cutoff = "8 A", kmax = 6}
    restriction = "intermolecular"

The possible values for ``restriction`` are:

* ``"intramolecular"`` or ``"intra-molecular"`` to act only on particles that
  are in the same molecule;
* ``"intermolecular"`` or ``"inter-molecular"`` to act only on particles that
  are **NOT** in the same molecule;
* ``"exclude12"`` to exclude particles directly bonded together;
* ``"exclude13"`` to exclude particles directly bonded together or forming an
  angle;
* ``"exclude14"`` to exclude particles directly bonded together; forming an
  angle or a dihedral angle;
* ``{scale14 = <scaling>}`` works like ``exclude13``, *i.e.* intramolecular
  interactions between three neighboring particles (connected by two bonds) will
  not be computed.  Additionally, interactions between the first and the forth
  (hence the ``14`` in ``scale14``) particle will be computed, but using scaled
  energies and forces. This simply means that the energies and forces are
  multiplied (linear scaling) by the given scaling factor, which must be between
  0 and 1.
