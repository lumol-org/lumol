.. _minimization:

Minimization
============

You can run a minimization by setting the propagator ``type`` to
``Minimization``. The unique needed key is the ``minimizer`` algorithm to use
for this simulation; you can also optionally set the criteria for minimization
convergence.

.. code::

    [simulations.propagator]
    type = "Minimization"
    minimizer = {type = "SteepestDescent"}
    criteria = {energy = "1e-5 kJ/mol", force2 = "1e-5 kJ^2/mol^2/A^2"}

The single minimization algorithm implemented is the steepest descent algorithm,
that updates the coordinates of the atom following the energy gradient.

The minimization stops when the energy difference between the previous and the
current step is lower than the ``energy`` criterion, or when the maximal squared
norm of the atomic force is lower than the ``force2`` criterion.
