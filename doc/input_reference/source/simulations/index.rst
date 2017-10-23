*****************************
The ``[simulations]`` section
*****************************

The way to propagate the system is defined in the ``[[simulations]]``
section of the input file. This section always contains at least two
keys: ``nsteps`` specify the number of steps in the simulation, and the
``simulations.propagator`` table specify which propagator to use.

Here is an example of NPT molecular dynamics:

.. code::

    [[simulations]]
    nsteps = 1_000_000

    [simulations.propagator]
    type = "MolecularDynamics"
    timestep = "1 fs"
    integrator = {type = "BerendsenBarostat", pressure = "100 bar", timestep = 1000}
    thermostat = {type = "Berendsen", temperature = "400 K", timestep = 100}

Two propagators are currently implemented: one for `molecular
dynamics <md.html>`__ and one for `Monte Carlo <mc.html>`__.


.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: Contents:

   output
   mc
   md
   min
