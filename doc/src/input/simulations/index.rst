``[[simulations]]`` section
***************************

The way to propagate the system is defined in the ``[[simulations]]`` section of
the input file. This section always contains at least two keys: ``nsteps``
specify the number of steps in the simulation, and the
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

Three propagators are currently implemented:

- A :ref:`minimization` propagator, to minimize energy of a system before
  running another propagator;
- A :ref:`molecular-dynamics` propagator;
- A :ref:`monte-carlo` propagator;


.. toctree::
   :maxdepth: 2

   output
   min
   md
   mc
