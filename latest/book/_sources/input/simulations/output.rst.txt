Specifying output information
=============================

Additionally, a simulation can also output the evolution of the system
properties. Which properties are needed is given in the ``outputs`` array:

**Example**

.. code::

    [[simulations]]
    nsteps = 1_000_000
    outputs = [
        {type = "Trajectory", file = "filename.xyz", frequency = 100},
        {type = "Energy", file = "energy.dat", frequency = 200},
        {type = "Custom", file = "custom.dat", template = "{vx[3] / mass[3]}"},
    ]

    [simulations.propagator]
    ...

This array is an array of tables, containing three keys:

- the ``type`` of output
- the ``file`` to write the output to
- the ``frequency`` of the output.

The ``file`` is the path where the output will be written to.  The frequency is
a number and the output will be written every ``frequency`` steps to the file.
Except for the ``Trajectory`` output, all files are formatted with header lines
starting with a ``#``, and containing information about the quantities and the
units used for the output followed by multiple lines containing the step and
associated quantities.  The available outputs are the following:

-  The ``Energy`` output will write the potential, kinetic and total energy;
-  The ``Cell`` output will write the unit cell parameters, lengths and angles;
-  The ``Properties`` output will write the volume, the instant pressure
   (computed from the virial equation) and the instant temperature of the
   system;
-  The ``Stress`` output will write all the components of the stress tensor
   (computed from the virial equation);
-  The ``Trajectory`` output should be used to write a trajectory. The format of
   the trajectory will be guessed from the ``file`` extension.  Supported
   formats are documented in `chemfiles`_ documentation.
-  The ``Custom`` output is the most powerful one, taking an user-provided
   template string and using it to output data. The template should be given as
   a string with the ``template`` key in the TOML input file.


.. _chemfiles: http://chemfiles.org/

Here are some examples of custom output templates:

- A constant string is reproduced as it is: ``some data``;
- Anything in braces is replaced by the corresponding values: ``{pressure}
  {volume}`` will write the pressure and volume;
- Mathematical operators are allowed in braces: ``{pressure / volume}`` will
  print an entry containing the quotient of pressure and volume.  You can use
  ``+``, ``-``, ``/``, ``*``, ``^`` for exponentiation and parentheses;
- Some properties are arrays of atomic properties: ``{x[0] + y[20]}`` will print
  the sum of the x position of the 0'th atom and the y position of the 20'th
  position;
- Finally, all the properties are given in the internal units but one can
  specify a different unit: ``x[0] / nm``.


Here is a list of all properties available to custom outputs:

- Atomic properties: `x`, `y` and `z` for cartesian coordinates, `vx`, `vy` and
  `vz` for cartesian components of the velocity, `mass` for the atomic mass,
  `charge` for the atomic charge.
- Physical properties: `pressure`, `volume`, `temperature`, `natoms`, stress
  tensor components: `stress.xx`, `stress.yy`, `stress.zz`, `stress.xy`,
  `stress.xz`, `stress.yz`, simulation `step`.
- Unit Cell properties: `cell.a`, `cell.b`, `cell.c` are the unit cell vector
  lengths; `cell.alpha`, `cell.beta` and `cell.gamma` are the unit cell angles.
