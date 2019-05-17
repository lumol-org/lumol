Input file reference
====================

This section describes how to use input files to run your simulations with
Lumol.  An input file contains all information that you need to run a simulation
and it is usually organized in four main sections: **input**, **log**,
**systems** and **simulations**.

-  The **input** section contains metadata about the input itself (i.e. a version number).
-  The **log** sections explains different methods about how Lumol reports information about your simulation such as warnings and errors.
-  The **systems** section contains information about the initial configuration, the interactions between atoms and the simulation cell.
-  The **simulations** section defines how your system will propagate. You can generally choose between molecular dynamics (MD), Monte-Carlo (MC) and energy minimization.

Lumol's input files use the `TOML <https://github.com/toml-lang/toml>`__
format, a simple and minimalist configuration format based on
``key = value`` pairs. You can read an introduction to the TOML format
`here <https://github.com/toml-lang/toml>`__.

**Example**:

.. code::

    # an example input file for a Monte Carlo simulation

    # input section
    [input]
    version = 1

    # log section
    [log]
    target = "lumol.log"

    # systems section
    [[systems]]
    file = "data/ethane.xyz"
    guess_bonds = true
    cell = 100.0

    [systems.potentials.global]
    cutoff = "14.0 A"
    tail_correction = true

    [systems.potentials.pairs.C-C]
    type = "lj"
    sigma = "3.750 A"
    epsilon = "0.814 kJ/mol"
    restriction = "InterMolecular"

    [systems.potentials.bonds]
    C-C = {type = "null"}

    # simulations section
    [[simulations]]
    nsteps = 1_000_000
    outputs = [
        {type = "Energy", file = "ethane_ener.dat", frequency = 500},
        {type = "Properties", file = "ethane_prp.dat", frequency = 500}
    ]

    [simulations.propagator]
    type = "MonteCarlo"
    temperature = "217.0 K"
    update_frequency = 500

    moves = [
        {type = "Translate", delta = "20 A", frequency = 50, target_acceptance = 0.5},
        {type = "Rotate", delta = "20 deg", frequency = 50, target_acceptance = 0.5},
        {type = "Resize", pressure = "5.98 bar", delta = "5 A^3", frequency = 2, target_acceptance = 0.5},
    ]


.. toctree::
   :hidden:

   input
   log
   systems/index
   interactions/index
   simulations/index
