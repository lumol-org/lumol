.. lumol_input_reference documentation master file, created by
   sphinx-quickstart on Thu Oct 19 16:57:49 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to the Lumol input reference!
=====================================

An input file contains all information that you need to run a
simulation. It is usually organized in three main sections: **input**,
**systems** and **simulations**.

-  The **input** section contains metadata about the input itself (i.e. a version number).
-  The **systems** section contains information about the initial configuration, the interactions between atoms and the simulation cell.
-  The **simulations** section defines how your system will propagate. You can generally choose between molecular dynamics (MD) or Monte-Carlo (MC).

Interactions (often also called *force field*) are the heart of every
simulation since they encapsulate the physical behaviour of atoms by
defining how they interact with each other. You can specify interactions
within the main input file as a part of the ``systems`` section, but
since your system can contain a huge number of interactions it is often
more convenient to create a separate input file. Doing so has two
advantages. First, it will keep your main input file short and readable
and second, you can simply reuse your force field input for different
simulations (you can even build your own force fields library).

Format
------

Lumol input files use the `TOML <https://github.com/toml-lang/toml>`__
format, a simple and minimalist configuration format based on
``key = value`` pairs. You can read an introduction to the TOML format
`here <https://github.com/toml-lang/toml>`__.

.. toctree::
   :hidden:

   input_section/index
   systems/index
   interactions/index
   simulations/index
   units

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
