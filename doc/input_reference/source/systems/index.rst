#######
Systems
#######

Let's talk about how you can set up your system. The system contains
information about: - the configuration, i.e. the positions (and
velocities) of your atoms; - which atoms are connected (bonded) to form
molecules; - how atoms will interact with each other; - and the
simulation cell (i.e. volume).

All these details are listed after the ``[[systems]]`` keyword. The
double brackets indicate arrays of tables in TOML. Don't get confused
too much, we will talk about these in more detail while we go through
the different parts of the input file.

.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: Contents:

   initial_configuration
   interactions

