Simulation
**********

A simulation in Lumol always contains the same steps:

1. Setup the system and initialize the algorithms using system specific
   information;
2. Propagate the system for one step using a Propagator;
3. Compute the physical properties of the system and output them to a file;
4. Check if the simulation is finished (either the required number of steps has
   been done, or a convergence criterion is reached) and return the updated
   system. If the simulation is not finished, go back to (2)

Propagators
-----------

Propagators are at the heart of a Simulation. They have the responsibility to
update the system at each simulation step. Currently, three propagators exists:
a molecular dynamics one, a Monte Carlo one and a minimizer, for energy
minimization.

Output algorithms
-----------------

Output algorithms have the responsibility to compute and output statistical data
about the simulated system. Temperature, energy, radial distribution functions
are some examples of output algorithms.
