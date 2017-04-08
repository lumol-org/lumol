# Simulation

The flow of a simulation in Lumol is represented below.

![Simulation flow graph](static/img/simulation.svg#center)

First, it uses a system to setup any algorithms which uses system specific
information.

Then it propagate the system for one step, compute the physical properties as
needed, and output these properties to the hard drive.

If the simulation is finished (either the required number of steps has been
done, or come convergence criterion is reached) it returns the updated system.
If the simulation is not finished, the system is propagated for one more step.

While propagating the simulation, Lumol uses two type of algorithms: propagators
and output algorithms.

## Propagators

Propagators are at the heart of a Simulation. They have the responsibility to
update the system at each simulation step. Currently, three propagators exists:
a molecular dynamics one, a Monte-Carlo one and a minimizer, for energy
minimization.

## Output algorithms

Output algorithms have the responsibility to compute and output statistical data
about the simulated system. Temperature, energy, radial distribution functions
are some examples of output algorithms. The full list of currently implemented
output algorithms is available [here][Output].

[Output]: input/simulations.html#outputs
