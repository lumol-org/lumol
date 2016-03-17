# Welcome to Cymbalum user manual!

Cymbalum is a flexible and extensible molecular simulation engine. That means:

* You can use it to drive classical molecular simulations; from Monte-Carlo to
  molecular dynamics and energy minimization;
* You have full control over the simulation steps, and the simulated system;
* You can *easily* extend the engine with new algorithms, and hook into its
  internals.

This user manual describes how to use Cymbalum either as a command line tool or
a library.

## Systems and Simulations

The two main concepts in Cymbalum are a `System` and a `Simulation`. A `System`
contains all the simulated system data: atomic names, positions and velocities,
force-field, simulation box. On the other hand, a `Simulation` contains
everything needed to run a simulation with a system, *i.e.* algorithms and
associated data.

A `System` is usually used in combination with one `Simulation` to explore the
properties of the system. Sometimes multiple systems can be associated with the
same simulation (in Gibbs ensemble Monte-Carlo or in parallel tempering); and
sometimes multiple simulations can be associated with the same system, when
running an energy minimization before a molecular dynamics run.
