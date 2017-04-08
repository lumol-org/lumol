# The Lumol user manual

Welcome to the Lumol user manual. In this book we teach you how to use Lumol to
set up and run classical molecular simulations. We designed Lumol to be
*flexible* and *extensible*; you are able to customize your simulation to suit
your needs and use it as a platform to implement your own algorithms and
customized potential functions in an *easy* way. You can use Lumol as a command
line tool as well as a library in your own code.

Excited? Then let's start with an overview of this book:
- [Installation](installation.html): Where to get and how to install Lumol.
- [The first simulations](): Step-by-step tutorials on
  how to perform basic molecular dynamics and Monte-Carlo simulations.
- [Simulation concepts](concepts/intro.html): Learn about the details. How are
  systems assembled? How do system propagators work?
- [Input files reference](input/intro.html): The complete reference for the
  input files needed to perform simulations;
- [Advanced Tutorials](): Extending and customizing Lumol.

**Note**: Lumol is actively developed and should be considered as alpha software.
As the code is likely to change so is this documentation.

## Systems and Simulations

Lumol mainly consists of two building blocks: the `System` and the `Simulation`.
The `System` contains general information about the system, such as atomic
names, positions and velocities, force-field, simulation box. On the other hand,
a `Simulation` contains everything needed to run a simulation with a system,
*i.e.* algorithms and associated data.

A `System` is usually used in combination with one `Simulation` to explore the
properties of the system. Sometimes multiple systems can be associated with the
same simulation (in Gibbs ensemble Monte-Carlo or in parallel tempering); and
sometimes multiple simulations can be associated with the same system, when
running an energy minimization before a molecular dynamics run.
