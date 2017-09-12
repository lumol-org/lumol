# Frequently Asked Questions

Here are some commons questions about Lumol. If you have more questions, please
contact us on [Gitter] to ask it, so that we can add it here!

- [What kind of simulation can I run with Lumol?](#what-kind-of-simulation-can-i-run-with-lumol)
- [Why should I use Lumol?](#why-should-i-use-lumol)
- [Why should I NOT use Lumol?](#why-should-i-not-use-lumol)
- [How can I build the initial configuration?](#how-can-i-build-the-initial-configuration)
- [Is the code parallel?](#is-the-code-parallel)
- [Why is Lumol written in Rust?](#why-is-lumol-written-in-rust)
- [Is there any graphical interface to Lumol?](#is-there-any-graphical-interface-to-lumol)

## What kind of simulation can I run with Lumol?

You should be able to run any kind of classical simulation, the only limit being
the number of atoms fitting in memory.

## Why should I use Lumol?

If any of these statement is true for you, you should consider using Lumol:

- You need to use a specific potential that is not yet available in other
  codes, or develop your own potential. Adding a new potential in Lumol is very
  simple and take less than 20 lines of code;
- You are developing new simulation algorithms, for example more efficient
  free-energy computations or better parallel scaling of Coulomb computations.
  Lumol allow you to write the specific algorithm, and reuse all the other part
  of the simulation engine;

Other nice goodies include:

- Nicely formatted and easy to read input files;
- *(and more to come ...)*

## Why should I NOT use Lumol?

Here are some reasons for you not to use Lumol:

- You need to get the fastest code for your simulations because you are working
  with a lot of atoms. Lumol is relatively young and is not yet fully optimized;
- You need to run your simulation on a cluster. Lumol can run on multiple cores
  (think OpenMP), but not yet on multiple nodes (think MPI).

## How can I build the initial configuration?

Lumol does not provide tools for building the initial simulation configuration.
There are already a lot of very good tools around, that you can use. Examples
include [VMD], [packmol], [gromacs], and many others. Because Lumol uses
[chemfiles] to read initial configuration, any [format supported by
chemfiles][formats] can be used.

[VMD]: http://www.ks.uiuc.edu/Research/vmd/
[packmol]: http://www.ime.unicamp.br/~martinez/packmol/
[gromacs]: http://gromacs.org/
[chemfiles]: http://chemfiles.org/
[formats]: http://chemfiles.org/chemfiles/latest/formats.html#list-of-supported-formats

## Is the code parallel?

Lumol can run in parallel on a single computer, using the multiple cores of the
processor (this is shared memory parallelism, like OpenMP). It is not yet
possible to run Lumol on multiple nodes in a cluster (message passing
parallelism, like MPI).

## Why is Lumol written in Rust?

[Rust] is a language created by Mozilla, and was released in 1.0 version
in may 2015. It is a modern language, that provides the same access to the bare
metal performances as C or C++, but prevents some programmer mistakes leading to
crashes and corruptions.

This allow to build better software faster, because the programmer does not need
to spend as much time debugging the code. At the same time, it also allow to
check at compile-time that a code is data-race free, and allow to build parallel
programs more easily.

## Is there any graphical interface to Lumol?

Not yet. But because Lumol is built as a library implementing all the simulation
algorithms, it should be relatively easy to create a graphical interface around
it. If you are interested in a graphical user interface (using it or building
it), please contact us!

[Gitter]: https://gitter.im/lumol-org/lumol
[Rust]: http://www.rust-lang.org/
