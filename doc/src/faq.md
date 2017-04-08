# FAQ

Here are some commons questions about Lumol. If you have more questions, please
contact us on [Gitter][Gitter] to ask it, so that we can add it here!

## Is there any graphical interface to Lumol?

Unfortunately, no. But because Lumol is built around a core library implementing
all the simulation algorithms, it should be easier to create a graphical
interface around it. If you are interested in a graphical user interface (using
it or building it), please contact us!

## What kind of simulation can I run with Lumol?

You should be able to run any kind of classical simulation, from 2 atom up to as
atoms fits in your computer memory. Lumol try very hard not to be biased toward
some systems, and is as flexible as possible.

## How can I build the initial configuration?

Lumol do not provide tools for building the initial simulation configuration.
There are already a lot of very good tools around, that you can use. Examples
include [VMD][VMD], [packmol][packmol], [gromacs][gromacs], and many others.

[VMD]: http://www.ks.uiuc.edu/Research/vmd/
[packmol]: http://www.ime.unicamp.br/~martinez/packmol/
[gromacs]: http://gromacs.org/

## Is the code parallel? If so, how well does it scale?

The code is currently exclusively running in serial mode. But is it written in
such a way porting it to use parallel algorithms should be very easy. If you
have trivially parallel computations to run (running the same simulation on
multiple systems, or on the same system but changing parameters; parallel
tempering); you should already be able to assemble the building block in Lumol
to run it. The usage of Rust guarantee us that the code is data-race free.

## Why should I use Lumol?

If any of these statement is true for you, you should consider using Lumol:

- You need to use a specific potential that is not yet available in other
  codes, or develop your own potential. Adding a new potential in Lumol is very
  simple and take less than 20 lines of code;
- You are developing new simulation algorithms, for example more efficient
  free-energy computations or better parallel scaling of Coulomb computations.
  Lumol allow you to write the specific algorithm, and reuse all the other part
  of the simulation engine;
- You are developing new simulation algorithms, for more efficient free-energy
  computations or better parallel scaling of Coulomb computations. Lumol allow
  you to write the specific algorithm, and reuse all the other part of the
  simulation engine;

Other nice goodies include:

- Nicely formatted and easy to read input files;
- *(and more to come ...)*

## Why should I NOT use Lumol?

Here are some reasons for you not to use Lumol:

- You need to get the fastest code for your simulations because you are working
  with a lot of atoms. Lumol is relatively young and is not yet fully optimized;
- You need to scale on a lot of processors, to make your computations faster.
  Lumol is [not yet parallel][parallel], and will be made so, but you should
  continue to use other codes in the meantime.

[parallel]: faq.html#is-the-code-parallel-if-so-how-well-does-it-scale

## Why is Lumol written in Rust?

[Rust][Rust] is a language created by Mozilla, and was released in 1.0 version
in may 2015. It is a modern language, that provides the same access to the bare
metal performances as C or C++, but prevents some programmer mistakes leading to
crashes and corruptions.

This allow to build better software faster, because the programmer does not need
to spend as much time debugging the code. At the same time, it also allow to
check at compile-time that a code is data-race free, and allow to build parallel
programs more easily.

[Gitter]: https://gitter.im/lumol-org/lumol
[Rust]: http://www.rust-lang.org/
