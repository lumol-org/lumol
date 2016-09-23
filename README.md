# Lumol: extensible molecular simulation engine

[![Build Status](https://travis-ci.org/lumol-org/lumol.svg?branch=master)](https://travis-ci.org/lumol-org/lumol)
[![Coverage](https://codecov.io/gh/lumol-org/lumol/branch/master/graph/badge.svg)](https://codecov.io/gh/lumol-org/lumol)
[![Documentation](https://img.shields.io/badge/documentation-latest-brightgreen.svg)](https://lumol-org.github.io/lumol/latest/index.html)
[![Gitter](https://badges.gitter.im/lumol-org/lumol.svg)](https://gitter.im/lumol-org/lumol)

Lumol is a classical molecular simulation engine that provides a solid
base for developing new algorithms and methods.

Using Lumol, you can customize the behavior of all the algorithms in a
simulation (from force fields to barostats and Monte-Carlo moves). Lumol is

- **Easy to extend**: the code is modular, object-oriented, well documented,
  open-source and readable;
- **Easy to use**: the user interface is nice, with human-oriented input files.
- **Stable**: it will never crash on a good input, and try to provide helpful
  error messages.

Lumol is actively developed, and should be considered as alpha software. If
you are interested, have some question or want to participate, you can open a
Github issue or go to the project [chat room][Gitter].

[Gitter]: https://gitter.im/lumol-org/lumol

## Features

Currently, the following algorithms are implemented:
- Pair and molecular interactions;
- Generic methods for interactions evaluation;
- Support for user-defined potentials in an ergonomic way;
- Electrostatic interactions using Ewald or Wolf method;
- Energy minimization;
- Molecular dynamics simulations in the NVE, NVT and NPT ensembles;
- Monte-Carlo simulations in the NVT ensemble;

In a short-term, I'd like to add:
- Monte-Carlo moves for the NPT and µVT ensembles;
- Extended ensemble (Nose-Hoover) integration for molecular dynamics;

## Getting started

Lumol is usable both as a Rust library to write your own simulation code
using the provided bricks, and as a ready to use command line tool for running
simulations.

### Installation as a command line tool

[Rust]: https://www.rust-lang.org/downloads.html

You will need a stable Rust compiler, [grab one][Rust] if you do not have one
yet. Then, you can download the code, build it and install it by running:

```bash
cargo install --git https://github.com/lumol-org/lumol
```

This will produce the a `lumol` binary in `~/.cargo/bin`.

### Usage as a library & developement

You can also clone a build lumol locally, for developing it or using it as a
library by running:

```bash
git clone https://github.com/lumol-org/lumol
cd lumol
cargo build
```

Various usage are in the `examples` directory. You can build them by running
`cargo test --release --no-run`; and run them from the `examples` directory.

We also have unit and integration tests, that you can run with:

```bash
# Run only the unit test
cargo test --lib
# Run all the tests in release mode.
cargo test --release
```

Any failing test is an issue, please [report it][NewIssue]!

[NewIssue]: https://github.com/lumol-org/lumol/issues/new

### Documentation

Documentation is hosted [here](http://lumol-org.github.io/lumol), and separated
in two main parts:

- The [user manual](http://lumol-org.github.io/lumol/latest/book/) contains
  informations to use lumol as a command line tool, and the complete input
  file documentation. Use this documentation if you want to use Lumol as a
  simulation engine — without writing code.
- To use Lumol as a library inside your own code, we have a [developer
  documentation](http://lumol-org.github.io/lumol/latest/lumol/), which
  should be complete with all functions documented.

## Contributing

If you want to contribute to Lumol, there are several ways to go: improving
the documentation and the usage of English language; testing the code on your
systems to find bugs; adding new algorithms and potentials; providing feature
requests. Please come by and [talk with us][Gitter] a bit before staring new
work, or open an [issue][NewIssue] to discuss improvements. We also have
[recomendation](https://github.com/lumol-org/lumol/blob/master/Contributing.md)
for contributors.

## License

This software is licensed under the BSD license, see the LICENSE file for legal
text.

Unless you explicitly state otherwise, any contribution intentionally submitted
for inclusion in the work by you, shall be licensed under the same BSD license,
without any additional terms or conditions.
