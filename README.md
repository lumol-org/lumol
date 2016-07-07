# Cymbalum: extensible molecular simulation engine

<div align="center">
[![Build Status](https://travis-ci.org/Luthaf/cymbalum.svg?branch=master)](https://travis-ci.org/Luthaf/cymbalum)
[![Coverage Status](https://codecov.io/github/Luthaf/cymbalum/coverage.svg?branch=master)](https://codecov.io/github/Luthaf/cymbalum?branch=master)
[![Documentation](https://img.shields.io/badge/documentation-latest-brightgreen.svg)](https://luthaf.github.io/cymbalum/latest/index.html)
[![Gitter](https://badges.gitter.im/Luthaf/cymbalum.svg)](https://gitter.im/Luthaf/cymbalum)
</div><br />

Cymbalum is a classical molecular simulation engine that provides a solid
base for developing new algorithms and methods.

Using Cymbalum, you can customize the behavior of all the algorithms in a
simulation (from force fields to barostats and Monte-Carlo moves). Cymbalum is

- **Easy to extend**: the code is modular, object-oriented, well documented,
  open-source and readable;
- **Easy to use**: the user interface is nice, with human-oriented input files.
- **Stable**: it will never crash on a good input, and try to provide helpful
  error messages.

Cymbalum is actively developed, and should be considered as alpha software. If
you are interested, have some question or want to participate, you can open a
Github issue or go to the project [chat room][Gitter].

[Gitter]: https://gitter.im/Luthaf/cymbalum

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

Cymbalum is usable both as a Rust library to write your own simulation code
using the provided bricks, and as a ready to use command line tool for running
simulations.

### Installation as a command line tool

[Rust]: https://www.rust-lang.org/downloads.html

You will need a stable Rust compiler, [grab one][Rust] if you do not have one
yet. Then, you can download the code, build it and install it by running:

```bash
cargo install --git https://github.com/Luthaf/cymbalum
```

This will produce the a `cymba` binary in `~/.cargo/bin`.

### Usage as a library & developement

You can also clone a build cymbalum locally, for developing it or using it as a
library by running:

```bash
git clone https://github.com/Luthaf/cymbalum
cd cymbalum
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

[NewIssue]: https://github.com/Luthaf/cymbalum/issues/new

### Documentation

Documentation is hosted [here](http://luthaf.github.io/cymbalum), and separated
in two main parts:

- The [user manual](http://luthaf.github.io/cymbalum/latest/book/) contains
  informations to use cymbalum as a command line tool, and the complete input
  file documentation. Use this documentation if you want to use Cymbalum as a
  simulation engine — without writing code.
- To use Cymbalum as a library inside your own code, we have a [developer
  documentation](http://luthaf.github.io/cymbalum/latest/cymbalum/), which
  should be complete with all functions documented.

## Contributing

If you want to contribute to Cymbalum, there are several ways to go: improving
the documentation and the usage of English language; testing the code on your
systems to find bugs; adding new algorithms and potentials; providing feature
requests. Please come by and [talk with us][Gitter] a bit before staring new
work, or open an [issue][NewIssue] to discuss improvements. We also have
[recomendation](https://github.com/Luthaf/cymbalum/blob/master/Contributing.md)
for contributors.

## License

This software is licensed under the BSD license, see the LICENSE file for legal
text.

Unless you explicitly state otherwise, any contribution intentionally submitted
for inclusion in the work by you, shall be licensed under the same BSD license,
without any additional terms or conditions.
