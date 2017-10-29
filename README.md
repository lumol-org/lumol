# Lumol molecular simulation engine

[![Build Status](https://travis-ci.org/lumol-org/lumol.svg?branch=master)](https://travis-ci.org/lumol-org/lumol)
[![Coverage](https://codecov.io/gh/lumol-org/lumol/branch/master/graph/badge.svg)](https://codecov.io/gh/lumol-org/lumol)
[![Documentation](https://img.shields.io/badge/documentation-latest-brightgreen.svg)](https://lumol-org.github.io/lumol/latest/index.html)
[![Gitter](https://badges.gitter.im/lumol-org/lumol.svg)](https://gitter.im/lumol-org/lumol)

Lumol is a classical molecular simulation engine that provides a solid base for
developing new algorithms and methods. Using Lumol, you can customize the
behavior of all the algorithms in a simulation. Adding a new force field,
customizing Monte Carlo moves or molecular dynamics integrators is easy and well
documented.

Lumol goals are to be flexible, reliable and extensible. For us, this means that
this software should be:

- **flexible**: the code can simulate all kind of systems, from proteins to
  crystals, using various methods: molecular dynamics, Monte Carlo, *etc.*
- **reliable**: the code is well tested, both at the function level; and at the
  simulation level, checking thermodynamic properties of the systems;
- **extendable**: the code is modular, object-oriented, well documented,
  open-source, and easy to read.

Lumol is actively developed, and should be considered as alpha software. If
you are interested, have some questions or want to participate, you can open a
[Github issue][issues] or go to the project [chat room][Gitter].

## Features

- Pair, molecular and electrostatic interactions (with Ewald or Wolf methods);
- Energy minimization;
- Molecular dynamics simulations in the NVE, NVT and NPT ensembles;
- Monte Carlo simulations in the NVT ensemble;
- and many others! Have a look at the [documentation](#documentation) for more
  information

## Getting started

Lumol provides both a command line tool for running simulations; and a Rust
library for writing your own simulations algorithms using the pre-existing
building blocks.

### Documentation

Documentation is hosted [here](http://lumol-org.github.io/lumol), and separated
in multiple parts:

- The [user manual][user_manual] contains information about the general
  concepts of systems and simulations used in Lumol. Additionally, it has
  tutorials on how to use and extend Lumol. Use this documentation if you want
  to know basic concepts and how they are used in Lumol.
- The [input reference][input_reference] contains information about - well,
  the input file system of Lumol.
  Use this document if you want to use Lumol as a command line tool
  without writing code.
- To use Lumol as a library inside your own code, we have a [developer
  documentation][devdoc], which contains documentation for all the library
  public functions, and examples for most of them.

### Installation as a command line tool

You will need a stable Rust compiler, [grab one][Rust] if you do not have one
yet. Then, you can download the code, build it and install it by running:

```bash
cargo install --git https://github.com/lumol-org/lumol
```

This will produce the a `lumol` binary in `~/.cargo/bin`.

### Usage as a library

You can add Lumol as a dependency in your project's `Cargo.toml`:

```toml
[dependencies]
lumol = {git = "https://github.com/lumol-org/lumol"}
```

A tutorial about how to implement new algorithms in Lumol is coming. While
waiting, you can ask your questions [here][Gitter].

## Contributing

If you want to contribute to Lumol, there are several ways to go: improving the
documentation and helping with language issues; testing the code on your systems
to find bugs; adding new algorithms and potentials; providing feature requests.
Please come by and [talk with us][Gitter] a bit before staring new work, or open
an [issue][issues] to discuss improvements. We also have
[recommendations][contributing] for contributors.

See the [AUTHORS](AUTHORS) file for a list of contributors to the code.

## License

This software is licensed under the BSD license, see the LICENSE file for legal
text.

Unless you explicitly state otherwise, any contribution intentionally submitted
for inclusion in the work by you, shall be licensed under the same BSD license,
without any additional terms or conditions.

[Rust]: https://www.rust-lang.org/downloads.html
[Gitter]: https://gitter.im/lumol-org/lumol
[issues]: https://github.com/lumol-org/lumol/issues/new
[contributing]: Contributing.md
[user_manual]: http://lumol-org.github.io/lumol/latest/book/
[input_reference]: http://lumol-org.github.io/lumol/latest/book/
[devdoc]: http://lumol-org.github.io/lumol/latest/lumol/
