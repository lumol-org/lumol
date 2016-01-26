# Cymbalum: molecular simulation framework

<div align="center">
[![Build Status](https://travis-ci.org/Luthaf/cymbalum.svg?branch=master)](https://travis-ci.org/Luthaf/cymbalum)
[![Coverage Status](https://codecov.io/github/Luthaf/cymbalum/coverage.svg?branch=master)](https://codecov.io/github/Luthaf/cymbalum?branch=master)
[![Documentation](https://img.shields.io/badge/documentation-latest-brightgreen.svg)](http://luthaf.github.io/cymbalum/cymbalum/index.html)
</div><br />

Cymbalum is a classical molecular simulation framework that provides a solid base for
developing new algorithms and methods.

Using Cymbalum, you can customize the behavior of all the algorithms in a simulation (from
force fields to barostats and Monte-Carlo moves). Cymbalum is
- **Easy to extend**: the code is modular, object-oriented, well documented, open-source
  and readable;
- **Easy to use**: the user interface is nice, with human-oriented input files.
- **Stable**: it will never crash on a good input, and try to provide helpful error
  messages.

It also try to be (or will be) fast, parallel and interactive. That's a lot of goals,
is'nt it? So let's define non-goals. Cymbalum will not provide analysis code, initial
simulation building or a graphical user interface. There are a lot of good existing tools
for that.

Cymbalum is actively developed, and should be considered as alpha software. If you want to
participate, please contact me (by email, in a Github issue, ...)!

## Features

Currently, the following algorithms are implemented:
- Pair interactions and molecular interactions;
- Generic methods for interactions evaluation;
- Support for user-defined potentials in an ergonomic way;
- Coulombic interactions using Ewald or Wolf method;
- Energy minimization;
- Molecular dynamics simulations in the NVE, NVT and NPT ensembles;

In a short-term, I'd like to add:
- Monte-Carlo simulations in the NVT, NPT, and µVT ensembles;
- Extended ensemble integration for molecular dynamics;

## Getting started

### Installation

You will need a stable [Rust](https://www.rust-lang.org) compiler, [grab
one](https://www.rust-lang.org/downloads.html) if you do not have one yet. Then, clone
this repository and build the code with

```bash
git clone https://github.com/Luthaf/cymbalum
cd cymbalum
cargo build --release
```

This will produce the library in `target/release/libcymbalum.rlib`, and binary examples
in `target/release/examples`.

In order to run the unit and integration tests, you can use these commands

```bash
# Run only the unit test
cargo test --lib

# Run all the tests. Use the --release flag for faster tests
cargo test
```

Any failing test is an issue, please [report
it](https://github.com/Luthaf/cymbalum/issues/new)!

### Documentation

The [API documentation](http://luthaf.github.io/cymbalum/cymbalum/index.html) should be
extensive, please fill [an issue](https://github.com/Luthaf/cymbalum/issues/new) if there
are parts of it which are not clear.

The user documentation is nonexistent at the time, you can have a look at the examples
folder to get a taste of the usage.

## Contributing

If you want to contribute to Cymbalum, there are several ways to go: improving the
documentation and the usage of English language; testing the code on your systems to find
bugs; adding new algorithms and potentials; providing feature requests. Feel free to
contact me to discuss about your contributions! 

## License

All this code brought to you under the terms of the Mozilla Public License v2.0.
The documentation is subject to CC-BY-SA license. By contributing to Cymbalum, you
agree that your contributions are released under the corresponding licenses.
