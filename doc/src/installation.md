# Installation

Lumol is written in [Rust][Rust] ([why?][why-rust]), and you will need a Rust
compiler to compile it. You can download one from [here][rust-download], or
using your package manager (Rust is available on Fedora, and with Homebrew on OS
X). Lumol also depends on some C++ libraries, so you will need a C++ compiler
and CMake to be installed.

Lumol is tested on Linux and OS X, and should build on Windows without any
issue. You will need a C++11 capable compiler on Windows (MSVC > 15 or Mingw
with gcc > 4.9). Be sure to pick the corresponding version of the Rust compiler.

When all the dependencies are installed on you system, you can install the
latest release of Lumol using:

```
cargo install lumol
```

and the latest development version with:

```
cargo install --git https://github.com/lumol-org/lumol
```

Both of these commands will download and install the command line interface for
Lumol in `~/.cargo/bin/lumol`, where `~` is your home directory. You may want to
add `~/.cargo/bin` to your PATH or move the `lumol` binary in another directory
accessible in your PATH.


[Rust]: https://www.rust-lang.org/
[why-rust]: faq.html#why-is-lumol-written-in-rust
[rust-download]: https://www.rust-lang.org/downloads.html
