# Installation

Lumol is written in the Rust programming language. To compile it, we need a
stable Rust compiler. We can simply install Lumol using *cargo*, the Rust
package manager. You can use Lumol as a command line tool or as a library.


## Installation as command line tool

Via cargo:

```
cargo install --git https://github.com/Luthaf/lumol
```

This will produce the binary named `cymba` in ~/.cargo/bin.

## Installation as library

To use Lumol as library the easiest way is to clone it from github.

```
git clone https://github.com/Luthaf/lumol
cd lumol
cargo build
```

We provide some examples in the examples directory. You can build them by
running `cargo test --release --no-run`. To run examples, change to the examples
directory and run the binaries.

We also have unit and integration tests that you can run.

```
# Run only the unit test
cargo test --lib
# Run all the tests in release mode.
cargo test --release
```

Any failing test is an issue, please report it!
