# Lumol documentation

Lumol has different documentations:

- The **user manual** contains information about the general concepts of systems
  and simulations used in Lumol. Additionally, it has tutorials on how to use
  and extend Lumol. Use this documentation if you want to know basic concepts
  and how they are used in Lumol.
- The **input reference** contains information about - well, the input file
  system of Lumol.  Use this document if you want to use Lumol as a command line
  tool without writing code.
- To use Lumol as a library inside your own code, we have a **developer
  documentation**, which contains documentation for all the library public
  functions, and examples for most of them.

Here is how you can build these documents. Both the user manual and the input
reference uses [sphinx] with [reStructuredText] markup, so you will need to
install some dependencies in order to build them:

```
pip install -r doc/requirements.txt
```

[sphinx]: http://www.sphinx-doc.org/en/stable/index.html
[reStructuredText]: www.sphinx-doc.org/en/stable/rest.html

## User manual

The user manual is located in the `book` directory.

```bash
cd book
make html
# HTML pages are in `build/html`
```

If you prefer other formats, *sphinx* has different options. You can inspect
all options by typing

```bash
make help
```

The documentation is found in `build/html/` (or `build/<out>/` where `<out>` is
the specified format for `make`).

## Input reference

The input reference is located in the `reference` directory.  It is also built
using sphinx:

```bash
cd reference
make html
# HTML pages are in `build/html`
```

## Programming interface documentation

In addition to the user manual, Lumol provides a complete programming interface
documentation, which you can build by running `cargo doc --open` in the root of
the repository.

This interface documentation is only useful if you want to use Lumol as a
library in you own code.
