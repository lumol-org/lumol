# Lumol documentation

## User manual

This directory contains the HTML user manual for the Lumol molecular
simulation engine. You need `cargo` to build it, like this:

```bash
cargo install mdbook
mdbook build doc
# HTML pages are in `doc/book/`
```

This book is written using standard [markdown](http://commonmark.org/help/), a
simple markup language.

If you want to edit the book, you can use `mdbook watch doc` to rebuild it at
every change.

## Programming interface documentation

In addition to the user manual, Lumol provides a complete programming
interface documentation, which you can build by running

```
cargo doc --open
```

This interface documentation is only useful if you want to use Lumol as a
library in you own code.
