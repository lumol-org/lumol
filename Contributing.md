So you want to contribute to Lumol and create the next generation of molecular
simulation engine? Great! This document explains how you can contribute, but in
case of doubt do not hesitate to contact us in [the chat room][Gitter]

# Bug reports and feature requests

Bug and feature requests should be reported as [Github issue][issues]. For bugs,
you should provide information so that we can reproduce it: what did you try?
What did you expect? What happened instead? Please provide any useful code
snippet or input file with your bug report.

# Contributing

A lot of contributions can be made to Lumol, and most of them do not involve
writing code:

- Improving the documentation and help with language issues;
- Testing the code on your systems to find bugs;
- Providing feedback on the interface;
- Adding new potentials;
- Adding new simulation algorithm;

All these contributions are very welcome. We accept contributions either via
Github pull request (have a look [here][PR] for Github model of pull request);
or you can send patches by email at luthaf@luthaf.fr.

If you want to work on the code and pick something easy to get started, have a
look at the [easy level issues][E-Easy].

If you want to add a new feature to Lumol, please create an issue so that we
can discuss it, and you have more chances to see your changes incorporated.

## Contribution check-list

Every item in this list is explained in the next section

- [ ] Fork Lumol;
- [ ] Create a local branch;
- [ ] Add code / correct typos / ...;
    - [ ] Add new tests with your code;
    - [ ] Add documentation for your code;
- [ ] Check that the tests still pass;
- [ ] Push to Github;
- [ ] Create a Pull-Request;
- [ ] Celebrate! :tada: :cake: :tada:

## Contribution tutorial

In this small tutorial, you should replace `<angle brackets>` as needed. If
anything is unclear, please [ask][Gitter] for clarifications! There are no dumb
questions.

---

Start by [forking Lumol][fork], and then clone and build your fork locally by
running:

```bash
git clone https://github.com/<YOUR USERNAME>/lumol
cd lumol
cargo build
```

Then create a new branch for your changes

```bash
git checkout -b <new-branch>
```

Do your changes, containing the documentation and associated unit tests for new
code.

Then, run all the test and check that they still pass:

```bash
# Unit tests
cargo test -all --lib
cargo test -all --doc
# integratop, tests, in relase mode to be faster
cargo test --all --tests --release
```

Finally, you can push your code to Github, and create a [Pull-Request][PR] to
the `lumol-org/lumol` repository.

We would appreciate if you run the rust formatting tool ([rustfmt]) on the code
before pushing the code to Github:

```bash
rustup toolchain add nightly
cargo +nightly install rustfmt
cargo +nightly fmt
```

[Gitter]: https://gitter.im/lumol-org/lumol
[issues]: https://github.com/lumol-org/lumol/issues/new
[PR]: https://help.github.com/articles/using-pull-requests/
[E-Easy]: https://github.com/lumol-org/lumol/labels/E-Easy
[fork]: https://help.github.com/articles/fork-a-repo/
[rustfmt]: https://github.com/rust-lang-nursery/rustfmt
