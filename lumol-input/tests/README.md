# Lumol-input tests

This directory contains the input tests for lumol. These tests check two things:

- that correct input files are parsed without errors;
- that incorrect input files generate the expeted error.

The test input files are either in `simulation` for the main input or in
`interactions` for the interactions inputs. In both directories, there is a
`good` and a `bad` directory of sample input files. The files in `bad` are
checked to check that the actual error message is the expeted one. The expected
error message can occur anywhere in the file, in a single line, starting with
`#^`:

```toml
[input]
version = 1

[[angles]]
atoms = ["A", "A"]
#^ Wrong size for 'atoms' array in angle potential. Should be 3, is 2
```

When multiple related errors are to be tested, files can contain multiple
standalone input, separated by `+++`:

```toml
[input]
version = 1

[[angles]]
atoms = ["A", "A"]
#^ Wrong size for 'atoms' array in angle potential. Should be 3, is 2

+++

[input]
version = 1

[[angles]]
atoms = ["A", "B", "A", "B"]
#^ Wrong size for 'atoms' array in angle potential. Should be 3, is 4
```
