# Input files

The input files in Cymbalum can be classified in three types:
 - Initial configuration of the system;
 - Interactions and force field description;
 - Simulation setup and configuration.

Currently, only the initial configuration of the system and the force field
input files are implemented. The initial configuration can be provided in
various formats. Check the [relevant documentation](input/initial.html) for a
list and description of the formats.

The force field input file uses the YAML format, which is described below. The
different entries in the input file are described in the [next
section](input/interactions.html).

## Yaml input format

The Yaml format is a configuration format which is easy to read by humans. It
uses significant whitespace (like Python), and is based on `key: value` pairs.
If you are not familiar with the Yaml syntax, [this
page](https://learnxinyminutes.com/docs/yaml/) will teach you the basics.

<!-- TODO: add an introduction to YAML -->

In addition to the Yaml syntax specification, Cymbalum uses some conventions:
all the keys must be in lower case, but string values are case insensitive. You
should use `CamelCase` for strings — *i.e.* you should use `LennardJones` over
`lennardjones` — it is easier to read and prevents issues with the string
`Null`.

```yaml
# This is OK and equivalent
bar: Foo
bar: FOO
bar: foo

# This is not OK
Bar: Foo
BAR: Foo
```

## Units in input

The unit of a value can be defined by a specific string, which will be parsed
and converted to the [internal unit system](input/units.html).
If there is no unit in the string, the internal unit for this type is used.
No consistency check is performed, and it is up to the code users to check the
given units.

```yaml
# cutoff is a distance
cutoff: 8 A  # OK
cutoff: "8"  # OK
cutoff: 8 m  # OK, but big !

cutoff: 8 ps # OK, but probably not what you want. Will be interpreted as 8000 A
cutoff: 8    # invalid, must be a Yaml string
```
