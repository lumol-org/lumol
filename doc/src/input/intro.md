# Input files

An input file contains all information that you need to run a simulation.
It is usually organized in three main sections: **input**, **systems** and
**simulations**.

- The [input](input/intro.html#input-metadata) section contains metadata about
  the input itself (i.e. a version number).
- The [systems](input/systems.html) section contains information about the
  initial configuration, the interactions between atoms and the simulation cell.
- The [simulations](input/simulations.html) section defines how your system will
  propagate. You can generally choose between molecular dynamics (MD) or Monte-
  Carlo (MC).

Interactions (often also called *force field*) are the heart of every simulation
since they encapsulate the physical behaviour of atoms by defining how they
interact with each other. You can specify interactions within the main input
file as a part of the `systems` section, but since your system can contain a
huge number of interactions it is often more convenient to create a separate
input file. Doing so has two advantages. First, it will keep your main input
file short and readable and second, you can simply reuse your force field input
for different simulations (you can even build your own force fields library). We
talk more about standalone input files for interactions on this
[page](input/interactions.html).

## Format

Lumol input files use the [TOML][TOML] format, a simple and minimalist
configuration format based on `key = value` pairs. You can read an introduction
to the TOML format [here][TOML].

[TOML]: https://github.com/toml-lang/toml

## Input metadata

All input files must contain an `[input]` section looking like this:

```toml
[input]
version = 1
```

Introducing a `version` key helps us to make changes to the input file format
while keeping compatibility with previous formats. Please note that
Lumol is not in version 1.0 yet and we currently cannot guarantee compatibility
for input files.

The input files can also contain a `[log]` section to control where should the
code output be printed. Please see the corresponding
[documentation](input/log.html) for more information.

## Units in input

The unit of a value can be defined by a specific string, which will be parsed
and converted to the [internal unit system](concepts/units.html).
If there is no unit in the string, the internal unit for this type is used.
No consistency check is performed, and it is up to you to check the given units.

```toml
# Here, 'cutoff' is a distance

# OK
cutoff = "8 A"
# OK, will use the internal unit of distance.
# This is not recommended. The internal unit may change, and the input convey less information
cutoff = "8"
# OK, but big !
cutoff = "8 km"

# OK, but probably not what you want. Will be interpreted as 8000 A
cutoff = "8 ps"
# invalid, 'cutoff' must be a string
cutoff = 8.0
```
