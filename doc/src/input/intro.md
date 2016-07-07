# Input files

Cymbalum input files uses the [TOML][TOML] format, a simple and minimalist
configuration format based on `key = value` pairs. You can read an introduction
to the TOML format [here][TOML].

[TOML]: https://github.com/toml-lang/toml

The input file describe everything needed for running a simulation: which system
to use, with which force field, and how to propagate the simulation. They
usually contain three main sections:
- The [input](input/intro.html#Input%20metadata) section describe metadata about
  the input itself;
- The [systems](input/systems.html) section describe the system to use in the
  simulation;
- The [simulations](input/simulations.html) section describe how to update the
  system during the simulation;

Interactions between particles in the system are a bit special: they can either
be specified in they own input file (to be reused in multiple simulations); or
be part of the `[[systems]]` section. This [page](input/interactions.html)
describe the standalone input for interactions.

## Input metadata

All input files must contain an `[input]` section looking like this:

```toml
[input]
version = 1
```

The purpose of the `version` key is to make changes to the input file format,
while keeping compatibility with the previous input format. Please note that
while Cymbalum have not reach version 1.0, no guarantee is made on input file
compatibility.

## Units in input

The unit of a value can be defined by a specific string, which will be parsed
and converted to the [internal unit system](input/units.html).
If there is no unit in the string, the internal unit for this type is used.
No consistency check is performed, and it is up to the code users to check the
given units.

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
