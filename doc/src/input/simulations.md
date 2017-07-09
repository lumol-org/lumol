# Simulations

The way to propagate the system is defined in the `[[simulations]]` section of
the input file. This section always contains at least two keys: `nsteps` specify
the number of steps in the simulation, and the `simulations.propagator` table
specify which propagator to use.

Here is an example of NPT molecular dynamics:

```toml
[[simulations]]
nsteps = 1_000_000

[simulations.propagator]
type = "MolecularDynamics"
timestep = "1 fs"
integrator = {type = "BerendsenBarostat", pressure = "100 bar", timestep = 1000}
thermostat = {type = "Berendsen", temperature = "400 K", timestep = 100}
```

Two propagators are currently implemented: one for [molecular dynamics][MD]; and
one for [Monte Carlo][MC].

[MD]: input/md.html
[MC]: input/mc.html

## Outputs

Additionally, a simulation can also output the evolution of the system
properties. Which properties are needed is given in the `outputs` array:

```toml
[[simulations]]
nsteps = 1_000_000
outputs = [
    {type = "Trajectory", file = "filename.xyz", frequency = 100},
    {type = "Energy", file = "energy.dat", frequency = 200},
    {type = "Custom", file = "custom.dat", template = "{vx[3] / mass[3]}"},
]

[simulations.propagator]
...
```

This array is an array of tables, containing three keys: the `type` of output,
the `file` to write the output to; and the `frequency` of the output. The file
is a path, and the output will be written to this path. The frequency is a
number, and the output will be written every `frequency` steps to the file.
Except for the `Trajectory` output, all files are formatted with header lines
starting with a `#`, and containing information about the quantities and the
units used for the output, and then multiple lines containing the step and the
quantities. The available outputs are the following:

- The `Energy` output will write the potential, kinetic and total energy;
- The `Cell` output will write the unit cell parameters, lengths and angles;
- The `Properties` output will write the volume, the instant pressure (computed
  from the virial equation) and the instant temperature of the system;
- The `Trajectory` output should be used to write a trajectory. The format of
  the trajectory will be guessed from the `file` extension. Supported formats
  are documented in [chemfiles](http://chemfiles.github.io/chemfiles/)
  documentation.
- The `Custom` output is the most powerful one, taking an user-provided
  template string and using it to output data. The template should be given
  as a string with the `template` key in the TOML input file.

  Here are some examples of custom output templates:
    - A constant string is reproduces as it is: `some data`;
    - Anything in braces is replaced by the corresponding values: `{pressure} {volume}`;
    - Mathematical operators are allowed in braces: `{pressure / volume}`. You
      can use `+`, `-`, `/`, `*`, `^` for exponentiation and parentheses;
    - Some properties are arrays of atomic properties `{x[0] + y[20]}`;
    - Finally, all the properties are given in the internal units. One can
      specify another unit: `x[0] / nm`.

    Here is a list of all properties available to custom outputs:

    - Atomic properties: `x`, `y` and `z` for cartesian coordinates, `vx`, `vy`
      and `vz` for cartesian components of the velocity, `mass` for the atomic
      mass, `charge` for the atomic charge.
    - Physical properties: `pressure`, `volume`, `temperature`, `natoms`
    - Unit Cell properties: `cell.a`, `cell.b`, `cell.c` are the unit cell
      vector lengths; `cell.alpha`, `cell.beta` and `cell.gamma` are the
      unit cell angles.
