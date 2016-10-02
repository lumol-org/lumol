# Monte-Carlo simulation of Argon

So let's run a simulation with Lumol. The easiest system to simulate is a
Lennard-Jones fluid, which is a good model for noble gases fluids. Here we will
simulate supercritical argon using Metropolis Monte-Carlo algorithm.

For this simulation, you will need the initial configuration available
[here][argon.xyz]. Download it and save it under the name `argon.xyz`.

You will also need the input file available [here][argon.toml]. Save it as
`argon.toml` and place it near to `argon.xyz`. You can then run the simulation
with:

```
lumol argon.toml
```

The simulation should complete in a few second (8s on my system), and produce
two files: `energy.dat` and `trajectory.xyz`.

## Input file anatomy

The input file is written using the TOML syntax, you can learn more about this
[syntax here][TOML]. The file starts with a header declaring the version of the
input file syntax used, here the version 1:

```toml
[input]
version = 1
```

Then, we declare which system we want to simulate, in the `systems` array. We
define this system using an XYZ file, and providing the unit cell size.

```toml
[[systems]]
file = "argon.xyz"
cell = 21.65
```

We also need to define the interactions potential between the atom in the
system, which we do in the `potential.pairs` section; using a Lennard-Jones
potentials with a cutoff distance of 10 A for all Ar-Ar pairs.

```toml
[[systems.potentials.pairs]]
atoms = ["Ar", "Ar"]
cutoff = "10 A"
lj = {sigma = "3.4 A", epsilon = "1.0 kJ/mol"}
```

Then we define how we want to simulate our system. Here we need to run the
simulation for `100000` steps, and output the energy to `energy.dat` every 100
steps, and the trajectory to `trajectory.xyz` every 100 steps too.

```toml
[[simulations]]
nsteps = 100000
outputs = [
    {type = "Energy", file = "energy.dat", frequency = 100},
    {type = "Trajectory", file = "trajectory.xyz", frequency = 100}
]
```

At the end we define how we propagate the system from one step to another. Here
we are using a Monte-Carlo simulation at 500 K, and the only Monte-Carlo move is
a translation of maximum amplitude of 1 A.

```toml
[simulations.propagator]
type = "MonteCarlo"
temperature = "500 K"
moves = [
    {type = "Translate", delta = "1 A"},
]
```

So now we know how to run simulations of Lennard-Jones fluids. How about we add
electrostatic interactions in the next example? Say hello to [Sodium
Chloride][nacl]!

[argon.xyz]: data/argon.xyz
[argon.toml]: data/argon.toml
[TOML]: https://github.com/toml-lang/toml
[nacl]: tutorial/nacl.html
