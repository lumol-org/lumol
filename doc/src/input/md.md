# Molecular dynamics

A molecular dynamics simulation is started by setting the propagator `type` to
`"MolecularDynamics"`. The only needed key is the `timestep`, which is the
time step to use in the integration of forces and velocities to positions.

```toml
[[simulations]]
nsteps = 1_000_000

[simulations.propagator]
type = "MolecularDynamics"
timestep = "1 fs"
integrator = {type = "BerendsenBarostat", pressure = "100 bar", timestep = 1000}
thermostat = {type = "Berendsen", temperature = "400 K", timestep = 100}
```

Other options are the `integrator` key to use another integration scheme, the
`thermostat` key to set a thermostat, and the `controls` key to add some
additional control algorithm to the simulation.

## Integrators

Integrators are algorithms that propagate the forces acting on the particles to
compute their motions. The simplest ones performs an NVE integration, but some
integrators allow to work in different ensembles. All NVE integrators can be
turned into NVT integrators by adding a [thermostat](input/md.html#thermostats)
to the simulation. In the input, if the `integrator` key is absent, the default
integrator is a Velocity-Verlet integrator.

### Velocity-Verlet integrator

Velocity-Verlet is the most common NVE integrator for molecular dynamics. See
this [page][VelocityVerlet] for more information about the algorithm.

In the input, it can be specified by using the `VelocityVerlet` integrator type:

```toml
[simulations.propagator]
type = "MolecularDynamics"
timestep = "1 fs"
integrator = {type = "VelocityVerlet"}
```

[VelocityVerlet]: https://en.wikipedia.org/wiki/Verlet_integration#Velocity_Verlet

### Verlet integrator

Verlet algorithm is another simple NVE integrator. See this [page][Verlet] for
more information. Most of the time, the Velocity-Verlet algorithm is
preferable, since it produces more precise velocities.

In the input, it can be specified by using the `Verlet` integrator type:

```toml
[simulations.propagator]
type = "MolecularDynamics"
timestep = "1 fs"
integrator = {type = "Verlet"}
```

[Verlet]: https://en.wikipedia.org/wiki/Verlet_integration#Basic_St.C3.B6rmer.E2.80.93Verlet

### Leap-Frog integrator

The Leap-Frog algorithm is a third NVE integrator. See this [page][LeapFrog] for
details about the algorithm.

In the input, it can be specified by using the `LeapFrog` integrator type:

```toml
[simulations.propagator]
type = "MolecularDynamics"
timestep = "1 fs"
integrator = {type = "LeapFrog"}
```

[LeapFrog]: https://en.wikipedia.org/wiki/Leapfrog_integration

### Berendsen barostat

The Berendsen barostat integrator algorithm use the Berendsen barostat with a
Velocity-Verlet integrator to achieve NPT integration. It must be use together
with a thermostat, preferentially the Berendsen thermostat. See this
[page][BerendsenBarostat] for more information about the algorithm.

This algorithm exists in two versions: an isotropic one and an anisotropic one.
The isotropic version of the barostat scale all the cell parameter by the same
value using the scalar pressure. The anisotropic version scale the different
cell parameters by different values, using the stress tensor instead.

In the input, the isotropic barostat can be specified by using the
`BerendsenBarostat` integrator type:

```toml
[simulations.propagator]
type = "MolecularDynamics"
timestep = "1 fs"
integrator = {type = "BerendsenBarostat", pressure = "100 bar", timestep = 1000}
thermostat = {type = "Berendsen", temperature = "400 K", timestep = 100}
```

The `pressure` key specify the target pressure for the simulation, and the
`timestep` is the relaxation time step of the barostat.

The anisotropic barostat can be specified by using the `AnisoBerendsenBarostat`
integrator type:

```toml
[simulations.propagator]
type = "MolecularDynamics"
timestep = "1 fs"
integrator = {type = "AnisoBerendsenBarostat", pressure = "100 bar", timestep = 1000}
thermostat = {type = "Berendsen", temperature = "400 K", timestep = 100}
```

The `pressure` key specify the target hydrostatic pressure for the simulation,
and the `timestep` is the relaxation time step of the barostat.

In both cases, the barostat time step is expressed in fraction of the main
integration time step. Using a main time step of 2 fs and a barostat time step
of 1000 will yield an effective relaxation time of 2000 fs or 2 ps.

[BerendsenBarostat]: http://www.sklogwiki.org/SklogWiki/index.php/Berendsen_barostat

## Thermostats

Thermostats are algorithms used to maintain the temperature of a system at a
given value. They are specified in the input by the `thermostat` key.

### Berendsen thermostat

The Berendsen thermostat is described [here][BerendsenThermostat], and provide a
simple exponential relaxation of the temperature to a target value. In the
input, it is declared with the `Berendsen` thermostat type, a target
`temperature` value, and a `timestep`.

```toml
[simulations.propagator]
type = "MolecularDynamics"
timestep = "1 fs"
thermostat = {type = "Berendsen", temperature = "400 K", timestep = 100}
```

The time step is expressed in fraction of the main integration time step, like
for the Berendsen barostat.

[BerendsenThermostat]: http://www.sklogwiki.org/SklogWiki/index.php/Berendsen_thermostat

### Rescaling thermostat

A rescaling thermostat is the simplest thermostat algorithm possible: it just
rescale all the velocities to set the temperature to the wanted value. It can be
useful for equilibration as it converges quickly. In the input, it is specified
by the `Rescale` thermostat type, a target `temperature` value, and a
`tolerance` value. The tolerance value is optional, and is used to let the
system fluctuate around the wanted temperature: while the instant temperature is
inside the `[temperature - tolerance : temperature + tolerance]` range, no
rescale happen.

```toml
[simulations.propagator]
type = "MolecularDynamics"
timestep = "1 fs"
thermostat = {type = "Rescale", temperature = "250 K", tolerance = "10 K"}
```

## Controls

Control algorithm are supplementary steps that modify the system to ensure some
invariant, or apply some constraint. They are specified in the `controls` array,
by giving a control `type`. The `every` key specifies that the algorithm should
only be run every `n` step of the simulation (optional, defaults to 1).

```toml
[simulations.propagator]
type = "MolecularDynamics"
timestep = "1 fs"
controls = [
    # Remove global rotation of the system every 4 timestep
    {type = "RemoveRotation", every = 4}
]
```

- The `RemoveTranslation` control removes the global system rotation;
- The `RemoveRotation` control removes the global system translation.
