# Monte-Carlo

A Monte-Carlo simulation is started by setting the propagator `type` to
`MonteCarlo`. The two needed keys are the `temperature` of the simulation, and
the Monte-Carlo `moves` to use.

```toml
[simulations.propagator]
type = "MonteCarlo"
temperature = "500 K"
moves = [
    {type = "Translate", delta = "1 A", frequency = 2},
    {type = "Rotate", delta = "20 deg", molecule = "CO2.xyz"},
]
```

## Moves

All the Monte-Carlo moves are specified as inline tables, with the `type` key
setting the type of move. They all accept an optional `frequency` value setting
at which frequency the move should be selected. By default the frequency is 1,
and the actual frequency will be the `frequency` value divided by the sum of all
the frequencies.

Some moves acting a single molecule also accept a `molecule` key, giving the
path to a file which [chemfiles][chemfiles] can read. The molecule associated
with the move will then be selected with the following algorithm:

- Read the first frame of the file;
- If the file does not contain any bonding information, try to guess the bonds;
- Use the first molecule of the frame.

[chemfiles]: chemfiles.github.io

The same move type can be present more than once in the input, to have different
amplitudes for different compounds for example.

```toml
[simulations.propagator]
type = "MonteCarlo"
temperature = "500 K"
moves = [
    {type = "Translate", delta = "2 A", molecule = "H2O.xyz"},
    {type = "Translate", delta = "0.5 A", molecule = "protein.pdb"},
]
```


### Translation

The `Translate` move type will do small translations of a single molecule at the
time. If the input contains the `molecule` key, the move will only apply to one
molecule. If not, the move will apply to all molecules in the system. The
`delta` key gives the amplitude of the translations, in distance unit.

```toml
[simulations.propagator]
type = "MonteCarlo"
temperature = "500 K"
moves = [
    {type = "Translate", delta = "1 A", frequency = 2},
    {type = "Translate", delta = "10 A", molecule = "He.xyz"},
]
```

### Rotation

The `Rotate` move type will do small rotations of a single molecule at the
time. If the input contains the `molecule` key, the move will only apply to one
molecule. If not, the move will apply to all molecules in the system. The
`delta` key gives the amplitude of the rotations, in angle unit.

```toml
[simulations.propagator]
type = "MonteCarlo"
temperature = "500 K"
moves = [
    {type = "Rotate", delta = "3 deg", frequency = 2},
]
```
