# Monte-Carlo

If you want to perform a Monte-Carlo simulation, you have to set the propagator
`type` to `"MonteCarlo"`. Every Monte-Carlo simulations needs a `temperature`
and a set of `moves` (this set can consist of a single move).

You can think of a "move" as a specific instruction to generate a new trial
configuration. For example: "translate a single molecule in the system", "rotate
a molecule", or "change the cell size". If a move is accepted (based on an
acceptance criterion), the trial configuration becomes the new configuration. If
the move is rejected, the system stays in its current configuration. You can add
multiple moves so that the ensemble of your choice is sampled.

Here is a list of all moves that are currently implemented in Lumol:
* [Translate](input/mc.html#translation): Change the center of mass position of
a molecule.
* [Rotate](input/mc.html#rotation): Perform a rotation of a molecule about its
center of mass.
* [Resize](input/mc.html#resize): Change the size of the simulation cell.

Currently, all Monte-Carlo simulations are carried out using Metropolis
acceptance criteria.

You can add all necessary information after the `[simulations.propagator]` label.

- Needed keys:
    * `type = "MonteCarlo"`
    * `temperature` (string): System temperature. The string contains the
temperature with unit.

- Optional keys:
    * `update_frequency` (positive integer): After this number of steps of a move, `delta` values
for this move are updated. Updates use statistics of a moves' acceptance ratio so it
is recommended to choose a sufficiently high number (>100).

Naturally, Monte-Carlo simulations are carried out at constant
temperature which is set using the `temperature` key.

Different from Molecular Dynamics, Monte-Carlo simulations don't carry
information about the velocities of particles.
As a consequence we cannot access temperature from the kinetic energy.

### Example
A sample input for a Monte-Carlo simulation (in the NPT ensemble)
can look like so:

```toml
[simulations.propagator]
type = "MonteCarlo"
temperature = "500 K"

moves = [
    {type = "Translate", delta = "1 A", frequency = 2},
    {type = "Rotate", delta = "20 deg", molecule = "CO2.xyz"},
    {type = "Resize", pressure = "10 bar", delta = "3 A^3", frequency = 0.001},
]
```

## Moves

All `moves` are specified as inline tables. You can add a move using the `type`
key with the name of the move.

`moves` accept an optional `frequency` parameter. During a Monte Carlo
simulation it is very important that a move is selected randomly from the whole
set. You can increase the chance to pick a certain move (compared to all other
moves) by assigning a high `frequency` to it. If you don't specify a
`frequency`, it is set to one.

Lumol normalizes frequencies after all moves are added. The easiest way to
handle frequencies is to use relative values. We will explain this below in the
given examples.

Some moves can be specified to act on a single molecule or particle type. These
moves accept a `molecule` key whose value is a path to a configuration file that
can be read by [chemfiles][chemfiles].

You can add the same move multiple times. For example, you can assign different
amplitudes for different species in a mixture to make sampling more efficient.
You can also use the `molecule` type to freeze a species by assigning a move for
all but the frozen species.

If you specify a molecule, it will be selected with the following algorithm:

- Read the first frame of the file;
- If the file does not contain any bonding information, try to guess the bonds;
- Use the first molecule of the frame.

`moves` that use a displacement (`delta`) can be added with the
`target_acceptance` key.
After a specific number of times a move was called (`update_frequency`),
we compute the acceptance ratio for the current `delta` value, *i.e.* how often
the move was accepted versus how often a move was attempted.
If the current acceptance is far away from the `target_acceptance`, we compute
a new value of `delta` based on the current acceptance.
A `target_acceptance` can only be used in conjunction with the
update_frequency` key that specifies the frequency between updates.

Sometimes a given acceptance value cannot be achieved. Either due to limits of
the adjusted `delta` value (it makes no sense to rotate a particle by more than 180° or to
translate it by multiple values of the cutoff range) or due to the nature of the
system.

To summarize, using an adjustable displacement, we can increase the efficiency of
our simulation, but strictly speaking we violate detailed balance and therefore
the Markov chain. To make sure you get correct results from your simulations, we
recommend to use adjustable displacements *only for equilibration runs*. You can
then take the resulting values for `delta` and use them for a production run,
where no further adjustments are made.

[chemfiles]: http://chemfiles.org/

### Example

```toml
# Equilibration of a protein in water.
[simulations.propagator]
type = "MonteCarlo"
temperature = "300 K"
# we update the maximum displacement `delta` after a move was called 500 times
update_frequency = 500

moves = [
    # we have much more water in the system so we want to move it more often
    # hence we set the `frequency = 100`
    # after 500 calls to this translation move, we adjust `delta` to get to approximately 50% acceptance
    {type = "Translate", delta = "2 A", molecule = "H2O.xyz", frequency = 100, target_acceptance = 0.5},
    # the single protein will be displaced only by a small distance `delta = "0.05 A"` during the whole run
    {type = "Translate", delta = "0.05 A", molecule = "protein.pdb", frequency = 1},
]
```

### Translation

The `Translate` move changes the position of a single, randomly selected  molecule
by adding a random displacement vector to its center of mass.

- Needed keys:
    * `type = "Translate"`
    * `delta` (string): Maximum amplitude for displacement.
- Optional keys:
    * `frequency` (float): Move frequency.
    * `molecule` (string): Select only the specified molecule type. The string
contains the path to the configuration file of the molecule.
    * `target_acceptance` (float): The target acceptance for this move. Value
has to be greater than zero and smaller than one. Can only be used in conjunction with `update_frequency`.

If the `molecule` key is used, the move will only apply to one molecule type. If
not, the move will apply to all molecule types in the system. The `delta` key is
the maximum magnitude of the translation vector. The conjugated string contains
the value with unit of distance.

#### Example

```toml
[simulations.propagator]
type = "MonteCarlo"
temperature = "500 K"
moves = [
    # Define a translation for all molecules in the system, including He.
    {type = "Translate", delta = "1 A", frequency = 2},
    # For He, pick a larger displacement with half the frequency of the
    # first move. Now there is a 66% chance to pick *any* molecule
    # and translate it by up to 1 A. There is a 33% chance to pick He (and only He)
    # and translate it by up to 10 A.
    {type = "Translate", delta = "10 A", molecule = "He.xyz"},
]
```

### Rotation

The `Rotate` move randomly rotates a single molecule around its
center of mass.

- Needed keys:
    * `type = "Rotate"`
    * `delta` (string): Maximum angle for rotation.
- Optional keys:
    * `frequency` (float): Move frequency.
    * `molecule` (string): Select only the specified molecule type. The string
contains the path to the configuration file of the molecule.
    * `target_acceptance` (float): The target acceptance for this move. Value
has to be greater than zero and smaller than one. Can only be used in conjunction with `update_frequency`.

If the `molecule` key is used, the move will only apply to one molecule type. If
not, the move will apply to all molecules in the system. The `delta` key is the
maximum angle. The conjugated string contains the value and the unit of either
radians or degrees (`rad` or `deg`).

#### Example

```toml
[simulations.propagator]
type = "MonteCarlo"
temperature = "500 K"
moves = [
    {type = "Rotate", delta = "3 deg", frequency = 2},
]
```

### Resize

The `Resize` move can be used to isotropically change the systems' volume.

- Needed keys:
    * `type = "Resize"`
    * `pressure` (string): Target pressure.
    * `delta` (string): Amplitude.
- Optional keys:
    * `frequency` (float): Move frequency.
    * `target_acceptance` (float): The target acceptance for this move. Value
has to be greater than zero and smaller than one. Can only be used in conjunction with `update_frequency`.

For a given `pressure`, the volume will fluctuate during the simulation. We can
use this move to sample an isobaric-isothermal ensemble. The `delta` key sets
the maximum amplitude of the volume change in units of cubic length.

By changing the volume, we effectively change all (center of mass) positions at
once. This makes `Resize` moves computationally expensive and we recommend to
use a comparatively low value for the `frequency`. As a rule of thumb, for a
system containing $N$ particles, every $N + 1$'th move should be a `Resize`
move, since a single volume change is approximately as expensive as $N$ particle
translations or rotations.

#### Example

```toml
# Simulation of 500 molecules.
[simulations.propagator]
type = "MonteCarlo"
temperature = "500 K"
moves = [
    {type = "Translate", delta = "1 A", frequency = 250},
    {type = "Rotate", delta = "20 deg", frequency = 250},
    {type = "Resize", pressure = "10 bar", delta = "3 A^3", frequency = 1},
]
```

As mentioned above, frequencies are normalized. In this example, 501 moves
consist of 250 translations, 250 rotations and a single resizing of the cell *on
average* (remember, moves are picked at random with their respective frequency).
Setting up a move set like we did in this example is very convenient and in
literature you'll often find the term "cycle" (here, 1 cycle = 501 moves) to
describe such a set of moves and respective frequencies.
