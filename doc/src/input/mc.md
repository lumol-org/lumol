# Monte-Carlo

If you want to perform a Monte-Carlo simulation, you have to set the propagator 
`type` to `"MonteCarlo"`.
Every Monte-Carlo simulations needs a `temperature` and a set of `moves` 
(actually, this set can consist of a single move).

You can think of a "move" as an instruction to change the system.
For example: "translate a single molecule in the system", "rotate a water 
molecule", or "change the cell size".
You can compose multiple moves such that you sample the desired ensemble.

Currently, all Monte-Carlo simulations are carried out using Metropolis 
acceptance criteria.

A sample input for a Monte-Carlo simulation can look like so:   

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

## Temperature

Naturally, Monte-Carlo simulations are carried out at constant 
temperature which is set using the `temperature` key. 

* `temperature`: system temperature  
    Format: String with units of temperature.

Different from Molecular Dynamics, Monte-Carlo simulations don't carry 
information about the velocities of particles.
(As a consequence we cannot access temperature from the kinetic energy.)

## Moves

All `moves` are specified as inline tables.
You can add a move using the `type` key with the name of the move. 

`moves` accept an optional `frequency` parameter.
During a Monte Carlo simulation it is very important that a move is selected 
randomly from the whole set. You can increase the chance to pick a certain move 
(compared to all other moves) by giving it a high frequency. 
If you don't specify a `frequency`, it is set to one. Frequencies 
are normalized so that they all add up to one but you don't have to do the math 
yourself: We recommend that you just use relative values and we will give you 
some examples below.



Some moves can be specified to act on a single molecule or particle type. These moves accept 
a `molecule` key whose value is a path to a file that can be read by [chemfiles][chemfiles]. 
You can add the same `type` of move multiple times.
For example you can assign different amplitudes for different species in a mixture to make sampling more efficient.
You can also use the `molecule` type to freeze a species by assigning a move 
for all other species, which can be usefull to simulation solid structures.

The molecule will then be selected with the following algorithm:

- Read the first frame of the file;
- If the file does not contain any bonding information, try to guess the bonds;
- Use the first molecule of the frame.

[chemfiles]: chemfiles.github.io

#### Example

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

The `Translate` move changes the position of a single molecule by adding 
a random displacement vector.

- Needed keys:
    * `type = "Translate"`
    * `delta`: Amplitude for displacement.  
        Format: String with units of length.
- Optional keys:
    * `frequency`: Move frequency.  
        Format: Float
    * `molecule`: Select only the specified molecule type.  
        Format: String containing path to configuration file.

If the input contains the `molecule` key, the move will only apply to one
molecule. If not, the move will apply to all molecules in the system. The
`delta` key gives the amplitude of the translations, in distance unit.

#### Example

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

The `Rotate` move randomly rotates a single molecule around its 
center-of-mass.

- Needed keys:
    * `type = "Rotate"`
    * `delta`: Amplitude for rotation.  
        Format: String with units of angle.
- Optional keys:
    * `frequency`: Move frequency.  
        Format: Float
    * `molecule`: Select only the specified molecule type.  
        Format: String containing path to configuration file.
 
If the input contains the `molecule` key, the move will only apply to one
molecule. If not, the move will apply to all molecules in the system. The
`delta` key gives the amplitude of the rotations, in units of angle.

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

The `Resize` move can be used to change the systems' volume.

- Needed keys:
    * `type = "Resize"`
    * `pressure`: Target pressure.   
        Format: String with unit of pressure.
    * `delta`: Amplitude.   
        Format: String with unit of cubic length.  
- Optional keys:
    * `frequency`: Move frequency.  
        Format: Float 

For a given `pressure`, the volume will fluctuate during the simulation such
that we sample an isobaric-isothermal ensemble.
The `delta` key sets the initial amplitude of the volume change in units of cubic length.

Note that by changing the volume, we effectively change all (center-of-mass) 
positions at once. This makes `Resize` moves expensive to use and it is 
recommended to use a comparatively low value for the `frequency`. 
As a rule of thumb, for a system containing $N$ particles, every $N + 1$'th 
move should be a `Resize` move, since a single volume change is as expensive as 
$N$ particle translations or rotations.

#### Example

For example, to simulate a system with 500 molecules, we could use the 
following input:

```toml
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
consist of 250 translations, 250 rotations and a single resizing of the 
cell *on average* (moves are picked at random with their respective frequency). 
Setting up a move set like we did above is very convenient and in literature 
you'll often find the term "cycle" (here, 1 cycle = 501 moves) to describe such a set of moves and 
frequencies. 