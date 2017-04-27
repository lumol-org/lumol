# A minimalist implementation

> What is a MC move again?  
> In MC, moves are used to propagate the system.
> That is, a move is a specific instruction about how to change the current configuration and to propose a new configuration.
> Configuration here means the positions of all atoms (and not the velocities - in contrast to molecular dynamics).
> A move can be accepted - then the system will adopt the proposed configuration - or rejected - then the system will not change.
> Moves can never happen simultaneously - they will take place always one after another

A simple translation of a molecule works as follows:
1. randomly select a molecule of the system
2. randomly translate the center of mass position of the molecule (which leads to a translation of all positions of the atoms of said molecule)
3. compute the energy difference due to moving the molecule. The energy difference is needed to decide wheter a move is accepted or rejected.
4. if the move is accepted, update the positions of the atoms in the molecule
5. if the move is rejected, leave the positions untouched

Because `McMove` is a trait, we have to generate a `struct` for our move that implements the `MCMove` trait. We first consider some fields in which we store the  information needed. Let's start with this:

```rust
/// Monte-Carlo move for translating a molecule
pub struct Translate {
    /// System state prior to the move
    old_system: System,
    /// Maximum displacement value
    dr: f64,
    /// Translation range for random number generation
    range: Range<f64>,
}
```

We store the current (old) system before changing it in `old_system`.
`dr` is called the "maximum displacement". When deciding how far we move a molecule, we pick three random (uniformly distributed) values in the `range` [`-dr, +dr`] and create a translation vector from it.
For larger `dr` values, a molecule moves further away from its initial position *on average*. 

> About the maximum displacement:  
> `dr` can be adjusted during the course of a simulation to increase efficiency.
> To illustrate this, think of a system containing densely packed hard spheres.
> Hard spheres are spheres that do not interact with each other but are not allowed to overlap.
> If we pick a sphere and translate it only a tiny bit, chances are good that it will not overlap with its neighbours.
> Accordingly, this new configuration is accepted.
> In contrast, if `dr` is large (and therefore the average displacement), it is very likely that the new configuration leads to an overlap and the proposed configuration is rejected.
> Now, if we set `dr` to a large value and never change it, our move produces
configurations that are rejected almost exclusively.
> To efficiently sample configurations it is convenient to change `dr` during the simulation based on the ratio of accepted and rejected moves.
> We will come back to this issue in our refined implementation.

We have our struct, now let's write a function to create an instance of it. In rust, it is idiomatic to name this function `new`.

```rust
// File: translate.rs
impl Translate {
    fn new(dr: f64) -> Translate {
        Translate {
            old_system: System::default(),
            dr: dr,
            range: Range::new(-dr, dr),
        }
    }
}
```

The `new` function takes a single argument, the maximum displacement `dr`.
The internal unit of length - the displacement is a length - is Angstrom.
Remember that in rust, we need to initialize all fields of our struct, so we use the default implementation for our system.

```rust
impl Default for Translate {
    fn default() -> Translate {
        Translate::new(1.0)
    }
}
``` 

We also implement `Default`, which simply calls `new` with a displacement of `dr = 1.0 angstrom`. 

## Implement the `MCMove` trait.

Good. We are able to create an instance of the struct.
Let's start implementing the `MCMove` trait.
All functions that follow are within the `impl` brackets.

```rust
impl MCMove for Translate {
    // here will be our functions.
}
```

First, we implement the `describe` function:

```rust
fn describe(&self) -> &str {
    "molecular translation"
}
```

Pretty self explanatory - a simple description of what the move does.

Next, the `setup`.
As mentioned earlier, this function is called once before the simulations starts.
Since it takes the system as parameter, we can compute quantities that may be important for the move - depending on the system we use.
For our initial implementation though, we don't use the `setup` function.

```rust
fn setup(&mut self, _: &System) { 
    // For now, do nothing.
}
```
The `prepare` function is a bit more complicated:

```rust
fn prepare(&mut self, system: &mut System, rng: &mut Box<Rng>) -> bool {
    // select a molecule randomly: pick a random number between 0 and the total number
    // of molecules in the system.
    let nmols = system.molecules().len();
    if nmols == 0 {
        return false // no molecule, so we can't perform the move
    } else {
        let molid = rng.gen_range(0, nmols)
    }

    // create a vector, containing the random displacements in every cartesian
    // direction.
    let delta = Vector3D::new(
        self.range.sample(rng), // x
        self.range.sample(rng), // y
        self.range.sample(rng), // z
    )

    self.old_system = system.clone(); // before modifying the system, create a copy
    // use the current positions and add the displacement
    for i in self.old_system.molecule(molid) {
        system[i].position += delta
    }
    true // we set everything up, so we can perform the move.
}
```

That is a lot of stuff.
Let's first have a look at the function arguments:
We take three mutable references: One to the `Translate` struct (`&mut self`), one
to the system and finally one to `Box<Rng>`, which is a **r**andom **n**umber **g**enerator.

Let us now go through the code blocks one by one.

First, we randomly select a molecule.
All molecules are stored in a vector inside system, `system.molecules`.
We cannot access this vector directly, since it is not public, but there is a public function that we can use: `system.molecules()` will return the vector of molecules.
The `len()` function will return the number of elements in the vector which equals the number of molecules.
To pick a random molecule, we draw a random number between 0 and the number of molecules.
Note that we also check if there are any molecules in the system.

Next, we create the random displacement.
In Lumol, positions are stored as `Vector3D` structures.
It is a vector with three elements, each of which stands for a cartesian coordinate.
We create a random displacement vector by drawing a random number for every cartesian coordinate.
In rust, we can do that by taking a sample from `range` which we set to be `[-dr, dr]` in our `new()` function.

To apply the displacement, we first store the current state of the system using the `clone` function.
In the `for` loop, we add the displacement by modifying the positions of the system.
To understand what's going on here, we have to take a closer look on how positions are stored in Lumol and how we can access them:
```rust
for i in self.old_system.molecule(molid) {
// ...
}
```
The `molecule(molid)` function is a function of `System`.
It takes an index (a number of type `usize`) as argument and returns a vector containing the indices of the particles that form the molecule.
Why is that?
Well, in Lumol, a molecule stores information about the bonds between particles but not the particles (or their positions) themselfs.
We want to add the displacement to all positions of a particle and to do that we first ask for the indices of the randomly selected molecule.
Using the `for i in ...` syntax, we now have a loop over the indices, that are stored (referenced to) in the binding `i`.
Finally, we can access the position of a particle with index  `i` using `system[i].position`.

The cost function is comparably short but don't get fooled as it is the most expensive function of the move.

```rust
fn cost(&self, system: &System, beta: f64, _: &mut EnergyCache) -> f64 {
    let new_energy = self.old_system.total_energy();
    let old_energy = system.total_energy();
    let delta_energy = new_energy - old_energy;
    beta * delta_energy // return the cost function
}
```

All the work is hidden in the `total_energy()` call that will compute the system's total energy *from scratch*.
The cost function in our move equals the dimensionless change of the potential energy - that is, the difference of the energy before and after the move.
The total energy of a system contains the potential as well as the kinetic energy, but since we are using Monte-Carlo, the kinetic energy is zero.
Beta is the inverse temperature.
Tt is defined as:

$$ \beta = \frac{1}{k_B T}$$

and its unit is the inversed energy. 

Since we apply the change to the system and store the old system in our struct, applying the move means just going to the next move using the new positions.

```rust
fn apply(&mut self, system: &mut System) {
    // We already changed the system. Nothing to do here.
}
```

On the other hand, when the move is rejected, we restore the old system by simply swapping it back.

```rust
fn restore(&mut self, system: &mut System) {
    // Restore the old system: swap the memory
    mem::swap(system, &mut self.old_system)
}
```

As mentioned before, changing `dr` is not in scope of this part of the tutorial.
Check out the refined version to learn how we can use updates of `dr` to increase efficiency.

```rust
fn update_amplitude(&mut self, scaling_factor: Option<f64>) {
    // For now, we don't add a rule to change `dr`
}
```

We did it!
That is all we have to do to implement our move.

## What we learned

In summary, we learned that to implement a move, we have to create a `struct` for which we have to implement the `MCMove` trait.
The `MCMove` trait contains a number of functions that are used by the MC propagator.
We talked a bit about what the different functions are and implemented them for a molecular translation.
Doing so, we learned that the central structure of the simulation is the `System` that holds information about molecules and particles and also algorithms to compute the energy.
To access molecules, positions and energies, we used functions provided by the system.

## Ready to move on ...

let's actually use our move in a simulation!


