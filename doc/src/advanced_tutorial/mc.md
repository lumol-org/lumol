# Adding a move for the Monte-Carlo propagator

In this tutorial we will learn how to implement and add a move for the Monte-Carlo propagator. As a hands-on example, we will implement an already existing move, namely the molecular `translation` move.

A good start for creating a Monte-Carlo (MC) move is to take a look at the module structure of the MC modules. Let's first navigate to the MC folder, which you will find at `lumol/src/core/src/sim/mc`.

At this moment, the folders' structure looks like so
```
$ tree .
.
├── mod.rs
├── monte_carlo.rs
└── moves
    ├── mod.rs
    ├── rotate.rs
    └── translate.rs
```

In rust terms, a move is a struct that implements the `MCMove` trait. The trait is defined in the `moves/mod.rs` file. Let's have a look:

```rust
// File: move/mod.rs
// Comments omitted.
pub trait MCMove {
    fn describe(&self) -> &str;
    fn setup(&mut self, system: &System);
    fn prepare(&mut self, system: &mut System, rng: &mut Box<Rng>) -> bool;
    fn cost(&self, system: &System, beta: f64, cache: &mut EnergyCache) -> f64;
    fn apply(&mut self, system: &mut System);
    fn restore(&mut self, system: &mut System);
    fn update_amplitude(&mut self, scaling_factor: Option<f64>);
}
```

In Lumol, MC propagation is divided into several steps. Basically, each of these steps makes use of one of the above shown methods of `MCMove`.

- The `describe` function returns a `&str` with a brief description of the move, e.g. "translation of a molecule" or "rotation of a molecule". Primarily, this function is used for log files and errors.
- `setup` is called *once before the actual simulation is run* for every move. You can use this function to set `System` dependent parameter.
- `prepare` can be used to perform small computations to set up the moves' parameter that depend on the `System`s state.
- The `cost` function is the core functionality of a move. It is the function where most (if not all) of the computational work takes place. The return value is used from the propagator to decide if a proposed move is accepted or rejected. 
- If a move is accepted, we update the `System` by invoking the `apply` function.
- If a move is rejected, we restore the `System` to its old state by invoking the `restore` function.
- The `update_amplitude` function can be used to increase the efficiency of a simulation run by adjusting how much a move actually changes the system.

These are the functions every move has to implement. But let's get our hands dirty! In the following we reimplement an already existing move to make all the mentioned stuff more graspable.

## Implementing a molecule translation move

A simple molecular translation works as follows:
1. randomly select a molecule
2. randomly translate the center of mass position of the molecule (which leads to a translation of all positions of the atoms of said molecule)
3. compute the energy difference due to moving the molecule
4. if the move is accepted, update the positions of the atoms in the molecule
5. if the move is rejected, restore the old positions

First, lets create a new file `my_translate.rs` in the `move/` directory. 


### The struct

So, what information do we want to store in our struct that will implement the `MCMove` trait? Let us start with this

```rust
/// Monte-Carlo move for translating a molecule
pub struct Translate {
    /// Index of the molecule to translate
    molid: usize,
    /// New positions of the atom in the translated molecule
    newpos: Vec<Vector3D>,
    /// Maximum displacement value
    dr: f64,
    /// Translation range for random number generation
    range: Range<f64>,
}
```

In `molid` we will store the index of the molecule that we will translate. Instead of directly modifying the state (i.e. the `position`) of the molecule in the system, we store the positions separately in `newpos`. Note that it is a vector of `Vector3D`. Every entry stands for the cartesian coordinates (x, y, z) of an atom of the molecule. 

`dr` is called the "maximum displacement". When deciding how far a molecule is moved, we pick three random (uniformly distributed) values in the `range` [`-dr, +dr`] and create a translation vector from it - for larger `dr` values, a molecule is moved further away from its initial position *on average*. 

#### Excursion: adjusting the maximum displacement

`dr` can be adjusted during the course of a simulation to increase efficiency. To illustrate, think of a system containing densely packed hard spheres. If we pick a sphere and translate it only a tiny bit, chances are good that it will not overlap with its neighbours. That configuration will be accepted. In contrast, if `dr` would be large (and therefore the average displacement), it is very likely that the new configuration will lead to an overlap which will result in a rejection. Now, if we set `dr` to a large value and never change it, our simulation will produce rejected configurations almost exclusively.

To efficiently sample configurations it is convenient to change `dr` during the simulation based on the ratio of accepted and rejected moves. The information about a moves' statistics is stored in a `MoveCounter`. We will come back to this issue later.

**Note:**
We have to be aware of the fact that modifying the displacement will violate the reversibility (see Markov chain) of our simulation. TODO: Add details.


### `impl Translate`: Adding a `new` method

We have our struct, now let's write a function to create it. In rust, it is idiomatic to name this function `new`.

```rust
// File: my_translate.rs
impl Translate {
    /// Create a new `Translate` move with maximum displacement of `dr`.
    pub fn new(dr: f64) -> Translate {
        assert!(dr > 0.0, "dr must be a positive value in Translation move.");
        let dr = dr / f64::sqrt(3.0);
        Translate {
            molid: usize::MAX,
            newpos: Vec::new(),
            dr: dr,
            range: Range::new(-dr, dr),
        }
    }
}

impl Default for Translate {
    fn default() -> Translate {
        Translate::new(1.0)
    }
}
``` 

We also implement `Default`, which will use `dr = 1.0 Angstrom`.


### `impl MCMove for Translate`: Building the move

```rust
fn describe(&self) -> &str {
    "molecular translation"
}
```

```rust
fn setup(&mut self, _: &System) { 
    // Do nothing.
}
```

```rust
fn prepare(&mut self, system: &mut System, rng: &mut Box<Rng>) -> bool {
    if let Some(id) = select_molecule(system, self.moltype, rng) {
        self.molid = id;
    } else {
        warn!("Can not translate molecule: no molecule of this type in the system.");
        return false;
    }

    let delta = Vector3D::new(
        self.range.sample(rng),
        self.range.sample(rng),
        self.range.sample(rng)
    );

    self.newpos.clear();
    for i in system.molecule(self.molid) {
        self.newpos.push(system[i].position + delta);
    }
    return true;
}
```

```rust
```

```rust
```
