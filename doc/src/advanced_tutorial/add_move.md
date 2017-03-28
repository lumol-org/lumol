# Adding a move for the Monte-Carlo propagator

In this tutorial we will learn how to implement and add a move for the Monte-Carlo (MC) propagator.
As a hands-on example, we will reimplement an already existing move, namely the molecular `translation` move. This will be a multiple-step process:
- First, we will implement the bare functionality.
- Then, we learn how to actually use the move by adding it to the propagator.
- We come back to our implementation and talk about how we can make the move more efficient.
- Finally, we expose the move to the input module of Lumol. This will allow us to add the move using Lumol input files.

> For beginners in the field of MC simulations, we will discuss conceptual topics as they appear (in a box like this one) or as link to another page of this document. 
> If you are proficient with MC and only want to learn how to extend Lumol, feel free to skip these parts. 

To follow this tutorial, you should know about the following basic concepts of rust:
- structs
- Traits
- functions

Ready? Let's start.

## Where are source files for the MC propagator located?

Let's start by taking a look at where MC related source files are located within the Lumol project. Navigate to the MC folder, which you will find at `lumol/src/core/src/sim/mc`.

The folders' structure looks like so:
```
$ tree .
.
├── mod.rs
├── monte_carlo.rs
└── moves
    ├── mod.rs
    ├── rotate.rs
    :
    └── translate.rs
```

Initially, all of our work will take place within the `moves` directory.

In rust terms, a move is a struct that implements the `MCMove` trait. The trait is defined in the `moves/mod.rs` file. Let's have a look:

```rust
// File: move/mod.rs
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

In Lumol, MC propagation (which is defined in `monte-carlo.rs`) is divided into several steps.
Basically, each of these steps makes use of one of the above shown methods of `MCMove`.
As you can see the `System` structure is involved in most functions.
The `System` is the heart of Lumol and is a quite sophisticated structure.
At the moment it is sufficient for us to notice that the system contains the positions (`system.positions`) and that we can access the energy of the current configuration directly (`system.total_energy()`).

Going through the functions one-by-one:

- The `describe` function returns a `&str` with a brief description of the move, e.g. "translation of a molecule" or "rotation of a molecule". Primarily, this function is used for log files and errors.
- `setup` is called *once before the actual simulation is run* for every move. You can use this function to set `System` dependent parameter.
- `prepare` can be used to perform small computations to set up the moves' parameter that depend on the state of the `System`.
- The `cost` function is the core functionality of a move. It is the function where most (if not all) of the computational work takes place. The return value is used from the propagator to decide if a proposed move is accepted or rejected. 
- If a move is accepted, we update the `System` by invoking the `apply` function.
- If a move is rejected, we restore the `System` to its old state by invoking the `restore` function (if we modified it beforehand).
- The `update_amplitude` function can be used to increase the efficiency of a simulation run by adjusting how much a move actually changes the system.

These are the functions every move has to implement. But let's get our hands dirty! In the following we reimplement the molecular translation move to make all the mentioned stuff more graspable.