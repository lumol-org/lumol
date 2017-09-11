# Adding a potential function

In this tutorial we will teach you how to implement a new potential function into Lumol.
A potential is a function that describes the energy and force between interaction sites.
In Lumol we differentiate between two types of potentials: First, there are potential functions that need the global state of the system, *i.e.* all positions, as input. The Ewald summation (which is used to compute electrostatic interactions) is an example for this kind of potential. In Lumol we call them `GlobalPotential`s.
And second, we have potentials that take a single geometric parameter as input. This geometric parameter can be for example a distance or an angle. Typical examples are Van-der-Waals potentials (*e.g.* Lennard-Jones) and potential functions describing covalent bonds (*e.g.* harmonic potential, cosine potential, torsion, *etc.*). There are plenty of potentials falling in this category, hence in Lumol we simply call them `Potential`s.

> There are two kind of potentials. `GlobalPotential` take the global systems' state as input
and `Potential` that need a single scalar value (distance, angle, ...) as input.

In this tutorial we will focus on the implementation of the latter. More specific, we will implement a potential to compute van-der-Waals interactions between pairs of particles. We will make liberal use of the API documentation of both the Rust standard library as well as the Lumol API. Please note that we will point to other references, such as the Rust book, concerning general Rust concepts to keep this tutorial brief. If you have questions concerning Rust or Lumol, please don't hesitate to file an issue on github or join the discussion on our gitter.

The tutorial is structured as follows: First we will have a look at how `Potentials` are represented (what data structures are used) and what functionalities we have to implement. Then we describe where you would add those functionalities in the code. Rust offers beautiful utilities to add documentation inside the code. We will have a look at the documentation of currently implemented potentials to guide us. We will then add some tests to make sure that our implementation is correct. This concludes the bulk the work, but to make our new potential function accessible to all Lumol users we will also add a parsing function. Doing so, our potential can be conveniently specified from an input file. We conclude the tutorial by adding a short documentation to the user manual.

We have a lot to do. Ready? Let's go.

## Rust prerequisites

We recommend reading these chapters in the rust book.

- [structs](https://doc.rust-lang.org/book/second-edition/ch05-00-structs.html)
- [traits](https://doc.rust-lang.org/book/second-edition/ch10-02-traits.html)

## Representation

As mentioned above, potential functions can be used to model all kinds of interactions between particles such as bonded  and non-bonded interactions. In Lumol, `Potential` is a [trait](https://doc.rust-lang.org/book/second-edition/ch10-02-traits.html). To further distinguish between bonded interactions (bond lengths, angles and dihedrals) and non-bonded interactions, we use another trait (often called marker traits). Your possible options to further specialise a `Potential` are

- [`PairPotential`](http://lumol.org/lumol/latest/lumol/energy/trait.PairPotential.html) for non-bonded two body interactions;
- [`BondPotential`](http://lumol.org/lumol/latest/lumol/energy/trait.BondPotential.html) for covalent bonds interactions;
- [`AnglePotential`](http://lumol.org/lumol/latest/lumol/energy/trait.AnglePotential.html) for covalent angles interactions;
- [`DihedralPotential`](http://lumol.org/lumol/latest/lumol/energy/trait.DihedralPotential.html) for covalent dihedral angles interactions.

In this tutorial, we will implement a potential to describe non-bonded pair interactions, namely the [**Mie potential**](http://www.sklogwiki.org/SklogWiki/index.php/Mie_potential). We will have to implement both the `Potential` as well as the `PairPotential` traits. If we wanted to implement a function that can be used as non-bonded as well as bond-length potential, we'd have to implement `Potential`, `PairPotential` as well as `BondPotential`. "Implementing a trait" means that we will define a [structure](https://doc.rust-lang.org/book/second-edition/ch05-00-structs.html) (a `struct`) for which we will add functions to satisfy the traits' requirements.

Let's start by having a look at the API documentation for `Potential`. Open the [API documentation](http://lumol.org/lumol/latest/lumol/index.html) and navigate to the `energy` module which is listed under the `modules` section at the bottom of the page. Traits are listed at the very bottom of the module documentation. Open the documentation for the `Potential` trait.

As you can see from the trait, a `Potential` defines the interface for two functions, `energy` and `force` (we ignore the `Sync + Send` statement for now):

```rust
pub trait Potential: Sync + Send {
    fn energy(&self, x: f64) -> f64;
    fn force(&self, x: f64) -> f64;
}
```

It's pretty self-explanatory what these to functions do: `energy` will compute the interaction energy between two sites as a function of `x`. The force is defined as the negative derivative of the energy function with respect to `x`.

Both functions take a single, scalar argument and return a single scalar value. In our case `x` stands for the distance between two interaction sites. Note that only the function definitions -- without a function body -- are specified. We will have to implement these functions for our potential.

Let's start the implementation.

## Adding functionality and documentation into the code

### Defining the struct

The energy function of the Mie potential reads

$u(x) = \frac{n}{n-m} \left(\frac{n}{m}\right)^{m/(n-m)}\epsilon \left[ \left( \frac{\sigma}{x}\right)^n - \left( \frac{\sigma}{x}\right)^m \right]$

where $x$ denotes the distance between two interaction sites $i, j$, $x = x_{ij} = | \mathbf{r}_j - \mathbf{r}_i |$. The parameters of the potential are

- $n, m$ the repulsive and attractive exponents, respectively,
- $\epsilon$ the energetic paramater,
- $\sigma$ the particle diameter or structural parameter.

We start by defining the `struct` for our potential. We will add our code to the `functions.rs` file located in `lumol/src/core/src/energy`. Add the following lines:

```rust
#[derive(Clone, Copy)]
pub struct Mie {
    /// Distance constant
    pub sigma: f64,
    /// Energy constant
    pub epsilon: f64,
    /// Exponent of repulsive contribution
    pub n: f64,
    /// Exponent of attractive contribution
    pub m: f64,
    /// Energetic prefactor computed from the exponents and epsilon
    pub prefac: f64,
}
```

Next, we implement a constructor function. This is usefull in this case since we want to compute the prefactor
of the potential once before we start our simulation.
In Rust we typically use `new` for the constructors' name.

```rust
impl Mie {
    fn new(sigma: f64, epsilon: f64, n: f64, m: f64) -> Mie {
        if m <= n {
            panic!("The repulsive exponent n has to be larger than the attractive exponent m")
        };
        let prefac = n/(n-m) * (n/m).powf(m/(n-m)) * epsilon;
        Mie {
            sigma: sigma,
            epsilon: epsilon,
            n: n,
            m: m,
            prefac: prefac
        }
    }
}
```

Our function takes the parameter set as input, computes the prefactor and returns a `Mie` struct. Notice that
it panics, `n` is smaller or equal than `m`.
The next step is to implement the `Potential` trait for `Mie`.

### Implementing `Potential`

Add the following lines below the structs implementation.

```rust
impl Potential for Mie {
    fn energy(&self, r: f64) -> f64 {
        let sigma_r = self.sigma/r;
        let repulsive = f64::powf(sigma_r, self.n);
        let attractive = f64::powf(sigma_r, self.m);
        self.prefac * (repulsive - attractive)
    }

    fn force(&self, r: f64) -> f64 {
        let sigma_r = self.sigma/r;
        let repulsive = f64::powf(sigma_r, self.n);
        let attractive = f64::powf(sigma_r, self.m);
        - self.prefac * (self.n*repulsive - self.m*attactive) / r
    }
}
```

`energy` is an implementation of the Mie potential equation shown above.
`force` is the negative derivative of the energy with respect to the distance, `r`.
To be more precise, the vectorial force can readily computed by multiplying the
result of `force` with the connection vector $\vec{r}$.

The next step is to make our `Potential` usable in Lumol's algortihms to compute
non-bonded energies and forces.
Therefore, we have to implement the `PairPotential` trait.

### Implementing `PairPotential`

Let's inspect the [documentation of `PairPotential`](http://lumol.org/lumol/latest/lumol/energy/trait.PairPotential.html).

```rust
pub trait PairPotential: Potential + BoxClonePair {
    fn tail_energy(&self, cutoff: f64) -> f64;
    fn tail_virial(&self, cutoff: f64) -> f64;

    fn virial(&self, r: &Vector3D) -> Matrix3 { ... }
}
```

First, we can see that `PairPotential` enforces the implementation of `Potential` which is denoted by `pub trait PairPotential: Potential ...` (we ignore `BoxClonePair` for now).
That makes sense from a didactive point of view since we said that `PairPotential` is a "specialization" of `Potential`
and furthermore, we can make use of all functions that we had to implement for `Potential`.

There are three functions in the `PairPotential` trait. The first two functions start with `tail_`.
These are functions to compute long range or tail corrections.
Often, we introduce a cutoff distance into our potential beyond which we set the energy to zero.
Doing so we intoduce an error which we can account for using a tail correction.
We need two of these corrections, one for the energy, `tail_energy`, and one for the pressure (which uses `tail_virial`
under the hood). For a beautiful derivation of tail corrections for cut potentials, [see here](https://engineering.ucsb.edu/~shell/che210d/Simulations_of_bulk_phases.pdf).

The third function, `virial`, already has its body implemented -- we don't have to write an implementation for
our potential.

We will omit the derivation of the formulae for tail corrections here but they are computed by solving these equations

energy: $\int_{r_c}^{\infty} u(r) r^2 \mathrm{d}r$

virial: $\int_{r_c}^{\infty} \frac{\partial u(r)}{\partial r} r^3 \mathrm{d}r$

The implementation looks like so

```rust
impl PairPotential for Mie {
    fn tail_energy(&self, cutoff: f64) -> f64 {
        if self.m <= 3.0 {
            panic!("No tail correction possible for Mie potential with an exponent lower than 3.0")
        };
        let sigma_rc = self.sigma/cutoff;
        let n_3 = self.n - 3.0;
        let m_3 = self.m - 3.0;
        let repulsive = f64::powf(sigma_rc, n_3);
        let attractive = f64::powf(sigma_rc, m_3);
        - self.prefac * self.sigma.powi(3) * (repulsive/n_3 - attactive/m_3)
    }

    fn tail_virial(&self, cutoff: f64) -> f64 {
        if self.m <= 3.0 {
            panic!("No tail correction possible for Mie potential with an exponent lower than 3.0")
        };
        let sigma_rc = self.sigma/cutoff;
        let n_3 = self.n - 3.0;
        let m_3 = self.m - 3.0;
        let repulsive = f64::powf(sigma_rc, n_3);
        let attractive = f64::powf(sigma_rc, m_3);
        - self.prefac * self.sigma.powi(3) * (repulsive * self.n / n_3 - attactive * self.m /m_3)
    }
}
```

Note that we cannot correct every kind of energy function.
In fact, the potential has to be a *short ranged* potential.
For our Mie potential, both the exponents have to be larger than 3.0 else our potential will be *long ranged* and
the integral that has to be solved to compute the tail corrections diverges.

### Adding API documentation

```rust
/// Mie potential.
///
/// The following expression of the Mie potential is used: `V(r) = n/(n-m)
/// * (n/m)^(m/(n-m)) * epsilon * ((sigma/r)^n - (sigma/r)^m)` where `sigma` is the Mie
/// distance constant, `epsilon` the Mie energetic constant and `n`, `m` are the
/// repulsive and attractive exponents, respectively.
///
/// # Examples
///
/// ```
/// use lumol::energy::Potential;
/// use lumol::energy::Mie;
///
/// let potential = Mie{sigma: 2.0, epsilon: 10.0, n: 12.0, m: 6.0, prefac: 4.0};
/// assert_eq!(potential.energy(2.0), 0.0);
/// assert_eq!(potential.energy(3.0), -3.203365942785746);
///
/// assert_eq!(potential.force(2.0), 120.0);
/// ```
#[derive(Clone, Copy)]
pub struct Mie {
    /// Distance constant
    pub sigma: f64,
    /// Energy constant
    // rest of struct is omitted
```

```rust
/// Returns a new Mie potential.
impl Mie {
    fn new(sigma: f64, epsilon: f64, n: f64, m: f64) -> Mie {
    // implementation omitted
```


## Adding tests


## Exposing the potential to the input module

## Writing documentation for the user manual


