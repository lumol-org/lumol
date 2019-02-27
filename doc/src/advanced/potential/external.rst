Implementation of the potential
===============================

We start by creating a new package using ``cargo``:

.. code-block:: bash

    cargo new potential
    cd potential

Open ``Cargo.toml`` and add the lines

.. code-block:: bash

    [dependencies]
    lumol = {git = "https://github.com/lumol-org/lumol"}

to add the ``lumol`` crate as a dependency to the package. To test if everything
works, run ``cargo build`` and check if an error occurs.

Defining the struct
-------------------

For the first part of the tutorial, the complete code will be written into the
``lib.rs`` file.

The energy function of the Mie potential reads

.. math::

    u(x) = \varepsilon \frac{n}{n-m} \left(\frac{n}{m}\right)^{\frac{m}{n-m}} \left[ \left( \frac{\sigma}{x}\right)^n - \left( \frac{\sigma}{x}\right)^m \right]

where :math:`x` denotes the distance between two interaction sites :math:`i, j`,
with :math:`x = x_{ij} = | \mathbf{r}_j - \mathbf{r}_i |`.  The parameters of
the potential are

-  :math:`n, m` the repulsive and attractive exponents, respectively,
-  :math:`\varepsilon` the energetic paramater,
-  :math:`\sigma` the particle diameter or structural parameter.

We start by defining the ``struct`` for our potential. Add the following lines
to ``lib.rs``:

.. literalinclude:: ../../../../tutorials/potential/src/lib.rs
    :language: rust
    :lines: 1-13

In the first two lines we define our imports from ``Lumol``, following with our
``Mie`` structure. Notice that we don't store the ``epsilon`` value, instead we
store an energetic prefactor that will make it easier to compute the potential.

.. math::

    \text{prefactor} = \varepsilon \frac{n}{n-m} \left(\frac{n}{m}\right)^{m/(n-m)}

Next, we implement a constructor function. That's usefull in this case since we
want to compute the prefactor of the potential once before we start our
simulation.

In Rust we typically use ``new`` for the constructors' name.

.. literalinclude:: ../../../../tutorials/potential/src/lib.rs
    :language: rust
    :lines: 15-28

Our function takes the parameter set as input, computes the prefactor and
returns a ``Mie`` struct. Notice that it panics, for ``n`` smaller than or equal
to ``m``. The next step is to implement the ``Potential`` trait for ``Mie``.

Implementing ``Potential``
--------------------------

Add the following lines below the structs implementation.

.. literalinclude:: ../../../../tutorials/potential/src/lib.rs
    :language: rust
    :lines: 30-44

``energy`` is the implementation of the Mie potential equation shown above.
``force`` is the negative derivative of the energy with respect to the distance,
``r``. To be more precise, the vectorial force can readily be computed by
multiplying the result of ``force`` with the connection vector :math:`\vec{r}`.

The next step is to make our ``Potential`` usable in Lumol's algorithms to
compute non-bonded energies and forces. Therefore, we have to implement the
``PairPotential`` trait.

Implementing ``PairPotential``
------------------------------

Let's inspect the `documentation <PairPotential_>`_  for ``PairPotential``.

.. _PairPotential: http://lumol.org/lumol/latest/lumol/energy/trait.PairPotential.html

.. code-block:: rust

    pub trait PairPotential: Potential + BoxClonePair {
        fn tail_energy(&self, cutoff: f64) -> f64;
        fn tail_virial(&self, cutoff: f64) -> f64;

        fn virial(&self, r: &Vector3D) -> Matrix3 { ... }
    }

First, we can see that ``PairPotential`` enforces the implementation of
``Potential`` which is denoted by ``pub trait PairPotential: Potential ...`` (we
ignore ``BoxClonePair`` for now, as it is automatically implemented for us if we
implement ``PairPotential`` manually). That makes sense from a didactive point
of view since we said that ``PairPotential`` is a "specialization" of
``Potential`` and furthermore, we can make use of all functions that we had to
implement for ``Potential``.

There are three functions in the ``PairPotential`` trait. The first two
functions start with ``tail_``. These are functions to compute long range or
tail corrections. Often, we introduce a cutoff distance into our potential
beyond which we set the energy to zero. Doing so we intoduce an error which we
can account for using a tail correction. We need two of these corrections, one
for the energy, ``tail_energy``, and one for the pressure (which uses
``tail_virial`` under the hood). For a beautiful derivation of tail corrections
for truncated potentials, `see here <tail-correction_>`_.

.. _tail-correction: https://engineering.ucsb.edu/~shell/che210d/Simulations_of_bulk_phases.pdf

The third function, ``virial``, already has its body implemented -- we don't
have to write an implementation for our potential.

We will omit the derivation of the formulae for tail corrections here but they
are computed by solving these equations

.. math::
    \text{tail energy} = \int_{r_c}^{\infty} u(r) r^2 \mathrm{d}r

.. math::
    \text{tail virial} = \int_{r_c}^{\infty} \frac{\partial u(r)}{\partial r} r^3 \mathrm{d}r

The implementation looks like so

.. literalinclude:: ../../../../tutorials/potential/src/lib.rs
    :language: rust
    :lines: 46-70

Note that we cannot correct every kind of energy function. In fact, the
potential has to be a *short ranged* potential. For our Mie potential, both the
exponents have to be larger than 3.0 else our potential will be *long ranged*
and the integral that has to be solved to compute the tail corrections diverges.
We return zero in that case.

Running a simulation
--------------------

That concludes the first part. To test your new and shiny potential, you can run
a small simulation. You'll find a minimal Monte Carlo simulation example in the
``tutorials/potential`` directory of the main lumol repository
where you will also find the ``src/lib.rs`` file we created in this tutorial.
You can then run the simulation via

.. code-block:: bash

    cargo run --release

Fantastic! You implemented a new potential and ran a simulation with it!

If you want to share your implementation with other Lumol users only some small
additional steps are neccessary. We will talk about them in the next part of
this tutorial (which is not yet written).
