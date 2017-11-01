Adding a potential
==================

In this tutorial we will teach you how to use your own potential function with
Lumol. This tutorial consists of two parts. In the first part, we will implement
the logic in a separate project using Lumol as an external library. The second
part describes the neccessary steps to implement the logic into the core package
of Lumol called ``lumol-core``.

A potential is a function that describes the energy and force between
interaction sites. In Lumol we differentiate between two types of potentials:
First, there are potential functions that need the global state of the system,
*i.e.* all positions, as input. The Ewald summation which is used to compute
electrostatic interactions is an example for this kind of potential. We call
them ``GlobalPotential``\ s. And second, we have potentials that take a single
geometric parameter as input. This geometric parameter can be for example a
distance or an angle. Typical examples are Van-der-Waals potentials (*e.g.*
Lennard-Jones) and potential functions describing covalent bonds (*e.g.*
harmonic potential, cosine potential, torsion, *etc.*). There are plenty of
potentials falling into this category, hence in Lumol we simply call them
``Potential``\ s.

.. note::

    There are two kind of potentials. A ``GlobalPotential`` takes the global
    systems' state as input and a ``Potential`` takes a single scalar value
    (distance, angle, ...) as input.

In this tutorial we will focus on the implementation of the latter. More
specific, we will implement a potential to compute van-der-Waals interactions
between pairs of particles. We will make liberal use of the API documentation of
both the Rust standard library as well as the Lumol API. Please note that we
will point to other references, such as the Rust book, concerning general Rust
concepts to keep this tutorial brief.  If you have questions concerning Rust or
Lumol, please don't hesitate to file an issue on `github`_ or join the
discussion on our `gitter`_.

.. _github: https://github.com/lumol-org/lumol/
.. _gitter: https://gitter.im/lumol-org/lumol

The tutorial is structured as follows: First we will have a look at how
``Potentials`` are represented (what data structures are used) and what
functionalities we have to implement. Then, we describe how we add those
functionalities. We will write a small simulation program that makes use of our
newly created potential. In the second part we talk about implementing the
potential into Lumol's core. Rust offers beautiful utilities to add
documentation inside the code. We will have a look at the documentation of
currently implemented potentials to guide us. We will then add some tests to
make sure that our implementation is correct. This concludes the bulk the work,
but to make our new potential function accessible to all Lumol users we will
also add a parsing function. Doing so, our potential can be conveniently
specified from an input file. We conclude the tutorial by adding a short
documentation to the user manual.

In this tutorial, we will implement a potential to describe non-bonded pair
interactions, namely the `Mie potential`_.

We have a lot to do. Ready? Let's go.

.. note::

    We recommend you read the chapters concerning `structs`_ and `traits`_ in
    the rust book.

.. toctree::

    intro
    external

.. _structs: https://doc.rust-lang.org/book/second-edition/ch05-00-structs.html
.. _traits: https://doc.rust-lang.org/book/second-edition/ch10-02-traits.html
.. _Mie potential: http://www.sklogwiki.org/SklogWiki/index.php/Mie_potential
