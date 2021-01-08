Traits used for potentials
--------------------------

As mentioned above, potential functions can be used to model all kinds of
interactions between particles such as bonded and non-bonded interactions. In
Lumol, ``Potential`` is a `trait <trait_>`_.  To further distinguish between
bonded interactions (bond lengths, angles and dihedrals) and non-bonded
interactions, we use another trait (often called marker traits). Your possible
options to further specialise a ``Potential`` are

-  `PairPotential`_ for non-bonded two body interactions;
-  `BondPotential`_ for covalent bonds interactions;
-  `AnglePotential`_ for covalent angles interactions;
-  `DihedralPotential`_ for covalent dihedral angles interactions.

.. _Potential: http://lumol.org/lumol/latest/lumol_core/energy/trait.Potential.html
.. _PairPotential: http://lumol.org/lumol/latest/lumol_core/energy/trait.PairPotential.html
.. _BondPotential: http://lumol.org/lumol/latest/lumol_core/energy/trait.BondPotential.html
.. _AnglePotential: http://lumol.org/lumol/latest/lumol_core/energy/trait.AnglePotential.html
.. _DihedralPotential: http://lumol.org/lumol/latest/lumol_core/energy/trait.DihedralPotential.html

For our Mie potential implementation, we will have to implement both the
``Potential`` as well as the ``PairPotential`` traits. If we wanted to implement
a function that can be used as non-bonded as well as bond-length potential, we'd
have to implement ``Potential``, ``PairPotential`` as well as ``BondPotential``.
"Implementing a trait" means that we will define a `struct`_ for which we will
add functions to satisfy the traits' requirements.

Let's start by having a look at the documentation for ``Potential``: open the
`API documentation  <Potential_>`_. As you can see from the trait, a
``Potential`` defines two functions, ``energy`` and ``force`` (we ignore the
``Sync + Send`` statement for now):

.. code-block:: bash

    pub trait Potential: Sync + Send {
        fn energy(&self, x: f64) -> f64;
        fn force(&self, x: f64) -> f64;
    }

``energy`` will compute the interaction energy between atoms as a function of
``x``. The force is defined as the negative derivative of the energy function
with respect to ``x``.

Both functions take a single, scalar argument and return a single scalar value.
In our case ``x`` stands for the distance between two interaction sites. Note
that only the function definitions -- without a function body -- are specified.
We will have to implement these functions for our potential.

.. _struct: https://doc.rust-lang.org/book/second-edition/ch05-00-structs.html
.. _trait: https://doc.rust-lang.org/book/second-edition/ch10-02-traits.html
