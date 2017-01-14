# Electrostatic interactions

When some particles in a system are charged, they interact with a Coulomb
potential: $$ V(x) = \frac{Z_i Z_j}{4 \pi \epsilon r_{ij}}, $$ where $Z_{i,j}$
are the net charges of the particles, $r_{ij}$ the distance between them and
$\epsilon$ the dielectric permittivity of the current medium (usually the one of
void). Because this potential goes to zero at infinity slower than $1/r^3$, it
can not be computed in periodic simulations using a cutoff distance. This
section present the available solvers for electrostatic interactions.

In the input files, electrostatic interactions are specified with two sections:
the `[charges]` section set the values of the charges of the atoms in the
system, and the `[coulomb]` section set the solver to use for the interaction.

## Charge section

Charges for the particles in the system are set in a `[charges]` section in the
potential input file. This section should contains multiple `name = <charge>`
entries, one for each charged particle in the system.

```toml
# Some salt here
[charges]
Na = 1
Cl = -1
```

## Ewald solver

Ewald idea for computing the electrostatic interactions is to split the
interaction in a short-range term which can be handled with a cutoff scheme; and
a long range term that ca be computed by using a [Fourier
transform](https://en.wikipedia.org/wiki/Fourier_transform). For more
information about the Ewald summation and its variants, see
[[Frenkel2002]](http://dx.doi.org/10.1063/1.881812).


The `[coulomb]` section for using an Ewald solver looks like this in the input
file:

```toml
[coulomb]
ewald = {cutoff = "9 A", kmax = 7}
```

The `cutoff` parameter specify the cutoff distance for the short-range and
long-range interactions splitting. The `kmax` parameter give the number of
vector to use in the reciprocal space (the long-range part of interactions).
Usually 7-8 is a good value for pure water, for very periodic charges a lower
value like 5 can do it, and for more heterogeneous system, higher values of
`kmax` are needed.

## Wolf solver

Wolf method is another method for computing electrostatic interactions presented
in [[Wolf1999]](http://dx.doi.org/10.1063/1.478738). This method replaces the
expensive computation in reciprocal space from Ewald by a corrective term, and
can be expressed as a converging sum over the charged pairs in the system.

It is accessible with the `wolf` keyword in the input files:

```toml
[coulomb]
wolf = {cutoff = "11 A"}
```

The only parameter is a `cutoff`, which should be taken bigger than the
corresponding Ewald cutoff. A value of 10 A is a good one in pure water.

---

[Frenkel2002] Frenkel, D. & Smith, B. *Understanding molecular simulation.*
(Academic press, 2002).

[Wolf1999] Wolf, D., Keblinski, P., Phillpot, S. R. & Eggebrecht, J. *Exact
method for the simulation of Coulombic systems by spherically truncated,
pairwise 1/r summation.* The Journal of Chemical Physics **110**, 8254 (1999).
