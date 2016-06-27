# Available potentials

This section is a list of all the available potentials in Cymbalum, with the
associated parameters. All potentials have to provide additional parameters in
there definition, as a TOML table. Using inline tables is the easiest way to do
so:

```toml
[[pairs]]
atoms = ["A", "B"]
# Additional parameters here are 'sigma' and 'epsilon'.
lj = {sigma = "3 A", epsilon = "123 kJ/mol"}
```

The same potential can be used for either pairs (at distance $r$); or for
angles (at angle $\phi$). In all the formulas, the $x$ parameter
represents either a distance or an angle.

## Null potential

This potential is always 0, for all values of $x$. It should be used to remove
interactions between atoms in a pair/bond/angle/dihedral that are
present in the system but should not be interacting.

This potential can be used by specifying the `null` key with an empty table `{}`
as value.

```toml
[[pairs]]
atoms = ["O", "O"]
null = {}
```

## Lennard-Jones potential

The Lennard-Jones potential is a classical potential for pair interactions
expressed as: $$ V(x) = 4 \epsilon \left[\left(\frac{\sigma}{x}\right)^{12} -
\left(\frac{\sigma}{x}\right)^6\right].$$

The Lennard-Jones potential is defined using the `lj` or `lennardjones` key. The
parameters are `sigma` ($\sigma$) and `epsilon` ($\epsilon$), which should be
provided as strings.

```toml
[[pairs]]
atoms = ["O", "O"]
lennardjones = {sigma = "3.16 A", epsilon = "0.155 kcal/mol"}
```

## Harmonic potential

The Harmonic potential is usually used for intramolecular interactions such as
bonds, angles or dihedrals. It is expressed as:
$$ V(x) = \frac 12 k \ (x - x_0)^2$$

The potential type keyword is `harmonic`, and the parameters are `k` and `x0`,
provided as strings.

```toml
[[bonds]]
atoms = ["O", "H"]
harmonic = {k = "1054.2 kcal/mol/A^2", x0 = "1.0 A"}

[[angles]]
atoms = ["H", "O", "H"]
harmonic = {k = "75.9 kcal/mol/rad^2", x0 = "109.5 deg"}
```

## Cosine-Harmonic potential

This potential is usually used for angles and dihedral angles interactions,
because it presents a $2\pi$ periodicity. It is expressed as: $$ V(x) = \frac 12
k \ (\cos x - \cos x_0)^2$$

The potential type keyword is `cosine-harmonic`, and the parameters `k` and `x0`
should be provided as strings.

```toml
[[angles]]
atoms = ["H", "C", "H"]
cosine-harmonic = {k = "67 kJ/mol", x0 = "120 deg"}
```

## Torsion potential

This potential is usually used for dihedral interactions. It is
expressed as: $$ V(x) = k \ (1 + \cos(n x - \delta))$$

The potential type keyword is `torsion`, and the parameters `k` and `delta`
($\delta$) should be provided as strings, and `n` should be provided as an
integer.

```toml
[[dihedrals]]
atoms = ["C", "C", "C", "C"]
torsion = {k = "40 kJ/mol", delta = "120 deg", n: 4}
```

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

# Restrictions

# Potential computations

## Direct computation

## Cutoff computation

## Table interpolation

---

[Frenkel2002] Frenkel, D. & Smit, B. *Understanding molecular simulation.*
(Academic press, 2002).

[Wolf1999] Wolf, D., Keblinski, P., Phillpot, S. R. & Eggebrecht, J. *Exact
method for the simulation of Coulombic systems by spherically truncated,
pairwise 1/r summation.* The Journal of Chemical Physics **110**, 8254 (1999).
