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

Some force fields define additional restrictions concerning which particles
should interact together. For example, sometimes bonded particles should not
interact through electrostatic potential, or some interactions should only be
taken in account for particles not in the same molecule. The way to specify this
is to use restrictions. Restrictions can be used in two places: in the
`[[pairs]]` section, and in the `[coulomb]` section. In both cases, they are
specified with the `restriction` keyword, and one of the possible values.

```toml
[[pairs]]
atoms = ["O", "O"]
lj = {sigma = "3 A", epsilon = "123 kJ/mol"}
restriction = {scale14 = 0.5}

[coulomb]
ewald = {cutoff = "8 A", kmax = 6}
restriction = "intermolecular"
```

The possible values for `restriction` are:
- `"intramolecular"` or `"intra-molecular"` to restrict the potential to only
  particles in the same molecule;
- `"intermolecular"` or `"inter-molecular"` to restrict the potential to only
  particles NOT in the same molecule;
- `"exclude12"` to exclude particles directly bonded together;
- `"exclude13"` to exclude particles directly bonded together or forming an
  angle;
- `"exclude14"` to exclude particles directly bonded together; forming an angle
  or a dihedral angle;
- `{scale14 = <scaling>}` to exclude particles directly bonded together or
   forming an angle and to scale interaction between particles at 3 bonds of
   distance by the given `scaling` factor. The factor must be between 0 and 1.

# Potential computations

The same potential function (Lennard-Jones, Harmonic, *etc.*) can be computed
with different method: directly, by sifting it at the cutoff distance, using a
table interpolation, *etc.* This is the purpose of computation. The default way
is to use the mathematical function corresponding to a potential to compute it.
To use another computation, the `computation` keyword can be used in the
`[[pairs]]` section.

## Table interpolation

A way to compute a potential is to pre-compute it on a regularly spaced grid,
and then to interpolate values for points in the grid. In some cases, this can
be faster than recomputing the function every time.

This can be done with the `table` computation, which does a linear interpolation
in regularly spaced values in the `[0, max)` segment. You need to provide the
`max` value, and the number of points `n` for the interpolation:

```toml
[[pairs]]
atoms = ["O", "O"]
lj = {sigma = "3 A", epsilon = "123 kJ/mol"}
computation = {table = {max = "8 A", n = 5000}}
```

---

[Frenkel2002] Frenkel, D. & Smit, B. *Understanding molecular simulation.*
(Academic press, 2002).

[Wolf1999] Wolf, D., Keblinski, P., Phillpot, S. R. & Eggebrecht, J. *Exact
method for the simulation of Coulombic systems by spherically truncated,
pairwise 1/r summation.* The Journal of Chemical Physics **110**, 8254 (1999).
