# Available potentials

This section is a list of all the available potentials in Lumol, with the
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
