# Available potentials

This section is a list of all the available potentials in Cymbalum, with the
associated `type` keyword and parameters.

The same potential can be used for either pairs (at distance $r$); or for
angles (at angle $\phi$). In all the formula, the $x$ parameter
represents either a distance or an angle.

## Null potential

This potential is always 0, for all values of $x$. This potential should be
used to remove interactions between atoms in a pair/bond/angle/dihedral which is
present in the system but should not be interacting together.

The potential type keyword is `Null`. Please be aware that this keyword can not
be written as `null`, because that a special value in YAML.

## Lennard-Jones potential

Lennard-Jones potential is a classical potential for pair interactions expressed
as: $$ V(x) = 4 \epsilon \left[\left(\frac{\sigma}{x}\right)^{12} -
\left(\frac{\sigma}{x}\right)^6\right].$$

The potential type keyword is `LennardJones`, and the `sigma` ($\sigma$) and
`epsilon` ($\epsilon$) parameters should be provided as keys for this potential.

```yaml
pairs:
    - atoms: [O, O]
      type: LennardJones
      sigma: 3.16 A
      epsilon: 0.155 kcal/mol
```

## Harmonic potential

The Harmonic potential is usually used for molecular interactions: bonds, angles
or dihedrals. It is expressed as: $$ V(x) = \frac 12 k \ (x - x_0)^2$$

The potential type keyword is `Harmonic`, and the `k` and `x0` parameters should
be provided as keys.

```yaml
bonds:
    - atoms: [O, H]
      type: Harmonic
      k: 1054.2 kcal/mol/A^2
      x0: 1.0 A
angles:
    - atoms: [H, O, H]
      type: Harmonic
      k: 75.9 kcal/mol/rad^2
      x0: 109.5 deg
```

## Cosine-Harmonic potential

This potential is usually used for angles and dihedrals angles interactions,
because it presents a $2\pi$ periodicity. It is expressed as: $$ V(x) = \frac 12
k \ (\cos x - \cos x_0)^2$$

The potential type keyword is `CosineHarmonic`, and the `k` and `x0` parameters
should be provided as keys.

```yaml
angles:
    - atoms: [H, C, H]
      type: CosineHarmonic
      k: 67 kJ/mol
      x0: 120 deg
```

## Torsion potential

This potential is usually used for dihedrals angles interactions. It is
expressed as: $$ V(x) = k \ (1 + \cos(n x - \delta))$$

The potential type keyword is `Torsion`, and the `k`, `n` and `delta` ($\delta$)
parameters should be provided as keys.

```yaml
dihedrals:
    - atoms: [C, C, C, C]
      type: Torsion
      k: 40 kJ/mol
      delta: 120 deg
      n: 4
```

# Electrostatic interactions

## Ewald method

## Wolf method

# Restrictions

# Potential computations

## Direct computation

## Cutoff computation

## Table interpolation
