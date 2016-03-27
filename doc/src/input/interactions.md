# Interactions

The interactions in a system can be specified using YAML syntax. An
interaction is the association of atomic types, a potential function and
optionally restrictions and mode of computation.

The `pairs`, `bonds`, `angles` and `dihedrals` sections are arrays of
interactions, listing all possible atomic type combinations and the
associated potentials. The `coulomb` section contains information about the
treatment of long-range electrostatic interactions and the associated
charges.

An example of an input file for the f-SPC model of water is given bellow:

```yaml
# Non-bonded pairs
pairs:
    - atoms: [O, O]
      type: LennardJones
      sigma: 3.16 A
      epsilon: 0.155 kcal/mol
    - atoms: [O, H]
      type: Null
    - atoms: [H, H]
      type: harmonic
      k: 79.8 kcal/mol/A^2
      x0: 1.633 A
      restriction:
          type: IntraMolecular

# Molecular interactions
bonds:
    - atoms: [O, H]
      type: Harmonic
      k: 1054.2 kcal/mol/A^2
      x0: 1.0 A
angles:
    - atoms: [H, O, H]
      type: Harmonic
      k: 75.9 kcal/mol/rad
      x0: 109.5 deg

# Coulombic interactions
coulomb:
    type: Ewald
    cutoff: 10 A
    kmax: 7
    charges:
        O: -0.82
        H: 0.41
```

## Pairs and molecular interactions

The `pairs`, `bonds`, `angles` and `dihedrals` sections are arrays, in which
every entry must contain the `atoms` and `type` keys. The `atoms` key specifies
the atom types to which the interaction should be applied and the `type` key
the [potential](input/potentials.html#Available%20potentials) to use. The
`atoms` key should contain two atoms for `pairs` and `bonds`, three atoms for
`angles` and four atoms for `dihedrals`.

For example, to use a harmonic bond potential for all `C-H` bonds:

```yaml
bonds:
    - atoms: [C, H]
      type: Harmonic
      # Theses keys are specific to the Harmonic potential
      x0: 3.405 A
      k: 2385 kcal/mol/A^2
```

It is also possible to specify an additional
[restrictions](input/potentials.html#Restrictions) using the `restriction` key;
or a [computation method](input/potentials.html#Potential%20computations) using
the `computation` key.

```yaml
pairs:
    - atoms: [Ar, Ar]
      type: LennardJones
      sigma: 3.405 A
      epsilon: 0.2385 kcal/mol
      computation:
          type: table
          numpoints: 2000
          max: 20.0 A
      restriction:
         type: IntraMolecular
```

## Coulombic interactions

The method for treatment of coulombic interactions is specified in the `coulomb`
section, and should contain the `type` and the `charges` keys. The `type` key
indicates the method to use for computing the electrostatic energy, and the
`charges` gives the partial charges for all the atom types. The available
methods are described in the [coulombic
potential](input/potentials.html#Electrostatic%20interactions) section. An
additional [restrictions](input/potentials.html#Restrictions) can be used with
the `restriction` key.

```yaml
coulomb:
    type: Ewald
    charges:
        Na: 1
        Cl: -1
    # These values are specific to Ewald method
    cutoff: 10 A
    kmax: 7
```
