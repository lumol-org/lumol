Units
=====

The unit system used internally by Lumol is the following:

- Angstrom (``A``) for distances;
- Femtosecond (``fs``) for time;
- Unified atomic mass unit (``u`` or ``Da``) for mass;
- Kelvin (``K``) for temperature;
- Number of particles for quantity of matter;
- Radian (``rad``) for angles;

Any other internal unit is derived from this set:

- The internal unit of energy is ``u A^2 fs^-2``;
- The internal unit of force is ``u A fs^-2``;
- The internal unit of pressure is ``u A^-1 fs^-2``;
- *etc.*

For convenience Lumol provides convertion facilities for any value in these
internal unit to any others units. The following table lists available units
that can be converted:

+------------+-----------------------------+
| Quantity   | Accepted units              |
+============+=============================+
| Distance   | A, Å, nm, pm, fm, m, bohr   |
+------------+-----------------------------+
| Time       | fs, ps, ns                  |
+------------+-----------------------------+
| Mass       | u, Da, kDa, g, kg           |
+------------+-----------------------------+
| Matter     | mol                         |
+------------+-----------------------------+
| Angle      | rad, deg                    |
+------------+-----------------------------+
| Energy     | J, kJ, kcal, eV, H, Ry      |
+------------+-----------------------------+
| Force      | N                           |
+------------+-----------------------------+
| Pressure   | Pa, kPa, MPa, bar, atm      |
+------------+-----------------------------+

In the input files, the units are specified as strings, and must be
spelled exactly as in the above table. They can be combined with other
units using ``*`` for multiplication, ``/`` for division, and ``^`` for
exponentiation. Parentheses can be used to group sub-units together. Some
valid unit strings are ``kcal/mol``, ``(J / mol) * A^-2``, and
``m*fs^-1``.
