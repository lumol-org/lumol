# Units

Cymbalum features a set of internal units and offers facilities to convert
from and to this set of internal units.

The internal unit system is the following:
- Angstrom ($A$) for distances;
- Femto-second ($fs$) for time;
- Unified atomic mass unit ($u$ or $Da$) for mass;
- Kelvin ($K$) for temperature;
- Number of particles for quantity of matter;
- Radian ($rad$) for angles;

Any other internal unit is derived from this set:
- The internal unit of energy is $u A^2 fs^{-2}$;
- The internal unit of forces is $u A fs^{-2}$;
- The internal unit of pressure is $u A^{-1} fs^{-2}$;
- *etc.*

Cymbalum knows how to convert any value in these internal unit to others units.
The accepted units are:

| Quantity    | Accepted units                    |
|-------------|-----------------------------------|
| Distance    | A, &Aring;, nm, pm, fm, m, bohr   |
| Time        | fs, ps, ns                        |
| Mass        | u, Da, kDa, g, kg                 |
| Matter      | mol                               |
| Angle       | rad, deg                          |
| Energy      | J, kJ, kcal, eV, H, Ry            |
| Force       | N                                 |
| Pressure    | Pa, kPa, MPa, bar, atm            |

In the input files, the units are specified as strings, and must be spelled
exactly as in the above table. They can be combined with other units using
either a whitespace, `.` or `*` for multiplication, the `/` character for
division, and the `^` character for expontents.
For example, some valid strings are `kcal/mol`, `J * mol * A^-2`, `m.fs^{-1}`.
