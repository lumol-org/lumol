# System

A `System` contains all the data about the physical system we are simulating.
It contains four types of data:

- A list of **Particles**, which are physical objects, with a position, a
  velocity, a mass and a name;
- A list of **Molecules** containing information about how the particles are
  bounded together;
- An **UnitCell**, *i.e.* the bounding box of the simulation.
- **Interactions**, sometimes called a force-field.

![System components](static/img/system.svg#center)

## Unit cells

Lumol knows about three types of unit cells:

* Infinite cells do not have any boundaries, and can be used to simulate
  systems in vacuum;
* Orthorhombic cells have up to three independent lengths, and all the angles
  of the cell are set to 90°;
* Triclinic cells have 6 independent parameters: 3 lengths and 3 angles.

Orthorhombic and Triclinic cells are used in combination with [periodic boundary
conditions](https://en.wikipedia.org/wiki/Periodic_boundary_conditions) to
simulate infinite systems.

## Interactions

Interactions associate a potential and some particles kind. Lumol provides
code for various potentials types:

* Non-bonding pair potentials;
* Bonds potentials in molecules;
* Angles potentials in molecules;
* Dihedral angles potentials in molecules;
* Long-ranges coulombic potentials (Ewald and Wolf methods);
* Arbitrary external potential applying on the whole system at once.
