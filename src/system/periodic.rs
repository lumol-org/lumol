// Cymbalum, an extensible molecular simulation engine
// Copyright (C) 2015-2016 G. Fraux — BSD license

//! The periodic module give access to information about "extended" elements.
//! These elements can be real element like Hydrogen or Carbon, or composed
//! elements like methane or water.
//!
//! There are two way to access this information. The `PERIODIC_TABLE` array
//! contains all the data about real elements. The `PeriodicTable` struct
//! provides helper functions, and a way to add new elements in the list.
use std::sync::RwLock;

/// Data about one "extended" chemical element.
pub struct ElementData {
    /// The symbol is used to identify the element
    pub symbol: &'static str,
    /// The element full name
    pub name: &'static str,
    /// The element mass, in internal units
    pub mass: f32,
    /// The covalent radius, in internal units
    pub covalent: f32,
    /// The Van der Waals radius, in internal units
    pub vdw: f32,
}

/// The full periodic table of elements is accessible here. This data comes from
/// the blue obelisk repository at http://svn.code.sf.net/p/bodr/code/trunk/bodr
static PERIODIC_TABLE: [ElementData; 119] = [
    ElementData{symbol: "Xx", name: "Dummy", mass: 0.0f32, covalent: 0.0f32, vdw: 0.0f32},
    ElementData{symbol: "H", name: "Hydrogen", mass: 1.008f32, covalent: 0.37f32, vdw: 1.2f32},
    ElementData{symbol: "He", name: "Helium", mass: 4.002602f32, covalent: 0.32f32, vdw: 1.4f32},
    ElementData{symbol: "Li", name: "Lithium", mass: 6.94f32, covalent: 1.34f32, vdw: 2.2f32},
    ElementData{symbol: "Be", name: "Beryllium", mass: 9.012182f32, covalent: 0.9f32, vdw: 1.9f32},
    ElementData{symbol: "B", name: "Boron", mass: 10.81f32, covalent: 0.82f32, vdw: 1.8f32},
    ElementData{symbol: "C", name: "Carbon", mass: 12.011f32, covalent: 0.77f32, vdw: 1.7f32},
    ElementData{symbol: "N", name: "Nitrogen", mass: 14.007f32, covalent: 0.75f32, vdw: 1.6f32},
    ElementData{symbol: "O", name: "Oxygen", mass: 15.999f32, covalent: 0.73f32, vdw: 1.55f32},
    ElementData{symbol: "F", name: "Fluorine", mass: 18.9984032f32, covalent: 0.71f32, vdw: 1.5f32},
    ElementData{symbol: "Ne", name: "Neon", mass: 20.1797f32, covalent: 0.69f32, vdw: 1.54f32},
    ElementData{symbol: "Na", name: "Sodium", mass: 22.98976928f32, covalent: 1.54f32, vdw: 2.4f32},
    ElementData{symbol: "Mg", name: "Magnesium", mass: 24.305f32, covalent: 1.3f32, vdw: 2.2f32},
    ElementData{symbol: "Al", name: "Aluminium", mass: 26.9815386f32, covalent: 1.18f32, vdw: 2.1f32},
    ElementData{symbol: "Si", name: "Silicon", mass: 28.085f32, covalent: 1.11f32, vdw: 2.1f32},
    ElementData{symbol: "P", name: "Phosphorus", mass: 30.973762f32, covalent: 1.06f32, vdw: 1.95f32},
    ElementData{symbol: "S", name: "Sulfur", mass: 32.06f32, covalent: 1.02f32, vdw: 1.8f32},
    ElementData{symbol: "Cl", name: "Chlorine", mass: 35.45f32, covalent: 0.99f32, vdw: 1.8f32},
    ElementData{symbol: "Ar", name: "Argon", mass: 39.948f32, covalent: 0.97f32, vdw: 1.88f32},
    ElementData{symbol: "K", name: "Potassium", mass: 39.0983f32, covalent: 1.96f32, vdw: 2.8f32},
    ElementData{symbol: "Ca", name: "Calcium", mass: 40.078f32, covalent: 1.74f32, vdw: 2.4f32},
    ElementData{symbol: "Sc", name: "Scandium", mass: 44.955912f32, covalent: 1.44f32, vdw: 2.3f32},
    ElementData{symbol: "Ti", name: "Titanium", mass: 47.867f32, covalent: 1.36f32, vdw: 2.15f32},
    ElementData{symbol: "V", name: "Vanadium", mass: 50.9415f32, covalent: 1.25f32, vdw: 2.05f32},
    ElementData{symbol: "Cr", name: "Chromium", mass: 51.9961f32, covalent: 1.27f32, vdw: 2.05f32},
    ElementData{symbol: "Mn", name: "Manganese", mass: 54.938045f32, covalent: 1.39f32, vdw: 2.05f32},
    ElementData{symbol: "Fe", name: "Iron", mass: 55.845f32, covalent: 1.25f32, vdw: 2.05f32},
    ElementData{symbol: "Co", name: "Cobalt", mass: 58.933195f32, covalent: 1.26f32, vdw: 2.0f32},
    ElementData{symbol: "Ni", name: "Nickel", mass: 58.6934f32, covalent: 1.21f32, vdw: 2.0f32},
    ElementData{symbol: "Cu", name: "Copper", mass: 63.546f32, covalent: 1.38f32, vdw: 2.0f32},
    ElementData{symbol: "Zn", name: "Zinc", mass: 65.38f32, covalent: 1.31f32, vdw: 2.1f32},
    ElementData{symbol: "Ga", name: "Gallium", mass: 69.723f32, covalent: 1.26f32, vdw: 2.1f32},
    ElementData{symbol: "Ge", name: "Germanium", mass: 72.63f32, covalent: 1.22f32, vdw: 2.1f32},
    ElementData{symbol: "As", name: "Arsenic", mass: 74.9216f32, covalent: 1.19f32, vdw: 2.05f32},
    ElementData{symbol: "Se", name: "Selenium", mass: 78.96f32, covalent: 1.16f32, vdw: 1.9f32},
    ElementData{symbol: "Br", name: "Bromine", mass: 79.904f32, covalent: 1.14f32, vdw: 1.9f32},
    ElementData{symbol: "Kr", name: "Krypton", mass: 83.798f32, covalent: 1.1f32, vdw: 2.02f32},
    ElementData{symbol: "Rb", name: "Rubidium", mass: 85.4678f32, covalent: 2.11f32, vdw: 2.9f32},
    ElementData{symbol: "Sr", name: "Strontium", mass: 87.62f32, covalent: 1.92f32, vdw: 2.55f32},
    ElementData{symbol: "Y", name: "Yttrium", mass: 88.90585f32, covalent: 1.62f32, vdw: 2.4f32},
    ElementData{symbol: "Zr", name: "Zirconium", mass: 91.224f32, covalent: 1.48f32, vdw: 2.3f32},
    ElementData{symbol: "Nb", name: "Niobium", mass: 92.90638f32, covalent: 1.37f32, vdw: 2.15f32},
    ElementData{symbol: "Mo", name: "Molybdenum", mass: 95.96f32, covalent: 1.45f32, vdw: 2.1f32},
    ElementData{symbol: "Tc", name: "Technetium", mass: 97.0f32, covalent: 1.56f32, vdw: 2.05f32},
    ElementData{symbol: "Ru", name: "Ruthenium", mass: 101.07f32, covalent: 1.26f32, vdw: 2.05f32},
    ElementData{symbol: "Rh", name: "Rhodium", mass: 102.9055f32, covalent: 1.35f32, vdw: 2.0f32},
    ElementData{symbol: "Pd", name: "Palladium", mass: 106.42f32, covalent: 1.31f32, vdw: 2.05f32},
    ElementData{symbol: "Ag", name: "Silver", mass: 107.8682f32, covalent: 1.53f32, vdw: 2.1f32},
    ElementData{symbol: "Cd", name: "Cadmium", mass: 112.411f32, covalent: 1.48f32, vdw: 2.2f32},
    ElementData{symbol: "In", name: "Indium", mass: 114.818f32, covalent: 1.44f32, vdw: 2.2f32},
    ElementData{symbol: "Sn", name: "Tin", mass: 118.71f32, covalent: 1.41f32, vdw: 2.25f32},
    ElementData{symbol: "Sb", name: "Antimony", mass: 121.76f32, covalent: 1.38f32, vdw: 2.2f32},
    ElementData{symbol: "Te", name: "Tellurium", mass: 127.6f32, covalent: 1.35f32, vdw: 2.1f32},
    ElementData{symbol: "I", name: "Iodine", mass: 126.90447f32, covalent: 1.33f32, vdw: 2.1f32},
    ElementData{symbol: "Xe", name: "Xenon", mass: 131.293f32, covalent: 1.3f32, vdw: 2.16f32},
    ElementData{symbol: "Cs", name: "Caesium", mass: 132.9054519f32, covalent: 2.25f32, vdw: 3.0f32},
    ElementData{symbol: "Ba", name: "Barium", mass: 137.327f32, covalent: 1.98f32, vdw: 2.7f32},
    ElementData{symbol: "La", name: "Lanthanum", mass: 138.90547f32, covalent: 1.69f32, vdw: 2.5f32},
    ElementData{symbol: "Ce", name: "Cerium", mass: 140.116f32, covalent: 1.69f32, vdw: 2.48f32},
    ElementData{symbol: "Pr", name: "Praseodymium", mass: 140.90765f32, covalent: 1.69f32, vdw: 2.47f32},
    ElementData{symbol: "Nd", name: "Neodymium", mass: 144.242f32, covalent: 1.69f32, vdw: 2.45f32},
    ElementData{symbol: "Pm", name: "Promethium", mass: 145.0f32, covalent: 1.69f32, vdw: 2.43f32},
    ElementData{symbol: "Sm", name: "Samarium", mass: 150.36f32, covalent: 1.69f32, vdw: 2.42f32},
    ElementData{symbol: "Eu", name: "Europium", mass: 151.964f32, covalent: 1.69f32, vdw: 2.4f32},
    ElementData{symbol: "Gd", name: "Gadolinium", mass: 157.25f32, covalent: 1.69f32, vdw: 2.38f32},
    ElementData{symbol: "Tb", name: "Terbium", mass: 158.92535f32, covalent: 1.69f32, vdw: 2.37f32},
    ElementData{symbol: "Dy", name: "Dysprosium", mass: 162.5f32, covalent: 1.69f32, vdw: 2.35f32},
    ElementData{symbol: "Ho", name: "Holmium", mass: 164.93032f32, covalent: 1.69f32, vdw: 2.33f32},
    ElementData{symbol: "Er", name: "Erbium", mass: 167.259f32, covalent: 1.69f32, vdw: 2.32f32},
    ElementData{symbol: "Tm", name: "Thulium", mass: 168.93421f32, covalent: 1.69f32, vdw: 2.3f32},
    ElementData{symbol: "Yb", name: "Ytterbium", mass: 173.054f32, covalent: 1.69f32, vdw: 2.28f32},
    ElementData{symbol: "Lu", name: "Lutetium", mass: 174.9668f32, covalent: 1.6f32, vdw: 2.27f32},
    ElementData{symbol: "Hf", name: "Hafnium", mass: 178.49f32, covalent: 1.5f32, vdw: 2.25f32},
    ElementData{symbol: "Ta", name: "Tantalum", mass: 180.94788f32, covalent: 1.38f32, vdw: 2.2f32},
    ElementData{symbol: "W", name: "Tungsten", mass: 183.84f32, covalent: 1.46f32, vdw: 2.1f32},
    ElementData{symbol: "Re", name: "Rhenium", mass: 186.207f32, covalent: 1.59f32, vdw: 2.05f32},
    ElementData{symbol: "Os", name: "Osmium", mass: 190.23f32, covalent: 1.28f32, vdw: 2.0f32},
    ElementData{symbol: "Ir", name: "Iridium", mass: 192.217f32, covalent: 1.37f32, vdw: 2.0f32},
    ElementData{symbol: "Pt", name: "Platinum", mass: 195.084f32, covalent: 1.28f32, vdw: 2.05f32},
    ElementData{symbol: "Au", name: "Gold", mass: 196.966569f32, covalent: 1.44f32, vdw: 2.1f32},
    ElementData{symbol: "Hg", name: "Mercury", mass: 200.592f32, covalent: 1.49f32, vdw: 2.05f32},
    ElementData{symbol: "Tl", name: "Thallium", mass: 204.38f32, covalent: 1.48f32, vdw: 2.2f32},
    ElementData{symbol: "Pb", name: "Lead", mass: 207.2f32, covalent: 1.47f32, vdw: 2.3f32},
    ElementData{symbol: "Bi", name: "Bismuth", mass: 208.9804f32, covalent: 1.46f32, vdw: 2.3f32},
    ElementData{symbol: "Po", name: "Polonium", mass: 209.0f32, covalent: 1.46f32, vdw: 2.0f32},
    ElementData{symbol: "At", name: "Astatine", mass: 210.0f32, covalent: 1.46f32, vdw: 2.0f32},
    ElementData{symbol: "Rn", name: "Radon", mass: 222.0f32, covalent: 1.45f32, vdw: 2.0f32},
    ElementData{symbol: "Fr", name: "Francium", mass: 223.0f32, covalent: 1.45f32, vdw: 2.0f32},
    ElementData{symbol: "Ra", name: "Radium", mass: 226.0f32, covalent: 1.45f32, vdw: 2.0f32},
    ElementData{symbol: "Ac", name: "Actinium", mass: 227.0f32, covalent: 1.45f32, vdw: 2.0f32},
    ElementData{symbol: "Th", name: "Thorium", mass: 232.03806f32, covalent: 1.45f32, vdw: 2.4f32},
    ElementData{symbol: "Pa", name: "Protactinium", mass: 231.03588f32, covalent: 1.45f32, vdw: 2.0f32},
    ElementData{symbol: "U", name: "Uranium", mass: 238.02891f32, covalent: 1.45f32, vdw: 2.3f32},
    ElementData{symbol: "Np", name: "Neptunium", mass: 237.0f32, covalent: 1.45f32, vdw: 2.0f32},
    ElementData{symbol: "Pu", name: "Plutonium", mass: 244.0f32, covalent: 1.45f32, vdw: 2.0f32},
    ElementData{symbol: "Am", name: "Americium", mass: 243.0f32, covalent: 1.45f32, vdw: 2.0f32},
    ElementData{symbol: "Cm", name: "Curium", mass: 247.0f32, covalent: 1.45f32, vdw: 2.0f32},
    ElementData{symbol: "Bk", name: "Berkelium", mass: 247.0f32, covalent: 1.45f32, vdw: 2.0f32},
    ElementData{symbol: "Cf", name: "Californium", mass: 251.0f32, covalent: 1.45f32, vdw: 2.0f32},
    ElementData{symbol: "Es", name: "Einsteinium", mass: 252.0f32, covalent: 1.45f32, vdw: 2.0f32},
    ElementData{symbol: "Fm", name: "Fermium", mass: 257.0f32, covalent: 1.45f32, vdw: 2.0f32},
    ElementData{symbol: "Md", name: "Mendelevium", mass: 258.0f32, covalent: 1.45f32, vdw: 2.0f32},
    ElementData{symbol: "No", name: "Nobelium", mass: 259.0f32, covalent: 1.45f32, vdw: 2.0f32},
    ElementData{symbol: "Lr", name: "Lawrencium", mass: 262.0f32, covalent: 1.45f32, vdw: 2.0f32},
    ElementData{symbol: "Rf", name: "Rutherfordium", mass: 267.0f32, covalent: 1.45f32, vdw: 2.0f32},
    ElementData{symbol: "Db", name: "Dubnium", mass: 270.0f32, covalent: 1.45f32, vdw: 2.0f32},
    ElementData{symbol: "Sg", name: "Seaborgium", mass: 271.0f32, covalent: 1.45f32, vdw: 2.0f32},
    ElementData{symbol: "Bh", name: "Bohrium", mass: 270.0f32, covalent: 1.45f32, vdw: 2.0f32},
    ElementData{symbol: "Hs", name: "Hassium", mass: 277.0f32, covalent: 1.45f32, vdw: 2.0f32},
    ElementData{symbol: "Mt", name: "Meitnerium", mass: 276.0f32, covalent: 1.45f32, vdw: 2.0f32},
    ElementData{symbol: "Ds", name: "Darmstadtium", mass: 281.0f32, covalent: 1.45f32, vdw: 2.0f32},
    ElementData{symbol: "Rg", name: "Roentgenium", mass: 282.0f32, covalent: 1.45f32, vdw: 2.0f32},
    ElementData{symbol: "Cn", name: "Copernicium", mass: 285.0f32, covalent: 1.45f32, vdw: 2.0f32},
    ElementData{symbol: "Uut", name: "Ununtrium", mass: 285.0f32, covalent: 1.45f32, vdw: 2.0f32},
    ElementData{symbol: "Fl", name: "Flerovium", mass: 289.0f32, covalent: 1.45f32, vdw: 2.0f32},
    ElementData{symbol: "Uup", name: "Ununpentium", mass: 289.0f32, covalent: 1.45f32, vdw: 2.0f32},
    ElementData{symbol: "Lv", name: "Livermorium", mass: 293.0f32, covalent: 1.45f32, vdw: 2.0f32},
    ElementData{symbol: "Uus", name: "Ununseptium", mass: 294.0f32, covalent: 1.45f32, vdw: 2.0f32},
    ElementData{symbol: "Uuo", name: "Ununoctium", mass: 294.0f32, covalent: 1.45f32, vdw: 2.0f32},
];

/// The `PeriodicTable` struct give access to elements information, and give a
/// way to register new elements in the list.
pub struct PeriodicTable;

lazy_static!(
    static ref OTHER_ELEMENTS: RwLock<Vec<ElementData>> = RwLock::new(Vec::new());
);

impl PeriodicTable {
    /// Add a new element in the elements list. The additional list is searched
    /// before the `PERIODIC_TABLE` array, so adding a new element here will
    /// make it overwrite the real one.
    ///
    /// This function use a RwLock to protect itself against race conditions,
    /// and will block. As such, it should only be called by one thread, and
    /// preferably before any parallel region.
    pub fn add_element(element: ElementData) {
        let mut others = OTHER_ELEMENTS.write().expect("RwLock is poisonned");
        others.push(element);
    }

    /// Get the name of an element from it's symbol.
    pub fn name<S: AsRef<str>>(name: S) -> Option<&'static str> {
        let symbol = name.as_ref();
        let others = OTHER_ELEMENTS.read().expect("RwLock is poisonned");
        for element in others.iter() {
            if element.symbol == symbol {
                return Some(element.name)
            }
        }
        for element in PERIODIC_TABLE.iter() {
            if element.symbol == symbol {
                return Some(element.name)
            }
        }
        return None
    }

    /// Get the mass of an element from it's symbol.
    pub fn mass<S: AsRef<str>>(name: S) -> Option<f32> {
        let symbol = name.as_ref();
        let others = OTHER_ELEMENTS.read().expect("RwLock is poisonned");
        for element in others.iter() {
            if element.symbol == symbol {
                return Some(element.mass)
            }
        }
        for element in PERIODIC_TABLE.iter() {
            if element.symbol == symbol {
                return Some(element.mass)
            }
        }
        return None
    }

    /// Get the covalent radius of an element from it's symbol.
    pub fn covalent<S: AsRef<str>>(name: S) -> Option<f32> {
        let symbol = name.as_ref();
        let others = OTHER_ELEMENTS.read().expect("RwLock is poisonned");
        for element in others.iter() {
            if element.symbol == symbol {
                return Some(element.covalent)
            }
        }
        for element in PERIODIC_TABLE.iter() {
            if element.symbol == symbol {
                return Some(element.covalent)
            }
        }
        return None
    }

    /// Get the Van der Waals radius of an element from it's symbol.
    pub fn vdw<S: AsRef<str>>(name: S) -> Option<f32> {
        let symbol = name.as_ref();
        let others = OTHER_ELEMENTS.read().expect("RwLock is poisonned");
        for element in others.iter() {
            if element.symbol == symbol {
                return Some(element.vdw)
            }
        }
        for element in PERIODIC_TABLE.iter() {
            if element.symbol == symbol {
                return Some(element.vdw)
            }
        }
        return None
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn values() {
        // Get existing values
        assert_eq!(PeriodicTable::mass("O"), Some(15.999f32));
        assert_eq!(PeriodicTable::name("O"), Some("Oxygen"));
        assert_eq!(PeriodicTable::covalent("O"), Some(0.73f32));
        assert_eq!(PeriodicTable::vdw("O"), Some(1.55f32));

        // Errors on non-existing values
        assert_eq!(PeriodicTable::mass("HOH"), None);
        assert_eq!(PeriodicTable::name("HOH"), None);
        assert_eq!(PeriodicTable::covalent("HOH"), None);
        assert_eq!(PeriodicTable::vdw("HOH"), None);
    }

    #[test]
    fn add_elements() {
        // Add a new element
        let element = ElementData{symbol: "Ooo", name: "Ooo", mass: 0.0f32, covalent: 0.0f32, vdw: 0.0f32};
        PeriodicTable::add_element(element);

        assert_eq!(PeriodicTable::mass("Ooo"), Some(0.0f32));
        assert_eq!(PeriodicTable::name("Ooo"), Some("Ooo"));
        assert_eq!(PeriodicTable::covalent("Ooo"), Some(0.0f32));
        assert_eq!(PeriodicTable::vdw("Ooo"), Some(0.0f32));

        // Overwrite existing element
        let element = ElementData{symbol: "H", name: "New H", mass: 0.0f32, covalent: 0.0f32, vdw: 0.0f32};
        PeriodicTable::add_element(element);

        assert_eq!(PeriodicTable::mass("H"), Some(0.0f32));
        assert_eq!(PeriodicTable::name("H"), Some("New H"));
        assert_eq!(PeriodicTable::covalent("H"), Some(0.0f32));
        assert_eq!(PeriodicTable::vdw("H"), Some(0.0f32));
    }
}
