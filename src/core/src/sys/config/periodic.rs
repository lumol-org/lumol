// Lumol, an extensible molecular simulation engine
// Copyright (C) Lumol's contributors — BSD license

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
    pub mass: f64,
    /// The covalent radius, in internal units
    pub covalent: f64,
    /// The Van der Waals radius, in internal units
    pub vdw: f64,
}

/// The full periodic table of elements is accessible here. This data comes from
/// the blue obelisk repository at http://svn.code.sf.net/p/bodr/code/trunk/bodr
static PERIODIC_TABLE: [ElementData; 119] = [
    ElementData{symbol: "Xx", name: "Dummy", mass: 0.0, covalent: 0.0, vdw: 0.0},
    ElementData{symbol: "H", name: "Hydrogen", mass: 1.008, covalent: 0.37, vdw: 1.2},
    ElementData{symbol: "He", name: "Helium", mass: 4.002602, covalent: 0.32, vdw: 1.4},
    ElementData{symbol: "Li", name: "Lithium", mass: 6.94, covalent: 1.34, vdw: 2.2},
    ElementData{symbol: "Be", name: "Beryllium", mass: 9.012182, covalent: 0.9, vdw: 1.9},
    ElementData{symbol: "B", name: "Boron", mass: 10.81, covalent: 0.82, vdw: 1.8},
    ElementData{symbol: "C", name: "Carbon", mass: 12.011, covalent: 0.77, vdw: 1.7},
    ElementData{symbol: "N", name: "Nitrogen", mass: 14.007, covalent: 0.75, vdw: 1.6},
    ElementData{symbol: "O", name: "Oxygen", mass: 15.999, covalent: 0.73, vdw: 1.55},
    ElementData{symbol: "F", name: "Fluorine", mass: 18.9984032, covalent: 0.71, vdw: 1.5},
    ElementData{symbol: "Ne", name: "Neon", mass: 20.1797, covalent: 0.69, vdw: 1.54},
    ElementData{symbol: "Na", name: "Sodium", mass: 22.98976928, covalent: 1.54, vdw: 2.4},
    ElementData{symbol: "Mg", name: "Magnesium", mass: 24.305, covalent: 1.3, vdw: 2.2},
    ElementData{symbol: "Al", name: "Aluminium", mass: 26.9815386, covalent: 1.18, vdw: 2.1},
    ElementData{symbol: "Si", name: "Silicon", mass: 28.085, covalent: 1.11, vdw: 2.1},
    ElementData{symbol: "P", name: "Phosphorus", mass: 30.973762, covalent: 1.06, vdw: 1.95},
    ElementData{symbol: "S", name: "Sulfur", mass: 32.06, covalent: 1.02, vdw: 1.8},
    ElementData{symbol: "Cl", name: "Chlorine", mass: 35.45, covalent: 0.99, vdw: 1.8},
    ElementData{symbol: "Ar", name: "Argon", mass: 39.948, covalent: 0.97, vdw: 1.88},
    ElementData{symbol: "K", name: "Potassium", mass: 39.0983, covalent: 1.96, vdw: 2.8},
    ElementData{symbol: "Ca", name: "Calcium", mass: 40.078, covalent: 1.74, vdw: 2.4},
    ElementData{symbol: "Sc", name: "Scandium", mass: 44.955912, covalent: 1.44, vdw: 2.3},
    ElementData{symbol: "Ti", name: "Titanium", mass: 47.867, covalent: 1.36, vdw: 2.15},
    ElementData{symbol: "V", name: "Vanadium", mass: 50.9415, covalent: 1.25, vdw: 2.05},
    ElementData{symbol: "Cr", name: "Chromium", mass: 51.9961, covalent: 1.27, vdw: 2.05},
    ElementData{symbol: "Mn", name: "Manganese", mass: 54.938045, covalent: 1.39, vdw: 2.05},
    ElementData{symbol: "Fe", name: "Iron", mass: 55.845, covalent: 1.25, vdw: 2.05},
    ElementData{symbol: "Co", name: "Cobalt", mass: 58.933195, covalent: 1.26, vdw: 2.0},
    ElementData{symbol: "Ni", name: "Nickel", mass: 58.6934, covalent: 1.21, vdw: 2.0},
    ElementData{symbol: "Cu", name: "Copper", mass: 63.546, covalent: 1.38, vdw: 2.0},
    ElementData{symbol: "Zn", name: "Zinc", mass: 65.38, covalent: 1.31, vdw: 2.1},
    ElementData{symbol: "Ga", name: "Gallium", mass: 69.723, covalent: 1.26, vdw: 2.1},
    ElementData{symbol: "Ge", name: "Germanium", mass: 72.63, covalent: 1.22, vdw: 2.1},
    ElementData{symbol: "As", name: "Arsenic", mass: 74.9216, covalent: 1.19, vdw: 2.05},
    ElementData{symbol: "Se", name: "Selenium", mass: 78.96, covalent: 1.16, vdw: 1.9},
    ElementData{symbol: "Br", name: "Bromine", mass: 79.904, covalent: 1.14, vdw: 1.9},
    ElementData{symbol: "Kr", name: "Krypton", mass: 83.798, covalent: 1.1, vdw: 2.02},
    ElementData{symbol: "Rb", name: "Rubidium", mass: 85.4678, covalent: 2.11, vdw: 2.9},
    ElementData{symbol: "Sr", name: "Strontium", mass: 87.62, covalent: 1.92, vdw: 2.55},
    ElementData{symbol: "Y", name: "Yttrium", mass: 88.90585, covalent: 1.62, vdw: 2.4},
    ElementData{symbol: "Zr", name: "Zirconium", mass: 91.224, covalent: 1.48, vdw: 2.3},
    ElementData{symbol: "Nb", name: "Niobium", mass: 92.90638, covalent: 1.37, vdw: 2.15},
    ElementData{symbol: "Mo", name: "Molybdenum", mass: 95.96, covalent: 1.45, vdw: 2.1},
    ElementData{symbol: "Tc", name: "Technetium", mass: 97.0, covalent: 1.56, vdw: 2.05},
    ElementData{symbol: "Ru", name: "Ruthenium", mass: 101.07, covalent: 1.26, vdw: 2.05},
    ElementData{symbol: "Rh", name: "Rhodium", mass: 102.9055, covalent: 1.35, vdw: 2.0},
    ElementData{symbol: "Pd", name: "Palladium", mass: 106.42, covalent: 1.31, vdw: 2.05},
    ElementData{symbol: "Ag", name: "Silver", mass: 107.8682, covalent: 1.53, vdw: 2.1},
    ElementData{symbol: "Cd", name: "Cadmium", mass: 112.411, covalent: 1.48, vdw: 2.2},
    ElementData{symbol: "In", name: "Indium", mass: 114.818, covalent: 1.44, vdw: 2.2},
    ElementData{symbol: "Sn", name: "Tin", mass: 118.71, covalent: 1.41, vdw: 2.25},
    ElementData{symbol: "Sb", name: "Antimony", mass: 121.76, covalent: 1.38, vdw: 2.2},
    ElementData{symbol: "Te", name: "Tellurium", mass: 127.6, covalent: 1.35, vdw: 2.1},
    ElementData{symbol: "I", name: "Iodine", mass: 126.90447, covalent: 1.33, vdw: 2.1},
    ElementData{symbol: "Xe", name: "Xenon", mass: 131.293, covalent: 1.3, vdw: 2.16},
    ElementData{symbol: "Cs", name: "Caesium", mass: 132.9054519, covalent: 2.25, vdw: 3.0},
    ElementData{symbol: "Ba", name: "Barium", mass: 137.327, covalent: 1.98, vdw: 2.7},
    ElementData{symbol: "La", name: "Lanthanum", mass: 138.90547, covalent: 1.69, vdw: 2.5},
    ElementData{symbol: "Ce", name: "Cerium", mass: 140.116, covalent: 1.69, vdw: 2.48},
    ElementData{symbol: "Pr", name: "Praseodymium", mass: 140.90765, covalent: 1.69, vdw: 2.47},
    ElementData{symbol: "Nd", name: "Neodymium", mass: 144.242, covalent: 1.69, vdw: 2.45},
    ElementData{symbol: "Pm", name: "Promethium", mass: 145.0, covalent: 1.69, vdw: 2.43},
    ElementData{symbol: "Sm", name: "Samarium", mass: 150.36, covalent: 1.69, vdw: 2.42},
    ElementData{symbol: "Eu", name: "Europium", mass: 151.964, covalent: 1.69, vdw: 2.4},
    ElementData{symbol: "Gd", name: "Gadolinium", mass: 157.25, covalent: 1.69, vdw: 2.38},
    ElementData{symbol: "Tb", name: "Terbium", mass: 158.92535, covalent: 1.69, vdw: 2.37},
    ElementData{symbol: "Dy", name: "Dysprosium", mass: 162.5, covalent: 1.69, vdw: 2.35},
    ElementData{symbol: "Ho", name: "Holmium", mass: 164.93032, covalent: 1.69, vdw: 2.33},
    ElementData{symbol: "Er", name: "Erbium", mass: 167.259, covalent: 1.69, vdw: 2.32},
    ElementData{symbol: "Tm", name: "Thulium", mass: 168.93421, covalent: 1.69, vdw: 2.3},
    ElementData{symbol: "Yb", name: "Ytterbium", mass: 173.054, covalent: 1.69, vdw: 2.28},
    ElementData{symbol: "Lu", name: "Lutetium", mass: 174.9668, covalent: 1.6, vdw: 2.27},
    ElementData{symbol: "Hf", name: "Hafnium", mass: 178.49, covalent: 1.5, vdw: 2.25},
    ElementData{symbol: "Ta", name: "Tantalum", mass: 180.94788, covalent: 1.38, vdw: 2.2},
    ElementData{symbol: "W", name: "Tungsten", mass: 183.84, covalent: 1.46, vdw: 2.1},
    ElementData{symbol: "Re", name: "Rhenium", mass: 186.207, covalent: 1.59, vdw: 2.05},
    ElementData{symbol: "Os", name: "Osmium", mass: 190.23, covalent: 1.28, vdw: 2.0},
    ElementData{symbol: "Ir", name: "Iridium", mass: 192.217, covalent: 1.37, vdw: 2.0},
    ElementData{symbol: "Pt", name: "Platinum", mass: 195.084, covalent: 1.28, vdw: 2.05},
    ElementData{symbol: "Au", name: "Gold", mass: 196.966569, covalent: 1.44, vdw: 2.1},
    ElementData{symbol: "Hg", name: "Mercury", mass: 200.592, covalent: 1.49, vdw: 2.05},
    ElementData{symbol: "Tl", name: "Thallium", mass: 204.38, covalent: 1.48, vdw: 2.2},
    ElementData{symbol: "Pb", name: "Lead", mass: 207.2, covalent: 1.47, vdw: 2.3},
    ElementData{symbol: "Bi", name: "Bismuth", mass: 208.9804, covalent: 1.46, vdw: 2.3},
    ElementData{symbol: "Po", name: "Polonium", mass: 209.0, covalent: 1.46, vdw: 2.0},
    ElementData{symbol: "At", name: "Astatine", mass: 210.0, covalent: 1.46, vdw: 2.0},
    ElementData{symbol: "Rn", name: "Radon", mass: 222.0, covalent: 1.45, vdw: 2.0},
    ElementData{symbol: "Fr", name: "Francium", mass: 223.0, covalent: 1.45, vdw: 2.0},
    ElementData{symbol: "Ra", name: "Radium", mass: 226.0, covalent: 1.45, vdw: 2.0},
    ElementData{symbol: "Ac", name: "Actinium", mass: 227.0, covalent: 1.45, vdw: 2.0},
    ElementData{symbol: "Th", name: "Thorium", mass: 232.03806, covalent: 1.45, vdw: 2.4},
    ElementData{symbol: "Pa", name: "Protactinium", mass: 231.03588, covalent: 1.45, vdw: 2.0},
    ElementData{symbol: "U", name: "Uranium", mass: 238.02891, covalent: 1.45, vdw: 2.3},
    ElementData{symbol: "Np", name: "Neptunium", mass: 237.0, covalent: 1.45, vdw: 2.0},
    ElementData{symbol: "Pu", name: "Plutonium", mass: 244.0, covalent: 1.45, vdw: 2.0},
    ElementData{symbol: "Am", name: "Americium", mass: 243.0, covalent: 1.45, vdw: 2.0},
    ElementData{symbol: "Cm", name: "Curium", mass: 247.0, covalent: 1.45, vdw: 2.0},
    ElementData{symbol: "Bk", name: "Berkelium", mass: 247.0, covalent: 1.45, vdw: 2.0},
    ElementData{symbol: "Cf", name: "Californium", mass: 251.0, covalent: 1.45, vdw: 2.0},
    ElementData{symbol: "Es", name: "Einsteinium", mass: 252.0, covalent: 1.45, vdw: 2.0},
    ElementData{symbol: "Fm", name: "Fermium", mass: 257.0, covalent: 1.45, vdw: 2.0},
    ElementData{symbol: "Md", name: "Mendelevium", mass: 258.0, covalent: 1.45, vdw: 2.0},
    ElementData{symbol: "No", name: "Nobelium", mass: 259.0, covalent: 1.45, vdw: 2.0},
    ElementData{symbol: "Lr", name: "Lawrencium", mass: 262.0, covalent: 1.45, vdw: 2.0},
    ElementData{symbol: "Rf", name: "Rutherfordium", mass: 267.0, covalent: 1.45, vdw: 2.0},
    ElementData{symbol: "Db", name: "Dubnium", mass: 270.0, covalent: 1.45, vdw: 2.0},
    ElementData{symbol: "Sg", name: "Seaborgium", mass: 271.0, covalent: 1.45, vdw: 2.0},
    ElementData{symbol: "Bh", name: "Bohrium", mass: 270.0, covalent: 1.45, vdw: 2.0},
    ElementData{symbol: "Hs", name: "Hassium", mass: 277.0, covalent: 1.45, vdw: 2.0},
    ElementData{symbol: "Mt", name: "Meitnerium", mass: 276.0, covalent: 1.45, vdw: 2.0},
    ElementData{symbol: "Ds", name: "Darmstadtium", mass: 281.0, covalent: 1.45, vdw: 2.0},
    ElementData{symbol: "Rg", name: "Roentgenium", mass: 282.0, covalent: 1.45, vdw: 2.0},
    ElementData{symbol: "Cn", name: "Copernicium", mass: 285.0, covalent: 1.45, vdw: 2.0},
    ElementData{symbol: "Uut", name: "Ununtrium", mass: 285.0, covalent: 1.45, vdw: 2.0},
    ElementData{symbol: "Fl", name: "Flerovium", mass: 289.0, covalent: 1.45, vdw: 2.0},
    ElementData{symbol: "Uup", name: "Ununpentium", mass: 289.0, covalent: 1.45, vdw: 2.0},
    ElementData{symbol: "Lv", name: "Livermorium", mass: 293.0, covalent: 1.45, vdw: 2.0},
    ElementData{symbol: "Uus", name: "Ununseptium", mass: 294.0, covalent: 1.45, vdw: 2.0},
    ElementData{symbol: "Uuo", name: "Ununoctium", mass: 294.0, covalent: 1.45, vdw: 2.0},
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
    pub fn mass<S: AsRef<str>>(name: S) -> Option<f64> {
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
    pub fn covalent<S: AsRef<str>>(name: S) -> Option<f64> {
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
    pub fn vdw<S: AsRef<str>>(name: S) -> Option<f64> {
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
        assert_eq!(PeriodicTable::mass("O"), Some(15.999));
        assert_eq!(PeriodicTable::name("O"), Some("Oxygen"));
        assert_eq!(PeriodicTable::covalent("O"), Some(0.73));
        assert_eq!(PeriodicTable::vdw("O"), Some(1.55));

        // Errors on non-existing values
        assert_eq!(PeriodicTable::mass("HOH"), None);
        assert_eq!(PeriodicTable::name("HOH"), None);
        assert_eq!(PeriodicTable::covalent("HOH"), None);
        assert_eq!(PeriodicTable::vdw("HOH"), None);
    }

    #[test]
    fn add_elements() {
        // Add a new element
        let element = ElementData{symbol: "Ooo", name: "Ooo", mass: 0.0, covalent: 0.0, vdw: 0.0};
        PeriodicTable::add_element(element);

        assert_eq!(PeriodicTable::mass("Ooo"), Some(0.0));
        assert_eq!(PeriodicTable::name("Ooo"), Some("Ooo"));
        assert_eq!(PeriodicTable::covalent("Ooo"), Some(0.0));
        assert_eq!(PeriodicTable::vdw("Ooo"), Some(0.0));

        // Overwrite existing element
        let element = ElementData{symbol: "H", name: "New H", mass: 0.0, covalent: 0.0, vdw: 0.0};
        PeriodicTable::add_element(element);

        assert_eq!(PeriodicTable::mass("H"), Some(0.0));
        assert_eq!(PeriodicTable::name("H"), Some("New H"));
        assert_eq!(PeriodicTable::covalent("H"), Some(0.0));
        assert_eq!(PeriodicTable::vdw("H"), Some(0.0));
    }
}
