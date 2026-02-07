//! Periodic table data and element lookup.

/// A chemical element from the periodic table.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Element {
    pub atomic_number: u8,
    pub symbol: &'static str,
    pub name: &'static str,
    pub atomic_weight: f64,
    pub valence: u8,
    pub max_bonds: u8,
}

/// Elements 1–54 (H through Xe).
static ELEMENTS: [Element; 54] = [
    Element { atomic_number: 1, symbol: "H", name: "Hydrogen", atomic_weight: 1.008, valence: 1, max_bonds: 1 },
    Element { atomic_number: 2, symbol: "He", name: "Helium", atomic_weight: 4.003, valence: 0, max_bonds: 0 },
    Element { atomic_number: 3, symbol: "Li", name: "Lithium", atomic_weight: 6.941, valence: 1, max_bonds: 1 },
    Element { atomic_number: 4, symbol: "Be", name: "Beryllium", atomic_weight: 9.012, valence: 2, max_bonds: 2 },
    Element { atomic_number: 5, symbol: "B", name: "Boron", atomic_weight: 10.81, valence: 3, max_bonds: 4 },
    Element { atomic_number: 6, symbol: "C", name: "Carbon", atomic_weight: 12.011, valence: 4, max_bonds: 4 },
    Element { atomic_number: 7, symbol: "N", name: "Nitrogen", atomic_weight: 14.007, valence: 3, max_bonds: 4 },
    Element { atomic_number: 8, symbol: "O", name: "Oxygen", atomic_weight: 15.999, valence: 2, max_bonds: 3 },
    Element { atomic_number: 9, symbol: "F", name: "Fluorine", atomic_weight: 18.998, valence: 1, max_bonds: 1 },
    Element { atomic_number: 10, symbol: "Ne", name: "Neon", atomic_weight: 20.180, valence: 0, max_bonds: 0 },
    Element { atomic_number: 11, symbol: "Na", name: "Sodium", atomic_weight: 22.990, valence: 1, max_bonds: 1 },
    Element { atomic_number: 12, symbol: "Mg", name: "Magnesium", atomic_weight: 24.305, valence: 2, max_bonds: 2 },
    Element { atomic_number: 13, symbol: "Al", name: "Aluminum", atomic_weight: 26.982, valence: 3, max_bonds: 4 },
    Element { atomic_number: 14, symbol: "Si", name: "Silicon", atomic_weight: 28.086, valence: 4, max_bonds: 4 },
    Element { atomic_number: 15, symbol: "P", name: "Phosphorus", atomic_weight: 30.974, valence: 3, max_bonds: 6 },
    Element { atomic_number: 16, symbol: "S", name: "Sulfur", atomic_weight: 32.06, valence: 2, max_bonds: 6 },
    Element { atomic_number: 17, symbol: "Cl", name: "Chlorine", atomic_weight: 35.45, valence: 1, max_bonds: 1 },
    Element { atomic_number: 18, symbol: "Ar", name: "Argon", atomic_weight: 39.948, valence: 0, max_bonds: 0 },
    Element { atomic_number: 19, symbol: "K", name: "Potassium", atomic_weight: 39.098, valence: 1, max_bonds: 1 },
    Element { atomic_number: 20, symbol: "Ca", name: "Calcium", atomic_weight: 40.078, valence: 2, max_bonds: 2 },
    Element { atomic_number: 21, symbol: "Sc", name: "Scandium", atomic_weight: 44.956, valence: 3, max_bonds: 6 },
    Element { atomic_number: 22, symbol: "Ti", name: "Titanium", atomic_weight: 47.867, valence: 4, max_bonds: 6 },
    Element { atomic_number: 23, symbol: "V", name: "Vanadium", atomic_weight: 50.942, valence: 5, max_bonds: 6 },
    Element { atomic_number: 24, symbol: "Cr", name: "Chromium", atomic_weight: 51.996, valence: 3, max_bonds: 6 },
    Element { atomic_number: 25, symbol: "Mn", name: "Manganese", atomic_weight: 54.938, valence: 2, max_bonds: 6 },
    Element { atomic_number: 26, symbol: "Fe", name: "Iron", atomic_weight: 55.845, valence: 3, max_bonds: 6 },
    Element { atomic_number: 27, symbol: "Co", name: "Cobalt", atomic_weight: 58.933, valence: 3, max_bonds: 6 },
    Element { atomic_number: 28, symbol: "Ni", name: "Nickel", atomic_weight: 58.693, valence: 2, max_bonds: 6 },
    Element { atomic_number: 29, symbol: "Cu", name: "Copper", atomic_weight: 63.546, valence: 2, max_bonds: 6 },
    Element { atomic_number: 30, symbol: "Zn", name: "Zinc", atomic_weight: 65.38, valence: 2, max_bonds: 4 },
    Element { atomic_number: 31, symbol: "Ga", name: "Gallium", atomic_weight: 69.723, valence: 3, max_bonds: 4 },
    Element { atomic_number: 32, symbol: "Ge", name: "Germanium", atomic_weight: 72.63, valence: 4, max_bonds: 4 },
    Element { atomic_number: 33, symbol: "As", name: "Arsenic", atomic_weight: 74.922, valence: 3, max_bonds: 5 },
    Element { atomic_number: 34, symbol: "Se", name: "Selenium", atomic_weight: 78.96, valence: 2, max_bonds: 6 },
    Element { atomic_number: 35, symbol: "Br", name: "Bromine", atomic_weight: 79.904, valence: 1, max_bonds: 1 },
    Element { atomic_number: 36, symbol: "Kr", name: "Krypton", atomic_weight: 83.798, valence: 0, max_bonds: 0 },
    Element { atomic_number: 37, symbol: "Rb", name: "Rubidium", atomic_weight: 85.468, valence: 1, max_bonds: 1 },
    Element { atomic_number: 38, symbol: "Sr", name: "Strontium", atomic_weight: 87.62, valence: 2, max_bonds: 2 },
    Element { atomic_number: 39, symbol: "Y", name: "Yttrium", atomic_weight: 88.906, valence: 3, max_bonds: 6 },
    Element { atomic_number: 40, symbol: "Zr", name: "Zirconium", atomic_weight: 91.224, valence: 4, max_bonds: 6 },
    Element { atomic_number: 41, symbol: "Nb", name: "Niobium", atomic_weight: 92.906, valence: 5, max_bonds: 6 },
    Element { atomic_number: 42, symbol: "Mo", name: "Molybdenum", atomic_weight: 95.95, valence: 6, max_bonds: 6 },
    Element { atomic_number: 43, symbol: "Tc", name: "Technetium", atomic_weight: 98.0, valence: 7, max_bonds: 7 },
    Element { atomic_number: 44, symbol: "Ru", name: "Ruthenium", atomic_weight: 101.07, valence: 4, max_bonds: 8 },
    Element { atomic_number: 45, symbol: "Rh", name: "Rhodium", atomic_weight: 102.906, valence: 3, max_bonds: 6 },
    Element { atomic_number: 46, symbol: "Pd", name: "Palladium", atomic_weight: 106.42, valence: 2, max_bonds: 6 },
    Element { atomic_number: 47, symbol: "Ag", name: "Silver", atomic_weight: 107.868, valence: 1, max_bonds: 4 },
    Element { atomic_number: 48, symbol: "Cd", name: "Cadmium", atomic_weight: 112.414, valence: 2, max_bonds: 4 },
    Element { atomic_number: 49, symbol: "In", name: "Indium", atomic_weight: 114.818, valence: 3, max_bonds: 4 },
    Element { atomic_number: 50, symbol: "Sn", name: "Tin", atomic_weight: 118.710, valence: 4, max_bonds: 4 },
    Element { atomic_number: 51, symbol: "Sb", name: "Antimony", atomic_weight: 121.760, valence: 3, max_bonds: 5 },
    Element { atomic_number: 52, symbol: "Te", name: "Tellurium", atomic_weight: 127.60, valence: 2, max_bonds: 6 },
    Element { atomic_number: 53, symbol: "I", name: "Iodine", atomic_weight: 126.904, valence: 1, max_bonds: 1 },
    Element { atomic_number: 54, symbol: "Xe", name: "Xenon", atomic_weight: 131.293, valence: 0, max_bonds: 0 },
];

/// Look up an element by its symbol (e.g. "C", "Fe").
pub fn element_by_symbol(symbol: &str) -> Option<&'static Element> {
    ELEMENTS.iter().find(|e| e.symbol == symbol)
}

/// Look up an element by its atomic number (1-based).
pub fn element_by_number(n: u8) -> Option<&'static Element> {
    if n >= 1 && n <= 54 {
        Some(&ELEMENTS[(n - 1) as usize])
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lookup_carbon_by_symbol() {
        let c = element_by_symbol("C").unwrap();
        assert_eq!(c.atomic_number, 6);
        assert_eq!(c.name, "Carbon");
        assert!((c.atomic_weight - 12.011).abs() < 0.001);
        assert_eq!(c.valence, 4);
    }

    #[test]
    fn lookup_nitrogen_by_number() {
        let n = element_by_number(7).unwrap();
        assert_eq!(n.symbol, "N");
        assert_eq!(n.name, "Nitrogen");
        assert_eq!(n.valence, 3);
    }

    #[test]
    fn unknown_returns_none() {
        assert!(element_by_symbol("Zz").is_none());
        assert!(element_by_number(0).is_none());
        assert!(element_by_number(55).is_none());
    }
}
