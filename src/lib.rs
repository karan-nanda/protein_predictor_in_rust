use std::collections::HashMap;

// Defining the amino acid parameters as static hashmaps

pub fn alpha_param() -> HashMap<char, f64> {
    let mut map = HashMap::new();
    map.insert('A', 1.45); // Alanine
    map.insert('R', 0.79); // Arginine
    map.insert('N', 0.73); // Asparagine
    map.insert('D', 0.98); // Aspartic acid
    map.insert('C', 0.77); // Cysteine
    map.insert('Q', 1.17); // Glutamine
    map.insert('E', 1.53); // Glutamic acid
    map.insert('G', 0.53); // Glycine
    map.insert('H', 1.24); // Histidine
    map.insert('I', 1.00); // Isoleucine
    map.insert('L', 1.34); // Leucine
    map.insert('K', 1.07); // Lysine
    map.insert('M', 1.2);  // Methionine
    map.insert('F', 1.12); // Phenylalanine
    map.insert('P', 0.59); // Proline
    map.insert('S', 0.79); // Serine
    map.insert('T', 0.82); // Threonine
    map.insert('W', 1.14); // Tryptophan
    map.insert('Y', 0.61); // Tyrosine
    map.insert('V', 1.14); // Valine
    map.insert('B', 0.86); // Average of Aspartate and Asparagine
    map.insert('Z', 1.35); // Average of Glutamine and Glutamate
    map.insert('X', 1.00); // Average for all amino acids
    map
}

pub fn alpha_inner_param() -> HashMap<char, f64> {
    let mut map = HashMap::new();
    map.insert('A', 1.59); // Alanine
    map.insert('R', 0.67); // Arginine
    map.insert('N', 0.53); // Asparagine
    map.insert('D', 0.53); // Aspartic acid
    map.insert('C', 0.33); // Cysteine
    map.insert('Q', 0.98); // Glutamine
    map.insert('E', 1.45); // Glutamic acid
    map.insert('G', 0.53); // Glycine
    map.insert('H', 0.87); // Histidine
    map.insert('I', 1.22); // Isoleucine
    map.insert('L', 1.91); // Leucine
    map.insert('K', 1.13); // Lysine
    map.insert('M', 1.25); // Methionine
    map.insert('F', 1.14); // Phenylalanine
    map.insert('P', 0.0);  // Proline
    map.insert('S', 0.7);  // Serine
    map.insert('T', 0.75); // Threonine
    map.insert('W', 1.33); // Tryptophan
    map.insert('Y', 0.58); // Tyrosine
    map.insert('V', 1.42); // Valine
    map.insert('B', 0.53); // Average of Aspartate and Asparagine
    map.insert('Z', 1.21); // Average of Glutamine and Glutamate
    map.insert('X', 0.95); // Average for all amino acids
    map
}

pub fn beta_param() -> HashMap<char, f64> {
    let mut map = HashMap::new();
    map.insert('A', 0.97); // Alanine
    map.insert('R', 0.9);  // Arginine
    map.insert('N', 0.65); // Asparagine
    map.insert('D', 0.8);  // Aspartic acid
    map.insert('C', 1.3);  // Cysteine
    map.insert('Q', 1.23); // Glutamine
    map.insert('E', 0.26); // Glutamic acid
    map.insert('G', 0.81); // Glycine
    map.insert('H', 0.71); // Histidine
    map.insert('I', 1.6);  // Isoleucine
    map.insert('L', 1.22); // Leucine
    map.insert('K', 0.74); // Lysine
    map.insert('M', 1.67); // Methionine
    map.insert('F', 1.28); // Phenylalanine
    map.insert('P', 0.62); // Proline
    map.insert('S', 0.72); // Serine
    map.insert('T', 1.2);  // Threonine
    map.insert('W', 1.19); // Tryptophan
    map.insert('Y', 1.29); // Tyrosine
    map.insert('V', 1.65); // Valine
    map.insert('B', 0.73); // Average of Aspartate and Asparagine
    map.insert('Z', 0.75); // Average of Glutamine and Glutamate
    map.insert('X', 1.04); // Average for all amino acids
    map
}

pub fn coil_param() -> HashMap<char, f64> {
    let mut map = HashMap::new();
    map.insert('A', 0.66); // Alanine
    map.insert('R', 1.2);  // Arginine
    map.insert('N', 1.33); // Asparagine
    map.insert('D', 1.09); // Aspartic acid
    map.insert('C', 1.07); // Cysteine
    map.insert('Q', 0.79); // Glutamine
    map.insert('E', 0.87); // Glutamic acid
    map.insert('G', 1.42); // Glycine
    map.insert('H', 0.92); // Histidine
    map.insert('I', 0.78); // Isoleucine
    map.insert('L', 0.66); // Leucine
    map.insert('K', 1.05); // Lysine
    map.insert('M', 0.61); // Methionine
    map.insert('F', 0.81); // Phenylalanine
    map.insert('P', 1.45); // Proline
    map.insert('S', 1.27); // Serine
    map.insert('T', 1.05); // Threonine
    map.insert('W', 0.82); // Tryptophan
    map.insert('Y', 1.19); // Tyrosine
    map.insert('V', 0.66); // Valine
    map.insert('B', 1.21); // Average of Aspartate and Asparagine
    map.insert('Z', 0.83); // Average of Glutamine and Glutamate
    map.insert('X', 0.99); // Average for all amino acids
    map
}
