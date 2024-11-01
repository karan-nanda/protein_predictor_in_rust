//use std::collections::HashMap;
use std::thread::sleep;
use std::time::Duration;
use colored::*;

use protein_predictor_in_rust::{alpha_param, alpha_inner_param, beta_param, coil_param};

fn main() {
    let _alpha = alpha_param();
    let _alpha_inner = alpha_inner_param();
    let _beta = beta_param();
    let _coil = coil_param();

    let sequence = "MKAILVGGVVFAGKNTDIGAVALKKLTFGKTQAKKAERADLIAYIKGEIRKASEGDLTFLGGPMQKIVGKADLVATELIEYKMALRKGALDVEIGLKKLGVEQMKALFVGKTNDRIAAFGNQVTIGTKKVEAYIGMKGADIRAKTEL";

    // parameters for beta_nucleation
    let beta_regions = find_nucleation_regions(sequence, "beta", 3,5,0.5);
    let beta_regions_extended = extend_regions(sequence, beta_regions, "beta");

    print_regions(sequence, &beta_regions_extended, "blue");

}


pub fn find_nucleation_regions(sequence: &str, conformation: &str, window_threshold: usize, window_size: usize, min_param: f64) -> Vec<(usize, usize)> {
    let mut nucleation_regions = Vec::new();
    let seq_len = sequence.len();

    for i in 0..(seq_len - window_size + 1) {
        let window = &sequence[i..i + window_size];
        if is_nucleation_region(window, conformation, window_threshold, min_param) {
            nucleation_regions.push((i, i + window_size));
            print_regions(sequence, &[(i, i + window_size)], "red"); // Color nucleation regions red
        } else {
            print_regions(sequence, &[(i, i + window_size)], "blue"); // Color non-nucleation regions blue
        }
        sleep(Duration::from_millis(100)); // Wait time for visualization
    }
    nucleation_regions
}


pub fn extend_regions(seq: &str, regions: Vec<(usize, usize)>, conformation: &str) -> Vec<(usize, usize)> {
    let mut extended_regions = Vec::new();

    for reg in regions {
        let mut n_shift = 0;
        let mut c_shift = 0;

        // Helix extension: grow in multiples of ~3 residues
        if conformation == "alpha" {
            while reg.0 > n_shift * 3 && average_param(conformation, &seq[reg.0 - n_shift * 3..reg.0]) >= 1.0 {
                n_shift += 1;
            }
            while reg.1 + c_shift * 3 < seq.len() && average_param(conformation, &seq[reg.1..reg.1 + c_shift * 3]) >= 1.0 {
                c_shift += 1;
            }
        }
        // Strand extension: extend residue by residue if it remains beta-favorable
        else if conformation == "beta" {
            let mut shift = (reg.0, reg.0 + 4);

            // Moving towards the N-Terminus
            while shift.0 > 0 && average_param(conformation, &seq[shift.0..shift.1]) >= 0.7 {
                n_shift += 1;
                shift = (shift.0 - 1, shift.1); 
            }

            // Moving towards the C-Terminus
            shift = (reg.1, reg.1 + 4);
            while shift.1 < seq.len() && average_param(conformation, &seq[shift.0..shift.1]) >= 0.7 {
                c_shift += 1;
                shift = (shift.0, shift.1 + 1);
            }
        }

        extended_regions.push((reg.0 - n_shift * 3, reg.1 + c_shift * 3));
    }

    merge_overlapping_indexes(extended_regions)
}





pub fn is_nucleation_region(window: &str, conformation: &str, window_threshold: usize, min_param: f64) -> bool {
    let qualified_count = window.chars().filter(|&res| conformation_parameter(conformation, res) > min_param).count();
    println!("Window: {}, Qualified Count: {}", window, qualified_count); // Debug print
    qualified_count >= window_threshold && conformation_has_highest_average(conformation, window)
}

pub fn conformation_has_highest_average(conformation: &str, sequence: &str) -> bool {
    let average = average_param(conformation, sequence);
    let alpha_average = average_param("alpha", sequence);
    let a_inner_average = average_param("alpha_inner", sequence);
    let beta_average = average_param("beta", sequence);
    let coil_average = average_param("coil", sequence);

    println!("Averages: {} - Alpha: {}, Alpha Inner: {}, Beta: {}, Coil: {}", average, alpha_average, a_inner_average, beta_average, coil_average); // Debug print

    match conformation {
        "alpha" => average > a_inner_average && average > beta_average && average > coil_average,
        "alpha_inner" => average > alpha_average && average > beta_average && average > coil_average,
        "beta" => average > alpha_average && average > a_inner_average && average > coil_average,
        "coil" => average > alpha_average && average > a_inner_average && average > beta_average,
        _ => false,
    }
}

pub fn average_param(conformation: &str, window: &str) -> f64 {
    let total: f64 = window.chars().map(|res| conformation_parameter(conformation, res)).sum();
    total / window.len() as f64
}


pub fn conformation_parameter(conformation: &str, res: char) -> f64 {
    match conformation {
        "alpha" => *protein_predictor_in_rust::alpha_param().get(&res).unwrap_or(&0.0),
        "alpha_inner" => *protein_predictor_in_rust::alpha_inner_param().get(&res).unwrap_or(&0.0),
        "beta" => *protein_predictor_in_rust::beta_param().get(&res).unwrap_or(&0.0),
        "coil" => *protein_predictor_in_rust::coil_param().get(&res).unwrap_or(&0.0),
        _ => 0.0,
    }
}

pub fn print_regions(seq: &str, regions: &[(usize, usize)], color: &str) {
    for &(start, end) in regions {
        let region_str = &seq[start..end]; // Keep this as &str
        let colored_output = match color {
            "red" => region_str.red(),
            "blue" => region_str.blue(),
            "green" => region_str.green(),
            "yellow" => region_str.yellow(),
            _ => region_str.normal(),
        };
        println!("Region ({}:{}): {}", start, end, colored_output);
    }
}

pub fn merge_overlapping_indexes(indexes: Vec<(usize, usize)>) -> Vec<(usize, usize)> {
    let mut merged_indexes = Vec::new();
    if indexes.len() < 2 {
        return indexes;
    }

    let mut i = 0;
    while i < indexes.len() - 1 {
        let current = indexes[i];
        let next = indexes[i + 1];
        if !merged_indexes.is_empty() && try_merging(merged_indexes.last().unwrap(), &current).is_some() {
            merged_indexes.last_mut().map(|last| *last = try_merging(last, &current).unwrap());
        } else if try_merging(&current, &next).is_some() {
            merged_indexes.push(try_merging(&current, &next).unwrap());
            i += 1; // iterate to the next index for merging
        } else {
            merged_indexes.push(current);
        }
        i += 1;
    }
    merged_indexes
}

pub fn try_merging(i: &(usize, usize), j: &(usize, usize)) -> Option<(usize, usize)> {
    if i.1 >= j.0 {
        Some((i.0, j.1))
    } else {
        None
    }
}