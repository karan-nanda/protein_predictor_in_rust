use std::fs::File;
use std::io::{BufRead, BufReader};

mod structure_predictor;
use structure_predictor::{find_nucleation_regions, extend_regions, print_regions}; 


#[derive(Debug)]
pub struct Protein{
    pub header : String,
    pub sequence: String,
}

pub fn parse_fasta(file_path: &str) -> Result<Vec<Protein>, std::io::Error> {
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);

    let mut proteins = Vec::new();
    let mut current_protein = Protein {
        header : String::new(),
        sequence : String::new(),
    };


    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            if !current_protein.header.is_empty() {
                proteins.push (current_protein);
                current_protein = Protein {
                    header : String::new(),
                    sequence : String::new(),
                };
            }
            current_protein.header = line[1..].to_string();
        } else {
            current_protein.sequence.push_str(&line);
        }
    }
    
    if !current_protein.header.is_empty() {
        proteins.push(current_protein);
    }

    Ok(proteins)
}

fn main() {
    let fasta_file = r"D:\VSCode\lung_cancer_FASTA\combined_fasta.fasta";
    match parse_fasta(fasta_file) {
        Ok(proteins) => {
            for protein in proteins {
                println!("Predicting structure for: {:?}", protein.header);
                // Use the protein's sequence for prediction
                let beta_regions = find_nucleation_regions(&protein.sequence, "beta", 3, 5, 0.5);
                let beta_regions_extended = extend_regions(&protein.sequence, beta_regions, "beta");

                // Print extended beta regions
                print_regions(&protein.sequence, &beta_regions_extended, "green");
            }
        }
        Err(e) => {
            println!("Error parsing FASTA file: {:?}", e);
        }
    }
}