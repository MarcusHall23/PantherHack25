# PantherHack25

A tool to detect genetic patterns and potential disease markers in DNA sequences.

---

## Overview

PantherHack25 was developed as part of the **Panther Hacks** hackathon, specifically for the **Healthcare Track**. This project aims to provide an accessible and efficient tool to analyze DNA sequences and identify genetic mutations linked to diseases.

The tool analyzes long DNA sequences (strings of A, T, C, and G) to find known mutation patterns associated with genetic conditions, helping advance personalized medicine and healthcare research.

---

## Features

- **Pattern Matching:**  
  Uses string search algorithms and regular expressions (e.g., for Indels) to detect mutations in DNA sequences.

- **Mutation Types Supported:**  
  - SNPs (Single Nucleotide Polymorphisms)  
  - Indels (Insertions/Deletions with replacement sequences)

- **Input:**  
  - DNA sequence strings  
  - FASTA file uploads
  - Gene input for targeted mutation analysis

- **Output:**  
  - Mutation name, type, gene symbol, position  
  - Clinical significance and associated conditions  
  - Structured findings for downstream analysis

---

## Usage Example

### Running the tool

To analyze a patient's DNA sequence and check for known mutations, run the `main.py` script. The script will ask for the following inputs:

1. **Gene Name**: Enter the gene you want to check (e.g., NF1).
2. **FASTA File Path**: Provide the path to the patient’s DNA sequence in FASTA format.
3. **Variant CSV File**: The script uses a hardcoded path to a CSV file with known pathogenic variants. The variants file can be found in `outputs/pathogenic_variants_sorted.csv`.

Example logic for loading the data:

```python
def main():
    chromosome_input = input("Enter the gene being presented: ")
    fasta_path = input("Enter Patient Sequence File: ")
    variants_csv_path = "C:\\Users\\rogue\\PantherHack2025\\PantherHack25\\outputs\\pathogenic_variants_sorted.csv"

    print("Loading patient DNA...")
    # dna_sequence = load_fasta_sequence(fasta_path)

    print("Loading known pathogenic variants...")
    variants = load_variants(variants_csv_path)

    print("Analyzing for mutations...")
    findings = analyze_variants(fasta_path, variants, chromosome_input)
    if findings:
        print(f"\nMutations found for Chromosome {chromosome_input}:")
        for finding in findings:
            print(finding)
    else:
        print(f"\nNo mutations found for Chromosome {chromosome_input}.")
