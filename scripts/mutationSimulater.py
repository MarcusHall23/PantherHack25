import csv
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def parse_mutation_name(mutation_name, mutation_type):
    """
    Parses the mutation name and returns the start and end coordinates.
    Additionally, for Indels, it will return the inserted sequence.
    """

    # Case 1: Single Nucleotide Variant (SNV)
    if mutation_type == 'single nucleotide variant':
        # Case 1.1: Standard SNV (e.g., c.694C>T, c.1234A>G)
        match = re.search(r'c\.(\d+)([A-Za-z])>([A-Za-z])', mutation_name)
        if match:
            start = int(match.group(1))  # Start position (e.g., 694 or 1234)
            ref_allele = match.group(2)  # Reference allele (e.g., C or A)
            alt_allele = match.group(3)  # Alternate allele (e.g., T or G)
            return start, start, alt_allele, ref_allele

        # Case 1.2: Splice site mutation (e.g., c.700+4G>T, c.581-176A>T)
        match = re.search(r'c\.(\d+)([+-]\d+)([A-Za-z])>([A-Za-z])', mutation_name)
        if match:
            exon_position = int(match.group(1))  # Exon base
            splice_offset = int(match.group(2))  # +4 or -176 (converted to int)
            final_position = exon_position + splice_offset  # Compute final position
            ref_allele = match.group(3)
            alt_allele = match.group(4)
            return final_position, final_position, alt_allele, ref_allele
        else:
            return -1, -1, -1 , -1

    # Case 2: Insertion (Indel)
    elif mutation_type == 'Indel':
        # Regular expression for Indels (e.g., c.80_83delinsTGCTGTAAACTGTAACTGTAAA)
        match = re.search(r'c\.(\d+)_?(\d+)(delins[A-Za-z]+)', mutation_name)
        if match:
            start = int(match.group(1))  # Start position (e.g., 80)
            end = int(match.group(2))    # End position (e.g., 83)
            inserted_sequence = match.group(3)[5:]  # Get the inserted bases after 'delins'
            return start, end, inserted_sequence
        else:
            return -1, -1, -1

    # Case 3: Deletion (e.g., c.1234_1235del)
    elif mutation_type == 'Deletion':
            match = re.search(r'c\.(\d+[+-]?\d*)_(\d+[+-]?\d*)del', mutation_name)
            if match:
                start = match.group(1)  # '361-5'
                end = match.group(2)    # '361-1'

                # Now, remove '+' or '-' for easier processing
                start_clean = int(re.sub(r'[+-]', '', start))
                end_clean = int(re.sub(r'[+-]', '', end))

                return start_clean, end_clean  # No sequence for deletion
            else:
                return -1, -1

    else:
        raise ValueError(f"Mutation type {mutation_type} is not recognized.")



def simulate_mutation(reference_sequence, mutation_name, mutation_type, mutation_position):
    """
    Simulate a mutation on the reference sequence based on mutation type.
    """
    # start, end, alt_allele, ref_allele = parse_mutation_name(mutation_name, mutation_type)

    # if start == -1 or end == -1:
    #     # Invalid mutation
    #     return reference_sequence

    # Case 1: SNV (Single Nucleotide Variant)
    if mutation_type == 'single nucleotide variant':
        start, end, alt_allele, ref_allele = parse_mutation_name(mutation_name, mutation_type)
        # if start == -1 or end == -1:
        #     # Invalid mutation
        #     return reference_sequence
        # Replace the reference allele with the alternate allele
        start += mutation_position
        end += mutation_position
        mutated_sequence = reference_sequence[:start] + alt_allele + reference_sequence[start + 1:]
        return mutated_sequence
    
    # Case 2: Indel (Insertion or Deletion)
    elif mutation_type == 'Indel':
        start, end, inserted_sequence = parse_mutation_name(mutation_name, mutation_type)
        # if start == -1 or end == -1:
        #     # Invalid mutation
        #     return reference_sequence
        # If it's an insertion
        start += mutation_position
        end += mutation_position
        if inserted_sequence is not None:
            mutated_sequence = reference_sequence[:start] + inserted_sequence + reference_sequence[start:]
            return mutated_sequence
            
        # If it's a deletion, we simply remove the base(s)
        else:
            mutated_sequence = reference_sequence[:start] + reference_sequence[end:]
            return mutated_sequence

    # Case 3: Deletion
    elif mutation_type == 'Deletion':
        start, end = parse_mutation_name(mutation_name, mutation_type)
        if start == -1 or end == -1:
            # Invalid mutation
            return reference_sequence
        start += mutation_position
        end += mutation_position
        # Deletion: Remove the region between start and end positions
        mutated_sequence = reference_sequence[:start] + reference_sequence[end + 1:]
        return mutated_sequence
    else:
        return reference_sequence

    

def apply_mutations(reference_sequence, mutations, chromo):
    """
    Apply a list of mutations to the reference sequence.
    """
    mutated_sequence = reference_sequence

    chromo = int(chromo)
    
    for mutation in mutations:
        if int(mutation["Chromosome"]) == chromo:
            print("hello world")
            mutation_name = mutation['MutationName'].strip()
            mutation_type = mutation['Type'].strip()
            mutation_position = mutation['Start']
            mutated_sequence = simulate_mutation(mutated_sequence, mutation_name, mutation_type, mutation_position)
        else:
            continue
    
    return mutated_sequence

# def load_variants(csv_path):
#     variants = []
#     with open(csv_path, 'r') as file:
#         reader = csv.DictReader(file)
#         # Clean up column names by stripping leading/trailing spaces
#         reader.fieldnames = [field.strip() for field in reader.fieldnames]
        
#         for row in reader:
#             variants.append(row)
#     return variants

def load_fasta_sequence(fasta_path):
    """Loads the whole FASTA sequence into one string."""
    sequence_lines = []
    with open(fasta_path, 'r') as file:
        for line in file:
            if not line.startswith('>'):
                sequence_lines.append(line.strip())
    return ''.join(sequence_lines)

def save_sequence(mutated_sequence, output_file, record_id):
    print(f"Saving to {output_file}...")
    seq_record = SeqRecord(Seq(mutated_sequence), id="chr"+record_id, description="Mutated Sequence")

    with open(output_file, "w") as file_handle:
        SeqIO.write(seq_record, file_handle, "fasta")
    print(f"Sequence saved to {output_file}")

def main():
    chromosome_input = input("Enter the chromosome being presented: ")
    fasta_path_Reference = input("Enter Reference Chromosome Sequence File: ")
    variants_csv_path = "C:\\Users\\rogue\\PantherHack2025\\PantherHack25\\outputs\\pathogenic_variants_sorted.csv"
    output_file_name = input("Output File Name: (MUTATEDchr???.fa)")
    mutations = [
    {"Gene": "RB1", "Type": "single nucleotide variant", "MutationName": "NM_000321.3(RB1):c.784C>T (p.Arg262Trp)", "Chromosome": 13, "Start": 48362880},
    {"Gene": "BRCA2", "Type": "Deletion", "MutationName": "NM_000059.4(BRCA2):c.7447del (p.Ser2483fs)", "Chromosome": 13, "Start": 32356438},
    {"Gene": "BRCA2", "Type": "Indel", "MutationName": "NM_000059.4(BRCA2):c.3645_3647delinsTA (p.Phe1216fs)", "Chromosome": 13, "Start": 32912137},
    {"Gene": "KMT2C", "Type": "single nucleotide variant", "MutationName": "NM_170606.3(KMT2C):c.7042C>T (p.Gln2348Ter)", "Chromosome": 7, "Start": 152180818},
    # {"Gene": "KMT2C", "Type": "Deletion", "MutationName": "NM_170606.3(KMT2C):c.1238_1254del (p.Tyr413fs)", "Chromosome": 7, "Start": 152263061},
    # {"Gene": "KRIT1", "Type": "Indel", "MutationName": "NM_194454.3(KRIT1):c.879_880delinsTT (p.Arg294Ter)", "Chromosome": 7, "Start": 91863872},
    {"Gene": "SDCCAG8", "Type": "Deletion", "MutationName": "NM_006642.5(SDCCAG8):c.1420del (p.Glu474fs)", "Chromosome": 1, "Start": 243507579},
    {"Gene": "SDCCAG8", "Type": "single nucleotide variant", "MutationName": "NM_006642.5(SDCCAG8):c.740+356C>T", "Chromosome": 1, "Start": 243468435},
    {"Gene": "UROD", "Type": "single nucleotide variant", "MutationName": "NM_000374.5(UROD):c.874C>G (p.Arg292Gly)", "Chromosome": 1, "Start": 45480507},
    {"Gene": "P3H1", "Type": "Indel", "MutationName": "NM_022356.4(P3H1):c.1365_1366delinsC (p.Glu455fs)", "Chromosome": 1, "Start": 42752644}
]
    reference_sequence = load_fasta_sequence(fasta_path_Reference)
    mutated_sequence = apply_mutations(reference_sequence, mutations, chromosome_input)
    save_sequence(mutated_sequence, output_file_name, chromosome_input)

if __name__ == "__main__":
    main()