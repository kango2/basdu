#Author: Hardip Patel, Kirat Alreja

import pandas as pd
import sys
print("Number of arguments:", len(sys.argv))
print("Arguments:", sys.argv)

trf_input = sys.argv[1]
trf = pd.read_csv(trf_input)
fasta = sys.argv[2]
output_csv = sys.argv[3]
copies = int(sys.argv[4])
permatch = int(sys.argv[5])

def parse_fasta_to_df(fasta_file_path):
    sequences = {}
    current_seq_name = ""
    current_seq = ""
    
    with open(fasta_file_path, 'r') as fasta_file:
        for line in fasta_file:
            line = line.strip()
            if line.startswith(">"):
                if current_seq_name != "":
                    sequences[current_seq_name] = len(current_seq)
                    current_seq = ""
                current_seq_name = line[1:]  # Remove the '>' symbol
            else:
                current_seq += line
        
        # Add the last sequence to the dictionary
        if current_seq_name != "":
            sequences[current_seq_name] = len(current_seq)
    
    # Convert the dictionary to a DataFrame
    df = pd.DataFrame(list(sequences.items()), columns=['Sequence_ID', 'Length'])
    return df

def get_reverse_complement(seq):
    """Return the reverse complement of a DNA sequence."""
    complement = {
        'A': 'T',
        'T': 'A',
        'G': 'C',
        'C': 'G'
    }
    return "".join(complement[base] for base in reversed(seq))

# Base telomeric repeat
base_repeat = "TTAGGG"

# Generate possible variations of the telomeric repeat
variations = [base_repeat[i:] + base_repeat[:i] for i in range(len(base_repeat))]

# Get reverse complements of each variation
reverse_complements = [get_reverse_complement(variation) for variation in variations]

# Combine the variations and their reverse complements
telomeric_variations = set(variations + reverse_complements)

seq_lens = parse_fasta_to_df(fasta)
result_df = trf.merge(seq_lens, on='Sequence_ID', how='left')
result_df['Relative Start'] = (result_df['Start']/result_df['Length']).round(2)
result_df['Relative End'] = (result_df['End']/result_df['Length']).round(2)
result_df = result_df[result_df['cons_seq'].isin(telomeric_variations)]
result_df = result_df[result_df['copies']>=copies]
result_df = result_df[result_df['perc_match']>=permatch]
result_df.to_csv(output_csv, index=False)
