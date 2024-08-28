# Import necessary libraries
import pandas as pd
import sys

# Accept command-line arguments
trf_input = sys.argv[1]  # Input CSV file path containing TRF data
fasta = sys.argv[2]      # Input FASTA file path containing sequence data
output_csv = sys.argv[3] # Output CSV file path for the filtered results
copies = int(sys.argv[4]) # Minimum number of repeat copies to filter
permatch = int(sys.argv[5]) # Minimum percentage match to filter

# Read TRF data from CSV file
trf = pd.read_csv(trf_input)

def parse_fasta_to_df(fasta_file_path):
    """
    Parse a FASTA file and convert it to a DataFrame containing sequence IDs and their lengths.

    Args:
    fasta_file_path (str): The path to the FASTA file.

    Returns:
    pd.DataFrame: A DataFrame with columns 'Sequence_ID' and 'Length'.
    """
    sequences = {}  # Dictionary to store sequence ID and length
    current_seq_name = ""  # Temporary variable to store current sequence ID
    current_seq = ""  # Temporary variable to store current sequence

    # Open and read the FASTA file
    with open(fasta_file_path, 'r') as fasta_file:
        for line in fasta_file:
            line = line.strip()
            if line.startswith(">"):  # Detect sequence ID lines
                if current_seq_name != "":
                    # Save the previous sequence and its length
                    sequences[current_seq_name] = len(current_seq)
                    current_seq = ""  # Reset the sequence for the next one
                current_seq_name = line[1:]  # Store the new sequence ID (without '>')
            else:
                # Append sequence data
                current_seq += line
        
        # Add the last sequence to the dictionary
        if current_seq_name != "":
            sequences[current_seq_name] = len(current_seq)
    
    # Convert the dictionary to a DataFrame
    df = pd.DataFrame(list(sequences.items()), columns=['Sequence_ID', 'Length'])
    return df

def get_reverse_complement(seq):
    """
    Return the reverse complement of a DNA sequence.

    Args:
    seq (str): A DNA sequence.

    Returns:
    str: The reverse complement of the input DNA sequence.
    """
    complement = {
        'A': 'T',
        'T': 'A',
        'G': 'C',
        'C': 'G'
    }
    # Compute the reverse complement
    return "".join(complement[base] for base in reversed(seq))

# Base telomeric repeat sequence
base_repeat = "TTAGGG"

# Generate all cyclic variations of the base telomeric repeat
variations = [base_repeat[i:] + base_repeat[:i] for i in range(len(base_repeat))]

# Get reverse complements of each variation
reverse_complements = [get_reverse_complement(variation) for variation in variations]

# Combine the variations and their reverse complements into a set (to remove duplicates)
telomeric_variations = set(variations + reverse_complements)

# Parse the FASTA file to get a DataFrame of sequence IDs and their lengths
seq_lens = parse_fasta_to_df(fasta)

# Merge TRF data with sequence lengths based on 'Sequence_ID'
result_df = trf.merge(seq_lens, on='Sequence_ID', how='left')

# Calculate relative start and end positions of repeats within the sequences
result_df['Relative Start'] = (result_df['Start'] / result_df['Length']).round(2)
result_df['Relative End'] = (result_df['End'] / result_df['Length']).round(2)

# Filter the DataFrame for rows where the consensus sequence matches telomeric variations
result_df = result_df[result_df['cons_seq'].isin(telomeric_variations)]

# Further filter for rows with a minimum number of repeat copies
result_df = result_df[result_df['copies'] >= copies]

# Further filter for rows with a minimum percentage match
result_df = result_df[result_df['perc_match'] >= permatch]

# Save the filtered DataFrame to a CSV file
result_df.to_csv(output_csv, index=False)
