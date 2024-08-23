import sys
import csv

def read_fasta(file_path):
    with open(file_path, 'r') as f:
        sequences = {}
        current_header = None
        
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                current_header = line[1:]
                sequences[current_header] = []
            else:
                sequences[current_header].append(line)
                
        for key, value in sequences.items():
            sequences[key] = ''.join(value)
            
        return sequences

def gc_content(seq, window=10000):  # changed to 10 Mb
    gc_counts = []
    
    for i in range(0, len(seq) - window + 1, window):
        subseq = seq[i:i+window]
        gc_count = subseq.count('G') + subseq.count('C')
        gc_counts.append(gc_count)
    
    return gc_counts

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Please provide the path to the FASTA file as an argument.")
        sys.exit(1)

    fasta_path = sys.argv[1]
    
    # Get CSV output path
    csv_path = fasta_path.rsplit('.', 1)[0] + '.csv'
    
    sequences = read_fasta(fasta_path)
    
    with open(csv_path, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(['Sequence Header', 'Position Start', 'Position End', 'GC Count'])

        for header, seq in sequences.items():
            gc_counts = gc_content(seq)
            
            for i, count in enumerate(gc_counts):
                csvwriter.writerow([header, i*10000, (i+1)*10000 - 1, count])  # changed to 10 Mb
                
    print(f"Results written to {csv_path}")
