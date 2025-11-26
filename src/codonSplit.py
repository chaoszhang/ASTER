import sys

def split_codon_positions(input_fasta):
    with open(input_fasta, 'r') as infile, \
         open(input_fasta + ".pos1", 'w') as out1, \
         open(input_fasta + ".pos2", 'w') as out2, \
         open(input_fasta + ".pos3", 'w') as out3:

        sequence_id = None
        sequences = {}

        # Read the input FASTA file
        for line in infile:
            line = line.strip()
            if line.startswith('>'):
                sequence_id = line[1:]
                sequences[sequence_id] = ''
            else:
                sequences[sequence_id] += line

        # Split sequences into codon positions
        for seq_id, sequence in sequences.items():
            codon1 = []
            codon2 = []
            codon3 = []

            for i in range(0, len(sequence), 3):
                if i < len(sequence):
                    codon1.append(sequence[i])
                if i + 1 < len(sequence):
                    codon2.append(sequence[i + 1])
                if i + 2 < len(sequence):
                    codon3.append(sequence[i + 2])

            codon1_seq = ''.join(codon1)
            codon2_seq = ''.join(codon2)
            codon3_seq = ''.join(codon3)

            # Write to the output files
            out1.write(f">{seq_id}_pos1\n{codon1_seq}\n")
            out2.write(f">{seq_id}_pos2\n{codon2_seq}\n")
            out3.write(f">{seq_id}_pos3\n{codon3_seq}\n")

# Run the function
split_codon_positions(sys.argv[1])
