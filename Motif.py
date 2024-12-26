def identify_ORFs(input_file, output_file):
    sequences = {}
    current_header = None

    # Read the input File
    with open(input_file, 'r') as f:
        for line in f:
            line = line.strip()


            if line.startswith('>'):
                current_header = line
                sequences[current_header] = ''
            else:
                sequences[current_header] += line

    def scoreMotif(seq):
        base_idx = {'A': 0, 'T': 1, 'C': 2, 'G': 3}
        motif = [
            [0.5, 0.5, 0.5, 0.5, 0, 0, 0, 0, 0, 2, -99, -99, 0.5],  # A
            [0, 0, 0, 0, 0, 0, 0, 0, 0, -99, 2, -99, 0],  # T
            [0, 0, 0, 0, 0, 0, 0, 0, 0, -99, -99, -99, 0],  # C
            [0.5, 0.5, 0.5, 0.5, 0, 0, 0, 0, 0, 0.5, -99, 2, 0]  # G
        ]

        score = 0
        for i in range(13):
            nucleotide = seq[i]
            idx = base_idx[nucleotide]
            score += motif[idx][i]
        return score

    def scanSeq(sequence):
        start_positions = []
        orf_lengths = []
        orf_sequences = []

        for i in range(len(sequence) - 13):
            subseq = sequence[i:i + 13]
            motif_score = scoreMotif(subseq)

            if motif_score > 7.25:
                orf_start = i
                potential_orf = sequence[i:i + 60]

                for j in range(13, 60, 3):
                    codon = potential_orf[j:j + 3]
                    if codon in ['TAA', 'TAG', 'TGA']:
                        orf_end = i + j + 3
                        orf_seq = sequence[orf_start:orf_end]
                        orf_length = orf_end - orf_start

                        # Saving the ORF details
                        start_positions.append(orf_start)
                        orf_lengths.append(orf_length)
                        orf_sequences.append(orf_seq)
                        break

        return start_positions, orf_lengths, orf_sequences

    # Writing results to output file
    with open(output_file, 'w') as out_f:
        contig_number = 1

        for header, sequence in sequences.items():
            start_positions, orf_lengths, orf_sequences = scanSeq(sequence)

            for i in range(len(start_positions)):
                orf_header = f"> contig {contig_number}|Length {orf_lengths[i]}|at position {start_positions[i]}"
                out_f.write(orf_header + "\n")
                out_f.write(orf_sequences[i] + "\n")

            contig_number += 1


# Call the function with the input and output file names
if __name__ == "__main__":
    input_file = "spaceSeq.fa"  # Your input file name
    output_file = "output.fa"  # Specify your desired output file name
    identify_ORFs(input_file, output_file)
