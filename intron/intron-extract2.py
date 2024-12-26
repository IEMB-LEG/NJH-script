import argparse
from Bio import SeqIO
from Bio.Seq import Seq

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Extract intron sequences from a genome based on a GFF3 file, reverse complementing sequences from negative strands.")
    parser.add_argument("genome_file", type=str, help="Path to the genome sequence file in FASTA format.")
    parser.add_argument("gff_file", type=str, help="Path to the intron GFF3 file.")
    parser.add_argument("output_file", type=str, help="Path to the output FASTA file for extracted intron sequences.")
    return parser.parse_args()

def extract_introns(genome_file, gff_file, output_file):
    """Extract intron sequences and reverse complement if on negative strand."""
    genome = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))
    
    with open(gff_file, 'r') as gff, open(output_file, 'w') as out_fasta:
        for line in gff:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            chrom, feature, start, end, strand = parts[0], parts[2], int(parts[3]), int(parts[4]), parts[6]
            if feature.lower() == 'intron':
                attributes = parts[8]
                intron_id = [attr for attr in attributes.split(';') if 'ID=' in attr][0].split('ID=')[1]
                intron_seq = genome[chrom].seq[start-1:end]
                if strand == '-':
                    intron_seq = intron_seq.reverse_complement()
                out_fasta.write(f">{intron_id}\n{intron_seq}\n")

def main():
    args = parse_args()
    extract_introns(args.genome_file, args.gff_file, args.output_file)
    print(f"Intron sequences have been extracted to {args.output_file}")

if __name__ == "__main__":
    main()

