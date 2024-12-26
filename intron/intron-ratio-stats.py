from Bio import SeqIO
import argparse
import pandas as pd

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Summarize intron statistics for each gene based on genome annotation and sequence.")
    parser.add_argument("genome_file", type=str, help="Path to the genome sequence file in FASTA format.")
    parser.add_argument("gff_file", type=str, help="Path to the GFF3 genome annotation file.")
    parser.add_argument("output_file", type=str, help="Path to the output CSV file for gene and intron statistics.")
    return parser.parse_args()

def parse_gff(gff_file):
    """Parse GFF3 file to get exons and calculate introns and gene lengths."""
    genes = {}
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            chrom, feature, start, end, attributes = parts[0], parts[2], int(parts[3]), int(parts[4]), parts[8]
            if feature.lower() == 'exon':
                parent_id = [attr for attr in attributes.split(';') if 'Parent=' in attr][0].split('Parent=')[1]
                if parent_id not in genes:
                    genes[parent_id] = {'chrom': chrom, 'exons': [], 'gene_length': 0, 'intron_count': 0, 'total_intron_length': 0}
                genes[parent_id]['exons'].append((start, end))
    
    for gene_id, gene_info in genes.items():
        gene_info['exons'].sort(key=lambda x: x[0])
        gene_info['gene_length'] = gene_info['exons'][-1][1] - gene_info['exons'][0][0] + 1
        intron_lengths = [gene_info['exons'][i+1][0] - gene_info['exons'][i][1] - 1 for i in range(len(gene_info['exons']) - 1)]
        gene_info['intron_count'] = len(intron_lengths)
        gene_info['total_intron_length'] = sum(intron_lengths)
        if gene_info['gene_length'] > 0:
            gene_info['intron_length_ratio'] = gene_info['total_intron_length'] / gene_info['gene_length']
        else:
            gene_info['intron_length_ratio'] = 0
    
    return genes

def write_to_csv(genes, output_file):
    """Write gene statistics to CSV."""
    with open(output_file, 'w') as out_file:
        out_file.write("GeneID,Chromosome,IntronCount,TotalIntronLength,GeneLength,IntronLengthRatio\n")
        for gene_id, info in genes.items():
            out_file.write(f"{gene_id},{info['chrom']},{info['intron_count']},{info['total_intron_length']},{info['gene_length']},{info['intron_length_ratio']}\n")

def main():
    args = parse_args()
    genes = parse_gff(args.gff_file)
    write_to_csv(genes, args.output_file)
    print(f"Gene and intron statistics have been saved to {args.output_file}")

if __name__ == "__main__":
    main()

