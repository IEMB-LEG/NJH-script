#delete >95%
import argparse
import matplotlib.pyplot as plt
import numpy as np
import os
from collections import defaultdict
from matplotlib.backends.backend_pdf import PdfPages

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Generate line plots for distributions of gene length, exon number, intron length, and exon length from multiple GFF3 files, focusing on the 95% data.")
    parser.add_argument("gff_files", nargs='+', help="Paths to the GFF3 files.")
    parser.add_argument("output_file", type=str, help="Path to the output PDF file.")
    return parser.parse_args()

def parse_gff(gff_file):
    """Parse GFF3 file to extract required information for plotting."""
    gene_exons = defaultdict(list)
    
    with open(gff_file, 'r') as file:
        for line in file:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            feature_type, start, end, attributes = parts[2], int(parts[3]), int(parts[4]), parts[8]
            if feature_type == 'exon':
                gene_id = [attr for attr in attributes.split(';') if attr.startswith('Parent=')][0].split('Parent=')[1]
                gene_exons[gene_id].append((start, end))
                
    return gene_exons

def calculate_statistics(gene_exons):
    """Calculate statistics for gene length, exon number, intron length, and exon length."""
    gene_lengths = []
    exon_counts = []
    intron_lengths = []
    exon_lengths = []
    
    for exons in gene_exons.values():
        sorted_exons = sorted(exons, key=lambda x: x[0])
        gene_lengths.append(sorted_exons[-1][1] - sorted_exons[0][0] + 1)
        
        exon_count = len(sorted_exons)
        exon_counts.append(exon_count)
        
        for i in range(1, len(sorted_exons)):
            intron_length = sorted_exons[i][0] - sorted_exons[i-1][1] - 1
            if intron_length > 0:
                intron_lengths.append(intron_length)
                
        for exon in sorted_exons:
            exon_length = exon[1] - exon[0] + 1
            exon_lengths.append(exon_length)
            
    return gene_lengths, exon_counts, intron_lengths, exon_lengths

def plot_statistics(gff_files, output_file):
    """Plot statistics for multiple GFF3 files, focusing on the 95% data."""
    with PdfPages(output_file) as pdf:
        fig, axs = plt.subplots(2, 2, figsize=(12, 10))
        colors = plt.cm.tab10.colors

        for i, gff_file in enumerate(gff_files):
            gene_exons = parse_gff(gff_file)
            gene_lengths, exon_counts, intron_lengths, exon_lengths = calculate_statistics(gene_exons)
            label = os.path.basename(gff_file).split('.')[0]

            # Calculate 95th percentile and filter data
            gene_length_95th = np.percentile(gene_lengths, 95)
            exon_count_95th = np.percentile(exon_counts, 95)
            intron_length_95th = np.percentile(intron_lengths, 95)
            exon_length_95th = np.percentile(exon_lengths, 95)

            # Filter data
            filtered_gene_lengths = [l for l in gene_lengths if l <= gene_length_95th]
            filtered_exon_counts = [c for c in exon_counts if c <= exon_count_95th]
            filtered_intron_lengths = [l for l in intron_lengths if l <= intron_length_95th]
            filtered_exon_lengths = [l for l in exon_lengths if l <= exon_length_95th]

            # Plot
            axs[0, 0].hist(filtered_gene_lengths, bins=50, histtype='step', label=label, color=colors[i % len(colors)])
            axs[0, 1].hist(filtered_exon_counts, bins=50, histtype='step', label=label, color=colors[i % len(colors)])
            axs[1, 0].hist(filtered_intron_lengths, bins=50, histtype='step', label=label, color=colors[i % len(colors)])
            axs[1, 1].hist(filtered_exon_lengths, bins=50, histtype='step', label=label, color=colors[i % len(colors)])

        titles = ['Gene Length Distribution', 'Exon Number Distribution', 'Intron Length Distribution', 'Exon Length Distribution']
        for ax, title in zip(axs.flat, titles):
            ax.set_title(title)
            ax.set_xlabel('Length or Count')
            ax.set_ylabel('Frequency')
            ax.legend()

        plt.tight_layout()
        pdf.savefig(fig)
        plt.close()

def main():
    args = parse_args()
    plot_statistics(args.gff_files, args.output_file)
    print(f"Distributions saved to {args.output_file}")

if __name__ == "__main__":
    main()

