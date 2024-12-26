import argparse
import matplotlib.pyplot as plt
import numpy as np
import os
from collections import defaultdict
from matplotlib.backends.backend_pdf import PdfPages

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Generate aligned line plots for distributions of gene length, exon number, intron length, and exon length from multiple GFF3 files, ensuring data alignment by applying common bins.")
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
    gene_lengths, exon_counts, intron_lengths, exon_lengths = [], [], [], []
    for exons in gene_exons.values():
        sorted_exons = sorted(exons, key=lambda x: x[0])
        gene_lengths.append(sorted_exons[-1][1] - sorted_exons[0][0] + 1)
        exon_counts.append(len(sorted_exons))
        for i in range(1, len(sorted_exons)):
            intron_length = sorted_exons[i][0] - sorted_exons[i-1][1] - 1
            if intron_length > 0:
                intron_lengths.append(intron_length)
        for exon in sorted_exons:
            exon_length = exon[1] - exon[0] + 1
            exon_lengths.append(exon_length)
    return gene_lengths, exon_counts, intron_lengths, exon_lengths

def determine_common_bins(all_data, num_bins=50):
    """Determine common bins based on the 95th percentile of all data."""
    data_95th_percentile = np.percentile(all_data, 95)
    bins = np.linspace(0, data_95th_percentile, num_bins)
    return bins

def plot_distributions(gff_files, output_file):
    """Plot distributions for multiple GFF3 files using common bins."""
    # Aggregate data from all files
    all_gene_lengths, all_exon_counts, all_intron_lengths, all_exon_lengths = [], [], [], []
    for gff_file in gff_files:
        gene_exons = parse_gff(gff_file)
        gene_lengths, exon_counts, intron_lengths, exon_lengths = calculate_statistics(gene_exons)
        all_gene_lengths.extend(gene_lengths)
        all_exon_counts.extend(exon_counts)
        all_intron_lengths.extend(intron_lengths)
        all_exon_lengths.extend(exon_lengths)

    # Determine common bins for each type of data
    gene_length_bins = determine_common_bins(all_gene_lengths)
    exon_count_bins = np.arange(1, max(all_exon_counts) + 2)  # Exon counts are discrete
    intron_length_bins = determine_common_bins(all_intron_lengths)
    exon_length_bins = determine_common_bins(all_exon_lengths)

    with PdfPages(output_file) as pdf:
        fig, axs = plt.subplots(2, 2, figsize=(12, 10))

        for i, gff_file in enumerate(gff_files):
            gene_exons = parse_gff(gff_file)
            gene_lengths, exon_counts, intron_lengths, exon_lengths = calculate_statistics(gene_exons)
            label = os.path.basename(gff_file).split('.')[0]

            # Plot using common bins
            axs[0, 0].hist(gene_lengths, bins=gene_length_bins, histtype='step', label=label)
            axs[0, 1].hist(exon_counts, bins=exon_count_bins, histtype='step', label=label)
            axs[1, 0].hist(intron_lengths, bins=intron_length_bins, histtype='step', label=label)
            axs[1, 1].hist(exon_lengths, bins=exon_length_bins, histtype='step', label=label)

        # Set titles and labels
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
    plot_distributions(args.gff_files, args.output_file)
    print(f"Distributions saved to {args.output_file}")

if __name__ == "__main__":
    main()

