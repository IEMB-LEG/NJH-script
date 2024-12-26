import argparse
import matplotlib.pyplot as plt
import numpy as np
import os
from collections import defaultdict
from matplotlib.backends.backend_pdf import PdfPages
#from scipy.interpolate import interp1d
from scipy.interpolate import make_interp_spline 
def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Generate line plots for distributions with specified max values for the x-axis.")
    parser.add_argument("gff_files", nargs='+', help="Paths to the GFF3 files.")
    parser.add_argument("output_file", type=str, help="Path to the output PDF file.")
    parser.add_argument("--max_gene_length", type=int, required=True, help="Maximum value for gene length.")
    parser.add_argument("--max_exon_count", type=int, required=True, help="Maximum value for exon count.")
    parser.add_argument("--max_intron_length", type=int, required=True, help="Maximum value for intron length.")
    parser.add_argument("--max_exon_length", type=int, required=True, help="Maximum value for exon length.")
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
        exon_counts.append(len(sorted_exons))
        for i in range(1, len(sorted_exons)):
            intron_lengths.append(sorted_exons[i][0] - sorted_exons[i-1][1] - 1)
        for exon in sorted_exons:
            exon_lengths.append(exon[1] - exon[0] + 1)        
    return gene_lengths, exon_counts, intron_lengths, exon_lengths

def plot_statistics(gff_files, output_file, max_values):
    """Plot statistics for multiple GFF3 files, applying max values for x-axis, smoothing lines for specific distributions, and annotating bin lengths."""
    max_gene_length, max_exon_count, max_intron_length, max_exon_length = max_values
    with PdfPages(output_file) as pdf:
        fig, axs = plt.subplots(2, 2, figsize=(12, 10))
        colors = plt.cm.tab10.colors

        for i, gff_file in enumerate(gff_files):
            gene_exons = parse_gff(gff_file)
            gene_lengths, exon_counts, intron_lengths, exon_lengths = calculate_statistics(gene_exons)
            label = os.path.basename(gff_file).split('.')[0]

            datasets = [gene_lengths, exon_counts, intron_lengths, exon_lengths]
            maxes = [max_gene_length, max_exon_count, max_intron_length, max_exon_length]
            titles = ['Gene Length Distribution', 'Exon Number Distribution', 'Intron Length Distribution', 'Exon Length Distribution']

            for dataset, ax, title, max_val in zip(datasets, axs.flat, titles, maxes):
                filtered_data = np.array([d for d in dataset if d <= max_val])

                # Compute histogram
                bins = 50
                counts, bin_edges = np.histogram(filtered_data, bins=bins, density=False)
                bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
                bin_width = bin_edges[1] - bin_edges[0]

                # For 'Gene Length Distribution' and 'Exon Length Distribution', smooth the line using bins
                if title in ['Gene Length Distribution', 'Exon Length Distribution']:
                    if len(bin_centers) > 1:  # Need at least two points to interpolate
                        spline = make_interp_spline(bin_centers, counts, k=3)  # BSpline object
                        smooth_bin_centers = np.linspace(bin_centers.min(), bin_centers.max(), 300)
                        smooth_counts = spline(smooth_bin_centers)
                        ax.plot(smooth_bin_centers, smooth_counts, label=label, color=colors[i % len(colors)])
                        # Annotate bin width
                        ax.annotate(f'Bin width: {bin_width}', xy=(0.7, 0.9), xycoords='axes fraction', fontsize=9, color='blue')
                    else:
                        ax.scatter(bin_centers, counts, label=label, color=colors[i % len(colors)], zorder=5)
                else:
                    # Smooth line for other distributions
                    values, counts = np.unique(filtered_data, return_counts=True)
                    if len(values) > 1:
                        spline = make_interp_spline(values, counts, k=3)
                        smooth_values = np.linspace(values.min(), values.max(), 300)
                        smooth_counts = spline(smooth_values)
                        ax.plot(smooth_values, smooth_counts, label=label, color=colors[i % len(colors)])
                    else:
                        ax.scatter(values, counts, label=label, color=colors[i % len(colors)], zorder=5)

                ax.set_title(title)
                ax.set_xlabel('Length or Count')
                ax.set_ylabel('Frequency')
                ax.legend()

        plt.tight_layout()
        pdf.savefig(fig)
        plt.close()

def main():
    args = parse_args()
    max_values = (args.max_gene_length, args.max_exon_count, args.max_intron_length, args.max_exon_length)
    plot_statistics(args.gff_files, args.output_file, max_values)
    print(f"Distributions saved to {args.output_file}")

if __name__ == "__main__":
    main()
