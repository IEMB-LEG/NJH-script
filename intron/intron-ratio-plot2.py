import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Plot the relationship between intron count and intron length ratio for each gene with density information and save the plot to a PDF file.")
    parser.add_argument("input_file", type=str, help="Path to the input CSV file containing gene and intron statistics.")
    parser.add_argument("output_file", type=str, help="Path to the output PDF file for the plot.")
    return parser.parse_args()

def plot_intron_statistics_with_density(input_file, output_file):
    """Read gene statistics, plot intron count vs. intron length ratio with density, and save to PDF."""
    df = pd.read_csv(input_file)

    plt.figure(figsize=(10, 3))
    sns.scatterplot(data=df, x='IntronCount', y='IntronLengthRatio', hue='IntronCount', size='IntronLengthRatio', sizes=(20, 200), alpha=0.5, palette="viridis")
    plt.title('Relationship between Intron Count and Intron Length Ratio with Density')
    plt.xlabel('Intron Count')
    plt.ylabel('Intron Length Ratio')
    plt.grid(True)
    plt.legend(title='Intron Count & Length Ratio', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    
    # Save the plot to a PDF file
    plt.savefig(output_file, format='pdf')

def main():
    args = parse_args()
    plot_intron_statistics_with_density(args.input_file, args.output_file)
    print(f"Density plot successfully saved to {args.output_file}")

if __name__ == "__main__":
    main()

