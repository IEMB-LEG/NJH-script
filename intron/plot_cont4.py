import pandas as pd
import matplotlib.pyplot as plt
import argparse
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Generate a stacked bar chart focusing on 'No' counts with custom colors, save as PDF.")
    parser.add_argument("input_file", type=str, help="Path to the input file.")
    parser.add_argument("output_file", type=str, help="Path to the output PDF file.")
    parser.add_argument("--interval", type=int, default=10, help="Interval for length categories.")
    parser.add_argument("--max_length", type=int, default=70, help="Maximum length for categorization before lumping into a final category.")
    parser.add_argument("--y_max", type=int, default=100, help="Maximum value for the y-axis.")  # Add y-axis maximum as a command line argument
    return parser.parse_args()

def read_data(input_file):
    """Read the input file into a DataFrame."""
    df = pd.read_csv(input_file, sep="\t", header=None, names=["Length", "FirstCol", "SecondCol"])
    return df

def categorize_data(df, interval, max_length):
    """Categorize length data and calculate counts focusing on 'No'."""
    bins = np.arange(1, max_length, interval).tolist() + [max_length, df['Length'].max() + 1]
    labels = [f"{i}-{i+interval-1}" for i in range(1, max_length, interval)] + [f"{max_length}+"]
    df['Category'] = pd.cut(df['Length'], bins=bins, labels=labels, right=False)

    total_no = df[df['FirstCol'] == "No"].groupby('Category').size()
    double_no = df[(df['FirstCol'] == "No") & (df['SecondCol'] == "No")].groupby('Category').size()

    summary = pd.DataFrame({'TotalNo': total_no, 'DoubleNo': double_no}).fillna(0)
    return summary

def plot_data(summary, output_file, y_max):
    """Plot the data with custom colors, save as PDF, set Y-axis limit, and display the plot."""
    with PdfPages(output_file) as pdf:
        colors = ['#fbb4ae', '#b3cde3']  # Custom colors: a darker and a lighter shade
        ax = summary.plot(kind='bar', stacked=True, figsize=(10, 6), color=colors)
        ax.set_ylim([0, y_max])  # Set Y-axis upper limit
        plt.title('Stacked Bar Chart of No Counts')
        plt.xlabel('Length Category')
        plt.ylabel('Count')
        plt.xticks(rotation=45)
        plt.legend(['Total No', 'Both Columns No'], loc='upper right')
        plt.tight_layout()
        pdf.savefig()
        plt.show()

def main():
    args = parse_args()
    df = read_data(args.input_file)
    summary = categorize_data(df, args.interval, args.max_length)
    plot_data(summary, args.output_file, args.y_max)  # Pass the Y-axis max to the plotting function
    print(f"Stacked bar chart with custom colors saved to {args.output_file}")

if __name__ == "__main__":
    main()

