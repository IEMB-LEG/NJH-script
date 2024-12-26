import argparse
from Bio import SeqIO
from Bio.SeqUtils import GC
import numpy as np
import scipy.stats

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Calculate and summarize statistics for each contig in a FASTA file, correctly checking for TGA codons, with summary at the beginning.")
    parser.add_argument("input_file", type=str, help="Path to the input FASTA file.")
    parser.add_argument("output_file", type=str, help="Path to the output file with statistics.")
    return parser.parse_args()

def has_tga_codon(sequence):
    """Check for TGA codon presence in the correct reading frame."""
    for i in range(0, len(sequence)-2, 3):
        if sequence[i:i+3] == "TGA":
            return True
    return False

def calc_stats_and_write(input_file, output_file):
    """Calculate statistics and write to the output file."""
    stats = {
        "lengths": [],
        "div_by_three": 0,
        "not_div_by_three": 0,
        "div_and_stop": 0,
        "not_div_and_stop": 0,
        "gc_contents": [],
        "gc_contents_trimmed": []
    }

    contig_details = []

    for record in SeqIO.parse(input_file, "fasta"):
        length = len(record.seq)
        divisible_by_three = length % 3 == 0
        contains_tga = has_tga_codon(record.seq)
        gc_content = GC(record.seq)
        gc_content_trimmed = GC(record.seq[2:-2]) if length > 4 else GC(record.seq)

        stats["lengths"].append(length)
        stats["gc_contents"].append(gc_content)
        stats["gc_contents_trimmed"].append(gc_content_trimmed)

        if divisible_by_three:
            stats["div_by_three"] += 1
            if contains_tga:
                stats["div_and_stop"] += 1
        else:
            stats["not_div_by_three"] += 1
            if contains_tga:
                stats["not_div_and_stop"] += 1

        contig_details.append((record.id, length, divisible_by_three, contains_tga, gc_content, gc_content_trimmed))

    # Calculate summary statistics
    mean_length, ci_length = np.mean(stats["lengths"]), scipy.stats.t.interval(0.95, len(stats["lengths"])-1, loc=np.mean(stats["lengths"]), scale=scipy.stats.sem(stats["lengths"]))
    mean_gc, ci_gc = np.mean(stats["gc_contents"]), scipy.stats.t.interval(0.95, len(stats["gc_contents"])-1, loc=np.mean(stats["gc_contents"]), scale=scipy.stats.sem(stats["gc_contents"]))
    mean_gc_trimmed, ci_gc_trimmed = np.mean(stats["gc_contents_trimmed"]), scipy.stats.t.interval(0.95, len(stats["gc_contents_trimmed"])-1, loc=np.mean(stats["gc_contents_trimmed"]), scale=scipy.stats.sem(stats["gc_contents_trimmed"]))

    # Write results to file
    with open(output_file, 'w') as out:
        out.write(f"# Summary Statistics:\n")
        out.write(f"Average Length: {mean_length:.2f}, 95% CI: [{ci_length[0]:.2f}, {ci_length[1]:.2f}]\n")
        out.write(f"Divisible by Three: {stats['div_by_three']}, Not Divisible by Three: {stats['not_div_by_three']}\n")
        out.write(f"Divisible and Contains TGA: {stats['div_and_stop']}, Not Divisible and Contains TGA: {stats['not_div_and_stop']}\n")
        out.write(f"Average GC Content: {mean_gc:.2f}, 95% CI: [{ci_gc[0]:.2f}, {ci_gc[1]:.2f}]\n")
        out.write(f"Average GC Content (Trimmed): {mean_gc_trimmed:.2f}, 95% CI: [{ci_gc_trimmed[0]:.2f}, {ci_gc_trimmed[1]:.2f}]\n\n")
        out.write("ContigID\tLength\tDivisibleByThree\tContainsTGA\tGCContent\tGCContentTrimmed\n")
        for detail in contig_details:
            out.write(f"{detail[0]}\t{detail[1]}\t{'Yes' if detail[2] else 'No'}\t{'Yes' if detail[3] else 'No'}\t{detail[4]:.2f}\t{detail[5]:.2f}\n")

def main():
    args = parse_args()
    calc_stats_and_write(args.input_file, args.output_file)
    print(f"Statistics and contig details, including summary, have been saved to {args.output_file}")

if __name__ == "__main__":
    main()

