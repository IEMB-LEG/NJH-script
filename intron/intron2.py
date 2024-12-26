import argparse
#generate intron gff3 from non intron annotation gff3. Will rename intron with NO..
def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Generate a GFF3 file with intron annotations based on exon annotations, with unique ID for each intron.")
    parser.add_argument("input_file", type=str, help="Path to the input GFF3 file.")
    parser.add_argument("output_file", type=str, help="Path to the output GFF3 file with intron annotations.")
    return parser.parse_args()

def read_gff3(input_file):
    """
    Reads a GFF3 file and returns a dictionary of mRNA features.
    Extracts exon and mRNA information, skipping non-mRNA features.
    """
    features = {}
    with open(input_file, 'r') as file:
        for line in file:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            attributes = {attr.split('=')[0]: attr.split('=')[1] for attr in parts[8].split(';') if '=' in attr}
            if parts[2] == 'mRNA':
                features[attributes['ID']] = {'exons': [], 'info': parts, 'ID': attributes['ID']}
            elif parts[2] == 'exon' and 'Parent' in attributes:
                parent_id = attributes['Parent']
                if parent_id in features:  # Ensure parent mRNA exists before adding exon
                    features[parent_id]['exons'].append((int(parts[3]), int(parts[4]), parts))
    return features

def write_introns(features, output_file):
    """
    Writes intron annotations to a GFF3 file based on exon information.
    Generates a unique ID for each intron based on its parent mRNA and an incrementing number.
    """
    with open(output_file, 'w') as out:
        out.write("##gff-version 3\n")
        for mRNA, data in features.items():
            data['exons'].sort(key=lambda x: x[0])  # Sort exons by start position
            intron_count = 1  # Initialize intron count for unique ID generation
            for i in range(len(data['exons']) - 1):
                start = data['exons'][i][1] + 1  # Intron start: end of the current exon + 1
                end = data['exons'][i + 1][0] - 1  # Intron end: start of the next exon - 1
                intron_id = f"{data['ID']}-intron-{intron_count}"
                intron_info = '\t'.join([
                    data['info'][0], "custom", "intron", str(start), str(end), ".", data['info'][6], data['info'][7],
                    f"ID={intron_id};Parent={mRNA}"
                ])
                out.write(f"{intron_info}\n")
                intron_count += 1  # Increment intron count for the next ID

def main():
    args = parse_args()
    features = read_gff3(args.input_file)
    write_introns(features, args.output_file)
    print(f"Intron annotations with unique IDs have been successfully written to {args.output_file}")

if __name__ == "__main__":
    main()

