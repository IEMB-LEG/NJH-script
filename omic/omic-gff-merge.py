import sys
import re

def load_cds_descriptions(cds_gff_file):
    descriptions = {}
    with open(cds_gff_file, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue  # 忽略不完整或无效的行
            id_match = re.search(r'ID=([^;]+)', parts[8])
            if id_match:
                id_value = id_match.group(1)
                # 为mRNA和对应的gene保存相同的描述
                if id_value.startswith('mRNA:'):
                    gene_id = 'gene:' + id_value[5:]
                    descriptions[gene_id] = ';'.join(parts[8].split(';')[1:])  # 获取除ID外的描述
                descriptions[id_value] = ';'.join(parts[8].split(';')[1:])  # 获取除ID外的描述
    return descriptions

def merge_descriptions(genome_gff3_file, cds_descriptions, output_file):
    with open(genome_gff3_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('#'):
                outfile.write(line)
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                outfile.write(line)  # 直接写入不完整或无效的行
                continue
            id_match = re.search(r'ID=([^;]+)', parts[8])
            if id_match and id_match.group(1) in cds_descriptions:
                # 将CDS描述合并到全基因组GFF3记录中
                parts[8] += ';' + cds_descriptions[id_match.group(1)]
            outfile.write('\t'.join(parts) + '\n')

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: python script.py <genome_gff3> <cds_gff> <output_gff3>")
        sys.exit(1)
    cds_gff_file = sys.argv[2]
    cds_descriptions = load_cds_descriptions(cds_gff_file)
    merge_descriptions(sys.argv[1], cds_descriptions, sys.argv[3])

