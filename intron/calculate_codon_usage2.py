import sys
from Bio import SeqIO
from Bio.SeqUtils import CodonUsage
import matplotlib.pyplot as plt
import pandas as pd
from collections import Counter

def parse_gff3_cds(gff3_file):
    """
    从GFF3文件解析CDS的位置信息。
    返回一个包含(序列ID, 起始位置, 终止位置, 链)的元组列表。
    """
    cds_locations = []
    with open(gff3_file, 'r') as gff:
        for line in gff:
            if line.startswith('#') or line.strip() == '':
                continue
            parts = line.split('\t')
            if len(parts) < 9:  # 确保行被正确分割成足够多的部分
                continue
            if parts[2] == 'CDS':
                seqid = parts[0]
                start = int(parts[3])
                end = int(parts[4])
                strand = parts[6]
                cds_locations.append((seqid, start, end, strand))
    return cds_locations

def extract_cds_sequences(fasta_file, cds_locations):
    """
    根据CDS位置信息，从FASTA文件中提取序列。
    返回CDS序列的列表。
    """
    cds_sequences = []
    fasta_sequences = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    for seqid, start, end, strand in cds_locations:
        if seqid in fasta_sequences:
            seq = fasta_sequences[seqid].seq
            if strand == '+':
                cds_seq = seq[start-1:end]
            else:
                cds_seq = seq[start-1:end].reverse_complement()
            cds_sequences.append(cds_seq)
    return cds_sequences

def count_codon_usage(cds_sequences):
    """
    统计CDS序列中的密码子使用频率。
    返回密码子使用频率的字典。
    """
    counter = Counter()
    for seq in cds_sequences:
        for i in range(0, len(seq) - 2, 3):
            codon = str(seq[i:i+3])
            if codon != 'TGA':  # 跳过终止密码子
                counter[codon] += 1
    return counter

import numpy as np

def plot_codon_usage(codon_usage, output_pdf):
    """
    根据密码子使用频率和GC含量，绘制并保存图表。
    使用四种固定颜色来表示不同的GC含量区间，并添加自定义图例。
    """
    df = pd.DataFrame(codon_usage.items(), columns=['Codon', 'Frequency'])
    df['GC_Content'] = df['Codon'].apply(lambda x: (x.count('G') + x.count('C')) / 3)
    df_sorted = df.sort_values(by='Frequency', ascending=False)

    # 定义四个GC含量区间对应的颜色
    gc_colors = {
        'Low': 'blue',
        'Moderate': 'green',
        'High': 'orange',
        'Very High': 'red'
    }
    # 分配颜色
    def assign_color(gc_content):
        if gc_content <= 0.33:
            return gc_colors['Low']
        elif gc_content <= 0.66:
            return gc_colors['Moderate']
        elif gc_content < 1.0:
            return gc_colors['High']
        else:
            return gc_colors['Very High']
    
    df_sorted['Color'] = df_sorted['GC_Content'].apply(assign_color)
    
    plt.figure(figsize=(10, 8))
    for color in gc_colors.values():
        plt.bar(df_sorted['Codon'][df_sorted['Color'] == color], df_sorted['Frequency'][df_sorted['Color'] == color], color=color, label=color)
    
    plt.xlabel('Codon')
    plt.ylabel('Frequency')
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.legend(title='GC Content')
    
    plt.savefig(output_pdf, format='pdf')


if __name__ == '__main__':
    fasta_file = sys.argv[1]
    gff3_file = sys.argv[2]
    output_pdf = sys.argv[3]

    cds_locations = parse_gff3_cds(gff3_file)
    cds_sequences = extract_cds_sequences(fasta_file, cds_locations)
    codon_usage = count_codon_usage(cds_sequences)
    plot_codon_usage(codon_usage, output_pdf)

