import pandas as pd
from upsetplot import plot
import matplotlib.pyplot as plt
import sys

def generate_upset_plot(input_file, output_file):
    df = pd.read_csv(input_file, sep='\t', usecols=lambda column: column not in ['Total'])

    # 准备用于UpSet图的数据结构
    data_contents = {}
    for species in df.columns[1:]:  # 跳过Orthogroup列
        # 获取该物种所有非零Orthogroup
        non_zero_orthogroups = df['Orthogroup'][df[species] > 0].tolist()
        if non_zero_orthogroups:
            data_contents[species] = set(non_zero_orthogroups)

    # 生成UpSet数据并绘图
    from upsetplot import from_contents
    upset_data = from_contents(data_contents)
    plot(upset_data, show_counts=True)
    plt.suptitle('UpSet Plot of Orthogroups')
    plt.savefig(output_file)
    plt.close()

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_file> <output_file.pdf>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    generate_upset_plot(input_file, output_file)

