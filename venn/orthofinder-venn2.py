import pandas as pd
import sys

def convert_orthofinder_output(input_file, output_file):
    # 读取OrthoFinder输出文件，忽略最后一列'Total'
    df = pd.read_csv(input_file, sep='\t', usecols=lambda column: column not in ['Total'])

    # 创建一个空列表来收集所有的行
    rows = []

    # 遍历DataFrame的每一行
    for index, row in df.iterrows():
        orthogroup = row['Orthogroup']
        # 对于每个物种，如果该物种在当前Orthogroup中的基因数不为0，则添加相应数量的行
        for species in df.columns[1:]:  # 跳过第一列'Orthogroup'
            # 检查是否为空值（NaN），是则替换为0
            gene_count = row[species] if pd.notnull(row[species]) else 0
            gene_count = int(gene_count)  # 确保将基因计数转换为整数
            for _ in range(gene_count):
                rows.append({'OTU': orthogroup, 'Sample': species})

    # 将收集到的行转换为DataFrame
    formatted_df = pd.DataFrame(rows, columns=['OTU', 'Sample'])

    # 保存转换后的数据到新文件
    formatted_df.to_csv(output_file, sep='\t', index=False)

    print(f'Converted data saved to {output_file}')

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_file> <output_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    convert_orthofinder_output(input_file, output_file)

