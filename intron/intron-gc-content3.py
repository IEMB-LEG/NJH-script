import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

# 设置命令行参数
parser = argparse.ArgumentParser(description='Plot GC content distribution from a TAB-separated file.')
parser.add_argument('input_file', type=str, help='Input file path')
parser.add_argument('output_file', type=str, help='Output plot file path')
args = parser.parse_args()

# 读取数据，忽略前7行，使用制表符作为字段分隔符
data = pd.read_csv(args.input_file, delimiter='\t', skiprows=7)

# 选取第6列和第7列的数据
gc_content_col_6 = data.iloc[:, 4]  # 第6列的索引是5
gc_content_col_7 = data.iloc[:, 5]

# 绘制两列数据的分布图
plt.figure(figsize=(12, 6))

# 绘制第6列的分布
sns.kdeplot(gc_content_col_6, color='blue', label='GC content')

# 绘制第7列的分布
sns.kdeplot(gc_content_col_7, color='red', label='GC content trimmed')

plt.xlabel('GC Content')
plt.ylabel('Density')
plt.title('GC Content Distribution')
plt.legend()

# 保存图像到文件
plt.savefig(args.output_file)
plt.close()

