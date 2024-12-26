import sys
import re

def load_correct_ids(filename):
    correct_ids = {}
    with open(filename, 'r') as file:
        for line in file:
            # 假设每行是一个正确的ID
            parts = line.strip().split('g')
            if len(parts) == 2:
                correct_ids[parts[1]] = line.strip()
    return correct_ids

def update_gff(input_gff, output_gff, correct_ids):
    with open(input_gff, 'r') as infile, open(output_gff, 'w') as outfile:
        for line in infile:
            if line.startswith('#'):
                # 直接写入注释行
                outfile.write(line)
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue  # 忽略不完整的行
            
            # 检查序列名称（第一列）是否需要更新或删除
            seq_id_update = False
            match = re.search(r'g(\d+)', parts[0])
            if match and match.group(1) in correct_ids:
                parts[0] = correct_ids[match.group(1)]
                seq_id_update = True
            
            # 检查ID名字（第九列）是否需要更新或删除
            id_update = False
            id_match = re.search(r'ID=([^;]+);', parts[8])
            if id_match:
                id_part = id_match.group(1)
                match = re.search(r'g(\d+)', id_part)
                if match and match.group(1) in correct_ids:
                    parts[8] = re.sub(r'ID=[^;]+;', f'ID={correct_ids[match.group(1)]};', parts[8])
                    id_update = True
            
            # 只有当序列名称和ID名字至少有一个被成功更新时，才写入这一行
            if seq_id_update or id_update:
                outfile.write('\t'.join(parts) + '\n')

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: python script.py <input_gff> <output_gff> <correct_ids>")
        sys.exit(1)
    correct_ids_filename = sys.argv[3]
    correct_ids = load_correct_ids(correct_ids_filename)
    update_gff(sys.argv[1], sys.argv[2], correct_ids)

