#!/usr/bin/env python3
"""
把vcf文件的第三列的ID的.替换为第1列和第2列拼接的字符
Author: chaimol@163.com
Date: 2023.11.9
"""

import sys
import gzip

def modify_vcf_file(input_file, output_file):
    with gzip.open(input_file, 'rt') if input_file.endswith('.gz') else open(input_file, 'r') as f_in:
        with open(output_file, 'w') as f_out:
            for line in f_in:
                if line.startswith('#'):  # 保留注释行
                    f_out.write(line)
                else:
                    fields = line.split('\t')
                    fields[2] = fields[2].replace('.', fields[0]+'_'+fields[1])
                    f_out.write('\t'.join(fields))

    print(f"修改完成，已保存为 {output_file}")

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("使用方法: python vcfaddID.py input_file output_file")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    modify_vcf_file(input_file, output_file)
