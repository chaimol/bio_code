#!/usr/bin/env python3
"""
tsv和pkl文件相互转换工具
Author: chaimol@163.com
Date: 2023.11.8
"""

import pandas as pd
import sys

# 检查命令行参数数量
if len(sys.argv) != 3:
    print("Usage: python script.py input_file output_file")
    sys.exit(1)

# 获取命令行参数
input_file = sys.argv[1]
output_file = sys.argv[2]

# 判断操作类型并进行相应操作
if input_file.endswith('.tsv') and output_file.endswith('.pkl'):
    # 将TSV文件转换为PKL文件
    df = pd.read_csv(input_file, sep='\t',index_col=0)
    df.to_pickle(output_file)
    print("TSV to PKL conversion completed!")

elif input_file.endswith('.pkl') and output_file.endswith('.tsv'):
    # 将PKL文件转换为TSV文件
    df = pd.read_pickle(input_file)
    df.to_csv(output_file, sep='\t', index=True)
    print("PKL to TSV conversion completed!")

else:
    print("Invalid file format!")
    sys.exit(1)
