#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
用于替换vcf的染色体字符串
"""
__version__     = "0.1"
__date__        = "2022-06-14"
__Author__      = "Chaimol@163.com"

import sys
import gzip

def replace_chromosomes(vcf_file, mapping_file, output_file):
    chromosome_mapping = {}

    # 读取第二个参数文件，建立染色体编号的映射关系
    with open(mapping_file, 'r') as file:
        for line in file:
            fields = line.strip().split('\t')
            chromosome_mapping[fields[0]] = fields[1]

    try:
        # 打开 VCF 文件（支持压缩文件）
        if vcf_file.endswith('.gz'):
            vcf_lines = gzip.open(vcf_file, 'rt')
        else:
            vcf_lines = open(vcf_file, 'r')

        new_vcf_lines = []
        for line in vcf_lines:
            if line.startswith('#'):
                # 如果是注释行，直接添加到新的 VCF 行列表中
                new_vcf_lines.append(line.rstrip())  # 移除行尾的换行符
            else:
                # 解析每一行的字段
                fields = line.strip().split('\t')
                chromosome = fields[0]

                # 判断染色体编号是否需要替换
                if chromosome in chromosome_mapping:
                    fields[0] = chromosome_mapping[chromosome]

                new_vcf_lines.append('\t'.join(fields))

        # 关闭文件
        vcf_lines.close()

        # 将新的 VCF 行列表写入到输出文件中
        with open(output_file, 'w') as file:
            file.write('\n'.join(new_vcf_lines))

        print("替换完成！")

    except Exception as e:
        print(f"发生错误：{e}")

# 读取命令行参数
args = sys.argv

# 判断用户输入
if len(args) > 3:
    vcf_file = args[1]
    mapping_file = args[2]
    output_file = args[3]
    replace_chromosomes(vcf_file, mapping_file, output_file)

elif len(args) == 2 and (args[1] == "-h" or args[1] == "--help"):
    print("Usage:\n"
          "用于替换VCF的染色体字符串\n"
          "参数1：输入vcf或vcf.gz文件\n"
          "参数2：两列id，使用tab分割，第1列是旧id，第2列是新id\n"
          "参数3：输出文件名\n"
          "python3 replaceVcfChr.py Input.vcf old2newidfile Output.vcf \n")

else:
    print("输入参数有误，请使用-h参数查看示例用法！")
