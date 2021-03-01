#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
用于修改fasta文件的头部，使用前先根据实际修改split_str函数的切割命令
"""

__author__      = "Frank Chai"
__credits__     = "Frank Chai"
__version__     = "0.1"
__email__       = "chaimol@163.com"
__date__        = "2021-02-26"
import sys
import os
import gzip
import Bio


#读取用户的输入
args = sys.argv

#定义对description切片的语法（每次运行前，先根据实际修改切片语法）
def split_head(description):
    #指定description的第4个，用=分割1次，取第1个
    split_str=description.split()[3].split("=",1)[1]
    return split_str

#读取普通fasta文件，修改文件头部
def re_header():
    '''
    #读取fasta,并更改header
    '''
    from Bio import SeqIO
    #读取fasta
    for seq_record in SeqIO.parse(InputFa,"fasta"):
        #print(seq_record.id)原始的fasta的id
        #print(seq_record.seq)
        #print(len(seq_record))
        #print(seq_record.description) 原始的fasta的id后的内容
        #指定description的第4个，用=分割1次，取第1个
        #此处的分割需要根据实际进行修改
        description_id=split_head(seq_record.description)
        liststr= ['>',description_id]
        new_id= ''.join(liststr)
        print(new_id)
        print(seq_record.seq)

#读取fasta.gz的文件,修改文件头部
def re_gz_header():
    '''
    read file.fasta.gz  
    '''
    from Bio import SeqIO
    with gzip.open(InputFa, "rt") as handle:
        for seq_record in SeqIO.parse(handle, "fasta"):
            #print(seq_record.id)
            #print(seq_record.seq)
            #print(len(seq_record))
            #指定description的第4个，用=分割1次，取第1个
            #此处的分割需要根据实际进行修改
            description_id=split_head(seq_record.description)
            liststr= ['>',description_id]
            new_id= ''.join(liststr)
            print(new_id)
            print(seq_record.seq)


#判断用户输入的是fasta还是fasta.gz
def judge_Input():
    if InputFa.endswith(".gz"):
        re_gz_header()
    else:
        re_header()

#判断用户输入
if len(args) >1:
    if args[1] == "-h" or args[1] == "--help":
        print("Usage:(输入文件支持fasta或fasta.gz)\n"
        "用于修改fasta文件的头部\n"
        "python3 farename.py test.fa #使用切片输出新的fasta.\n"
        "python3 farename.py test.fa.gz #使用切片输出新的fasta.\n"
        "Notes:切片前先检测你的fasta文件的头部的格式，之后修改函数split_head里的切片参数")
    elif len(args) == 2:
        InputFa = args[1]
        judge_Input()
    else:
        print("输入参数有误，请使用-h参数查看示例用法！")
else:
    print("输入参数有误，请使用-h参数查看示例用法！")
