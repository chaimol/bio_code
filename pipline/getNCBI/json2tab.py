#!/usr/bin/env python3
##用法：python3 json2tab.py NCBI.json
#输入文件是ffq输出的json文件
#输出是"projectID","sampleID","fq_url","size","md5","organism"
import sys
import pandas as pd
import json
#jsonfile = "NCBI.json"
jsonfile = sys.argv[1]
#filename = jsonfile.split('.')[0]
#解析输入的json文件
with open(jsonfile,"r") as load_f:
    load_dict = json.load(load_f)
#输出文件头部信息
#print("projectID","sampleID","fq_url","size","md5","organism")
for i in load_dict.keys(): #i是项目的编号，当有多个NCBI的号的时候
    idlist=[]
    for j in load_dict[i].keys():
        idlist.append(j) #把每一个的keys输出到idlist,如果keys里有files选项，则直接输出，如果有runs，需要多解析一层
    if "files" in idlist:
        data_ID=load_dict[i] #此时有files直接解析即可
        for list_url in data_ID['files']:
            print(load_dict[i]['accession'],data_ID['accession'],list_url['url'],list_url['size'],list_url['md5'],data_ID['sample']['organism'],sep='\t')
    elif "runs" in idlist: #此时是需要多解析一层
        for id in load_dict[i]['runs'].keys(): #id是每个项目里，对应的文件的号
            data_ID=load_dict[i]['runs'][id]
            for list_url in data_ID['files']:
                print(load_dict[i]['accession'],data_ID['accession'],list_url['url'],list_url['size'],list_url['md5'],data_ID['sample']['organism'],sep='\t')
    else: #此时可能是其他情况，自行解决
        print("不符合已知的字段，需要自行解析如下内容：")
        print(load_dict[i])