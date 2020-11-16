#!/bin/bash
##文件重命名脚本
##在下载的clean文件夹运行该脚本
##自动备份数据到origin，同时重命名为RNA-seq要求的文件格式

##### awk -F , 'NR!=1{print $2 >$1".txt"}' geneinfo.csv
#####逐行分列输出第二列内容到第一列文件名的文件。
# mkdir origin
# cp *.gz ./origin


##清除相关文件
rm -rf name R_info sample_name name0 name3 new_name

###获取所有测序数据名称
ls *.gz >name

###存储材料名
cat name |awk -F "_" '{print $1}' >sample_name

###存储测序的端
cat name |awk -F "." '{print $2}' >R_info
###获取所有文件的数量
linecount=`sed -n '$=' name`
for ((i=0;i<$linecount;i++))
do
	echo "mv" >>name0
	echo "." >>name2
	echo ".fq.gz" >>name3
done
###组装材料名的信息（材料名+R1/R2+.fq.gz）
paste -d "" sample_name name2 R_info name3 >new_name 
paste -d \  name0 name new_name  >1.bash

rm -rf name R_info sample_name name0 name2 name3 new_name

#bash 1.bash

echo "complete exchange file name!"
