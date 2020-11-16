#!/bin/bash

#SBATCH -J QC #作业名
#SBATCH -N 1 #节点数量
#SBATCH -n 4 # 运行的总核数量
#SBATCH -t 192:00:00  #最大运行时间
##SBATCH --partition=compute*  #设置运行的分区，不同的分区的硬件不同（中心的服务器有2个分区，compute*是小节点nodes=02-09，big是大节点,nodes=01,）
# set batch script's standard output
#SBATCH --output=QC.out

echo " my job id is $START:$SLURM_JOBID "| tee QC.log
echo run nodes is following: | tee -a QC.log

echo begin time is `date` | tee -a QC.log
id=`echo $PBS_JOBID|awk -F. '{print $1}' `
NP=`cat $PBS_NODEFILE|wc -l`

cd /share/home/chaimao/wes_cancer/project/
# run the application
srun sh QC.sh
echo end time is `date` | tee -a QC.log

##slurm资源管理系统命令
#运行命令方式：
#sbatch sh QC.sh