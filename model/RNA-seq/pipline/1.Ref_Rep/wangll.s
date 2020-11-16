##PBS任务系统提交命令程序
#PBS -N wangll
#PBS -l nodes=1:ppn=24
#PBS -j oe
#PBS -l walltime=756:00:00
#PBS -V
#PBS -o /public/home/tang/wangll/wangll.out
#PBS -e /public/home/tang/wangll/wangll.err

echo my job id is $PBS_JOBID | tee wangll.log
echo run nodes is following: | tee -a wangll.log
cat $PBS_NODEFILE | tee  -a wangll.log

echo begin time is `date` | tee -a wangll.log
id=`echo $PBS_JOBID|awk -F. '{print $1}' `
NP=`cat $PBS_NODEFILE|wc -l`

cd /public/home/tang/wangll/
bash wangll.bash
echo end time is `date` | tee -a wangll.log
