#!/bin/bash
##使用SimpleM自动分析计算最佳P的阈值的脚本
#判断是否是-h或者-v,加载conda信息比较慢
if [ $# -lt 1 ] || [[ $1 = "-h" ]] || [[ $1 = "-V" ]] || [ $1 = "--help" ] || [ $1 = "--version" ];
then
	echo -e "Usage:bash getPvalue.bash vcfile abbr \n
			输入两个参数：1是vcf文件，2是输出前缀 \n
			输出的结果在 abbr.SimpleM.pvalue \n
			结果文件的第三行即为P的阈值 \n
			"
	exit
fi

input_vcf=$1
abbr=$2


#用于获取脚本所在的路径，保存为变量path1,调用其他脚本都依赖这个路径。
path1="$(cd "$(dirname ${BASH_SOURCE[0]})";pwd)"

##软件路径配置
beagle="java -jar ${path1}/beagle.22Jul22.46e.jar"
plink="${path1}/plink"
SimpleM=".Rscript ${path1}/SimpleM.R"
kggsee="java -Xmx10g -jar ${path1}/kggsee.jar"

#使用beagle对vcf文件进行填充
$beagle gt=${input_vcf} out=${abbr}.beagle
#输出文件是${abbr}.beagle.vcf.gz

#使用plink过滤输入的数据
$plink --vcf ${abbr}.beagle.vcf.gz --maf 0.05 --geno 0.1 --recode vcf-iid --out ${abbr}.filter
#输出文件是${abbr}.filter.vcf

#把vcf文件转为numeric格式的矩阵只有0,1,2
cat ${abbr}.filter.vcf | grep -v "^##" | cut -f 10- | sed 's/0\/0/0/g' | sed 's/1\/1/2/g' | sed 's/0\/1/1/g' | sed 's/1\/0/1/g'| sed 's/\.\/\./-1/g' | tr "\t" " "|sed '1d' >${abbr}.genotype_transpose.txt

##使用SimpleM.R计算P的阈值
${SimpleM} ${abbr}.genotype_transpose.txt 0.05 >${abbr}.SimpleM.pvalue