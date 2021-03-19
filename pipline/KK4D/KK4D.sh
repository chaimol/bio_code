#!/bin/bash
#此程序是主程序，运行此程序会调用其他脚本。
source config.ini #配置文件
source getKaKs.sh #函数定义脚本
cd ${WorkPath} #进入工作路径


#判断是否是-h或者-v,加载conda信息比较慢
if [ $# -lt 1 ] || [[ $1 = "-h" ]] || [[ $1 = "-V" ]] || [ $1 = "--help" ] || [ $1 = "--version" ];
then
	echo -e "Welcome use KK4D.sh ! \n"
else
	#激活conda环境(shell脚本里，激活conda比较麻烦，需要先source)
	condapath=`conda info | grep 'base environment'|cut -d : -f2|cut -d " " -f2`
	source ${condapath}/etc/profile.d/conda.sh
	conda deactivate
	conda activate mmdetection
	#检测用户的输入文件是否存在
	if [ -e $gff3file1 ] && [ -e $cds1 ] && [ -e $protein1 ] && [ -e $gff3file2 ] && [ -e $cds2 ] && [ -e $protein2 ];then
		#echo "Please check the input file path or make sure the file is exist !"
		echo "Pass the file check !"
	else
		echo "May be you are run not from the first step!"
	fi
fi

function runbed(){
	echo "Begin run analysis of bed in `date "+%Y-%m-%d %H:%M:%S"`"
	#先从gff3获取bed
	if [ $group -eq 2 ];then
		getbed $gff3file1 $prefix1 $type1 $key1
		getbed $gff3file2 $prefix2 $type2 $key2
	else
		getbed $gff3file1 $prefix1 $type1 $key1
	fi
	echo "End run analysis of bed in `date "+%Y-%m-%d %H:%M:%S"`"
}

function runcds(){
	#先判断是否存在bed，否，则运行bed
	if [ ! -e $prefix1.bed ];then
		runbed
	fi
	#开始cds
	echo "Begin run analysis of cds in `date "+%Y-%m-%d %H:%M:%S"`"
	if [ $group -eq 2 ];then
		getcds $cds1 $prefix1
		getcds $cds2 $prefix2
	else 
		getcds $cds1 $prefix1
	fi
	echo "End run analysis of cds in `date "+%Y-%m-%d %H:%M:%S"`"	
}

function runpep(){
	#先判断是否存在bed，否，则运行bed
	if [ ! -e $prefix1.bed ];then
		runbed
	fi
	#再运行pep
	echo "Begin run analysis of pep in `date "+%Y-%m-%d %H:%M:%S"`"
	if [ $group -eq 2 ];then
		getpep $protein1 $prefix1
		getpep $protein2 $prefix2
	else 
		getpep $protein1 $prefix1
	fi
	echo "End run analysis of pep in `date "+%Y-%m-%d %H:%M:%S"`"
}

function runcoline(){
	#先判断是否存在cds，否，则运行cds
	if [ ! -e $prefix1.cds ];then
		runcds
	fi
	#先判断是否存在pep，否，则运行pep
	if [ ! -e $prefix1.pep ];then
		runpep
	fi
	#判断之前运行是否有错误。如果没有错误，则后续出错的概率比较低。
	if [ ! $? -eq 0 ];then
		echo "There is an ERROR before run getcoline!"
		exit 1
	fi
	echo "Begin run analysis of coline in `date "+%Y-%m-%d %H:%M:%S"`"
	if [ $group -eq 2 ];then
		getcoline $prefix1 $prefix2
		#可视化，可能会失败。
		VisualColine $prefix1 $prefix2 $chrnum1 $chrnum2
	else 
		getcoline $prefix1 $prefix1
	fi
	echo "End run analysis of coline in `date "+%Y-%m-%d %H:%M:%S"`"
}

function runprepare(){
	#先判断是否存在共线性文件，否，则运行coline
	if [ ! -e ${prefix1}.${prefix2}.anchors ];then
		runcoline
	fi
	echo "Begin analysis of prepareResult in `date "+%Y-%m-%d %H:%M:%S"`"
	prepareResult $prefix1 $prefix2 $threads
	echo "End analysis of prepareResult in `date "+%Y-%m-%d %H:%M:%S"`"
}

function run4DTV(){
	#判断是否已经生成上一步的输出文件夹result_dir
	if [ ! -d ${prefix1}_${prefix2}.result_dir ];then
		runprepare
	fi
	echo "Begin analysis of 4DTv in `date "+%Y-%m-%d %H:%M:%S"`"
	get4DTv $prefix1 $prefix2
	echo "End analysis of 4DTv in `date "+%Y-%m-%d %H:%M:%S"`"
}

function runKaKs(){
	#判断是否已经生成上一步的输出文件夹result_dir
	if [ ! -d ${prefix1}_${prefix2}.result_dir ];then
		runprepare
	fi
	echo "Begin run analysis of KaKS in `date "+%Y-%m-%d %H:%M:%S"`"
	getkaks $prefix1 $prefix2
	echo "End run analysis of KaKS in `date "+%Y-%m-%d %H:%M:%S"`"
}

function runAll(){

	#判断是否已经生成上一步的输出文件kaks
	if [ ! -e ${prefix1}_${prefix2}.all-kaks.results ];then
		runKaKs
	fi
	#判断是否已经生成上一步的输出文件4dtv
	if [ ! -e ${abbr1}_${abbr2}.all-4dtv.results ];then
		run4DTV
	fi
	echo "Begin run analysis of All begin in `date "+%Y-%m-%d %H:%M:%S"`"
	getkaks4DTv $prefix1 $prefix2
	echo "End run analysis of All in `date "+%Y-%m-%d %H:%M:%S"`"
}

case $1 in
	-h|--help)
	echo "Usage:
	
	run.sh bed/cds/pep/coline/kaks/4DTv/all
	
	Rember:First you must modify the config.ini file.
	"
	;;
	-V|--version)
		  echo -e "
	  Version:${Version}\n
	  Author:${Author} \n
	  Email:${Email} \n
	  Github:${Github} \n
	  Builddate:${Builddate}
	  "
	;;
	bed)
	runbed;;
	cds)
	runcds;;
	pep)
	runpep;;
	coline)
	runcoline;;
	kaks|KaKs)
	runKaKs;;
	4DTv|4DTV)
	run4DTV;;
	all|All)
	runAll;;
	*)
	echo "usage:-h or --help for help!"
esac
