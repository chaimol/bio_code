#!/bin/bash
#此程序是主程序，运行此程序会调用其他脚本。
source getKaKs.sh
#激活conda环境(shell脚本里，激活conda比较麻烦，需要先source)
condapath=`conda info | grep 'base environment'|cut -d : -f2|cut -d " " -f2`
source ${condapath}/etc/profile.d/conda.sh
conda deactivate
conda activate mmdetection

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
	if [! -e $prefix1.bed ];then
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
	if [! -e $prefix1.bed ];then
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
	if [! -e $prefix1.cds ];then
		runcds
	fi
	#先判断是否存在pep，否，则运行pep
	if [! -e $prefix1.pep ];then
		runpep
	fi
	echo "Begin run analysis of coline in `date "+%Y-%m-%d %H:%M:%S"`"
	if [ $group -eq 2 ];then
		getcoline $prefix1 $prefix2
		#可视化，可能会失败。
		VisualColine $abbr1 $abbr2 $chrnum1 $chrnum2
	else 
		getcoline $prefix1 $prefix1
	fi
	echo "End run analysis of coline in `date "+%Y-%m-%d %H:%M:%S"`"
}

function run4DTV(){
	#判断是否已经生成上一步的输出文件夹result_dir
	if [ ! -d ${prefix1}_${prefix2}.result_dir ];then
		echo "Begin analysis of prepareResult in `date "+%Y-%m-%d %H:%M:%S"`"
		prepareResult $prefix1 $prefix2 $threads
		echo "End analysis of prepareResult in `date "+%Y-%m-%d %H:%M:%S"`"
	fi
	echo "Begin analysis of 4DTv in `date "+%Y-%m-%d %H:%M:%S"`"
	get4DTv $prefix1 $prefix2
	echo "End analysis of 4DTv in `date "+%Y-%m-%d %H:%M:%S"`"
}

function runKaKs(){
	#判断是否已经生成上一步的输出文件夹result_dir
	if [ ! -d ${prefix1}_${prefix2}.result_dir ];then
		echo "Begin analysis of prepareResult in `date "+%Y-%m-%d %H:%M:%S"`"
		prepareResult $prefix1 $prefix2 $threads
		echo "End analysis of prepareResult in `date "+%Y-%m-%d %H:%M:%S"`"
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
	echo "Usage:run.sh bed/cds/pep/coline/kaks/4DTv/all
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
