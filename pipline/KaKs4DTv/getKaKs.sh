#!/bin/bash
####此程序是函数定义脚本。主运行脚本是run.sh

##读取配置文件
#source config.ini
if [ $group -eq 2 ];then
	prefix1=${abbr[0]}
	prefix2=${abbr[1]}
	gff3file1=${gff3[0]}
	gff3file2=${gff3[1]}
	latin1=${sample[0]}
	latin2=${sample[1]}
	protein1=${protein[0]}
	protein2=${protein[1]}
	cds1=${cds[0]}
	cds2=${cds[1]}
	key1=${key[0]}
	key2=${key[1]}
	type1=${type[0]}
	type2=${type[1]}
	chrnum1=${chrnum[0]}
	chrnum2=${chrnum[1]}
else
	prefix1=${abbr[0]}
	gff3file1=${gff3[0]}
	latin1=${sample[0]}
	protein1=${protein[0]}
	cds1=${cds[0]}
	key1=${key[0]}
	type1=${type[0]}
	chrnum1=${chrnum[0]}
	prefix2=${abbr[0]}
	gff3file2=${gff3[0]}
	latin2=${sample[0]}
	protein2=${protein[0]}
	cds2=${cds[0]}
	key2=${key[0]}
	type2=${type[0]}
	chrnum2=${chrnum[0]}
fi

##获取输入文件

#test.cds #每个基因最长的转录本的DNA序列
#test.pep #每个基因最长的蛋白序列


#从gff3文件获取bed,用法：getbed gff3file output前缀 第三列的type 第9列的前缀字符
function getbed(){
	if [ $# -lt 2 ];then
		echo "usage:
		getbed inputgff3 outputprefix type key
		
		inputgff3 can be gff3 or gff3.gz .(Required)
		
		outputprefix is the output file prefix.Preferably a 3-character abbr.(Required)
		
		type is gfffile the 3rd cloumn string (Value:mRNA ,gene or other,Default:mRNA) 
		
		key is gfffile the Prefix for column 9.(Value:ID or other,Default:ID)
		"
		return 1
		exit
	else
		inputgff3=$1
		prefix=$2
		if [ $# -eq 3 ];then
			type=$3
		elif [ $# -eq 4 ];then
			type=$3
			key=$4
		else
			echo "Usage: -h /-help "
		fi
	fi
	python -m jcvi.formats.gff bed --type=${type:=mRNA} --key=${key:=ID} ${inputgff3} -o ${prefix}.bed
	python -m jcvi.formats.bed uniq ${prefix}.bed
	mv ${prefix}.uniq.bed ${prefix}.bed
}


function getcds(){
	if [ $# -lt 2 ];then
		echo "usage:
		getbed input_cdsfa prefix 
		input_cdsfa can be fa or fa.gz .(Required)
		prefix is the input bed file prefix.Preferably a 3-character abbr.(Required)
		"
		return 1
		exit
	else
		input_cdsfa=$1
		prefix=$2
	fi
	seqkit grep -f <(cut -f4 ${prefix}.bed) ${input_cdsfa} | seqkit seq -i >${prefix}.cds
}

function getpep(){
	if [ $# -lt 2 ];then
		echo "usage:
		getbed input_proteinfa prefix 
		input_proteinfa can be fa or fa.gz .(Required)
		prefix is the input bed file prefix.Preferably a 3-character abbr.(Required)
		"
		return 1
		exit
	else
		input_proteinfa=$1
		prefix=$2
	fi
	seqkit grep -f <(cut -f4 ${prefix}.bed) ${input_proteinfa}  | seqkit seq -i >${prefix}.pep
}


function getcoline(){
	if [ $# -lt 2 ];then
		echo "usage:
		getcoline abbr1 abbr2 
		"
		return 1
		exit
	fi
	species1="$1"
	species2="$2"
	## 运行代码
	python -m jcvi.compara.catalog ortholog --dbtype prot --no_strip_names $species1 $species2
	python -m jcvi.compara.synteny screen --minspan=30 --simple $species1.$species2.anchors $species1.$species2.anchors.new
}

function VisualColine(){
	if [ $# -lt 2 ];then
		echo "usage:
		VisualColine abbr1 abbr2 chrnum1 chrnum2
		"
		return 1
		exit
	fi
	abbr1=$1
	abbr2=$2
	chrnum1=$3
	chrnum2=$4
	##可视化
	cat $abbr1.bed|cut -f1|sort |uniq |head -$chrnum1 |rev|cut -d " " -f1|rev >$abbr1.id
	cat $abbr1.id|awk 'BEGIN{c=0;} {for(i=1;i<=NF;i++) {num[c,i] = $i;} c++;} END{ for(i=1;i<=NF;i++){str=""; for(j=0;j<NR;j++){ if(j>0){str = str","} str= str""num[j,i]}printf("%s\n", str)} }' >$abbr1.ids
	cat $abbr2.bed|cut -f1|sort |uniq |head -$chrnum2 |rev|cut -d " " -f1|rev >$abbr2.id
	cat $abbr2.id|awk 'BEGIN{c=0;} {for(i=1;i<=NF;i++) {num[c,i] = $i;} c++;} END{ for(i=1;i<=NF;i++){str=""; for(j=0;j<NR;j++){ if(j>0){str = str","} str= str""num[j,i]}printf("%s\n", str)} }' >$abbr2.ids
	cat $abbr1.ids $abbr2.ids >seqids
	
	# 设置颜色，长宽等
	echo "
	# y, xstart, xend, rotation, color, label, va, bed
	 .6,    .1,    .8,    0,    red,    $latin1,    top,     $abbr1.bed
	 .4,    .1,    .8,    0,    blue,    $latin2,    top,    $abbr2.bed
	# edges
	e, 0, 1, $abbr1.$abbr2.anchors.simple
	" >layout
	#生成共线性图片，很可能运行失败。注意：修改layout的细节就好，python3对文件要求比较严格。
	python -m jcvi.graphics.karyotype seqids layout
	echo "karyotype.pdf is the coline picture!"
}

#准备kaks和4DTv的文件
function prepareResult(){
	if [ $# -lt 2 ];then
		echo "Usage:prepareResult abbr1 abbr2 threads 
		threads should be a number 32 or 64 or other
		"
	elif [ $# -eq 2 ];then
		abbr1=$1
		abbr2=$2
	else
		abbr1=$1
		abbr2=$2
		thread=$3
	fi
	#判断旧版本的输出目录是否存在
	if [ -d ${abbr1}_${abbr2}.result_dir ];then
		read -p "The fold ${abbr1}_${abbr2}.result_dir is exist.Delete the old version ?(Y/N):" -n 1 answer
		case $answer in
			Y|y)
				echo -e "\n ok!Delete the fold ${abbr1}_${abbr2}.result_dir！"
				;;
			N|n)
				echo -e "\n The old version will be rename ${abbr1}_${abbr2}.result_dir.old!"
				mv ${abbr1}_${abbr2}.result_dir ${abbr1}_${abbr2}.result_dir.old
				;;
			*)
				echo -e "\n The old version will be rename ${abbr1}_${abbr2}.result_dir.old!"
				mv ${abbr1}_${abbr2}.result_dir ${abbr1}_${abbr2}.result_dir.old
				;;
		esac
	fi
	echo ${thread:=32} >proc
	cat ${abbr1}.${abbr2}.anchors|grep -v ^#|cut -f 1-2 >${abbr1}_${abbr2}.homolog
	cat ${abbr1}.cds ${abbr2}.cds >${abbr1}_${abbr2}.cds
	cat ${abbr1}.pep ${abbr2}.pep >${abbr1}_${abbr2}.pep
	#此程序需要依赖较多
	ParaAT.pl -h ${abbr1}_${abbr2}.homolog -n ${abbr1}_${abbr2}.cds -a ${abbr1}_${abbr2}.pep -p proc -m mafft -f axt -g -k -o ${abbr1}_${abbr2}.result_dir
}

#输出结果在result_dir目录
function getkaks(){
	if [ $# -lt 2 ];then
		echo "usage:
		getkaks abbr1 abbr2
		"
		return 1
		exit
	fi
	abbr1=$1
	abbr2=$2
	#判断是否存在result_dir,不存在则需要先运行prepareResult
	if [ ! -d ${abbr1}_${abbr2}.result_dir ];then
		echo "请先运行prepareResult函数，以生成准备文件!"
		return 1
		exit
	fi
	#合并所有同源基因对的kaks值
	find ${abbr1}_${abbr2}.result_dir -name "*.axt.kaks"|xargs cat | cut -f 1,3,4,5 | grep -v 'Sequence'|sort|uniq >${abbr1}_${abbr2}.all-kaks.results
	cat ${abbr1}_${abbr2}.all-kaks.results|sed '1i\Seq\tKa\tKs\tKa/Ks'|tr "\t" "," >${abbr1}_${abbr2}.all-kaks.csv
}

function get4DTv(){
	if [ $# -lt 2 ];then
		echo "usage:
		get4DTv abbr1 abbr2
		"
		return 1
		exit
	fi
	abbr1=$1
	abbr2=$2
	#判断是否存在result_dir,不存在则需要先运行prepareResult
	if [ ! -d ${abbr1}_${abbr2}.result_dir ];then
		echo "请先运行prepareResult函数，以生成准备文件!"
		return 1
		exit
	fi
	##获取4DTv的值
	#将多行axt文件转换成单行
	for i in `find ${abbr1}_${abbr2}.result_dir -name "*.axt"`;do axt2one-line.py $i ${i}.one-line;done
	#使用calculate_4DTV_correction.pl脚本计算4dtv值
	find ${abbr1}_${abbr2}.result_dir -name "*.axt.one-line"|while read id;do calculate_4DTV_correction.pl $id >${id%%one-line}4dtv;done
	#合并所有同源基因对的4dtv
	find ${abbr1}_${abbr2}.result_dir -name "*.4dtv" |xargs cat| cut -f 1,3| grep -v '4dtv_raw'|sort|uniq >${abbr1}_${abbr2}.all-4dtv.results
	cat ${abbr1}_${abbr2}.all-4dtv.results| sed '1i\Seq\t4dtv_corrected'|tr "\t" "," >${abbr1}_${abbr2}.all-4dtv.csv
}

function getkaks4DTv(){
	if [ $# -lt 2 ];then
		echo "usage:
		getkaks4DTv abbr1 abbr2
		"
		return 1
		exit
	fi
	abbr1=$1
	abbr2=$2
	#判断是否存在result_dir,不存在则需要先运行prepareResult
	if [ ! -e ${abbr1}_${abbr2}.all-4dtv.results ];then
		echo "请先运行get4DTv函数，以生成4DTv!"
		return 1
		exit
	fi
	if [ ! -e ${abbr1}_${abbr2}.all-kaks.results ];then
		echo "请先运行getkaks函数，以生成kaks!"
		return 1
		exit
	fi
	#将kaks结果和4Dtv结果合并
	join -a 1 -a 2 -1 1 -2 1 ${abbr1}_${abbr2}.all-4dtv.results ${abbr1}_${abbr2}.all-kaks.results |sed '1i\Seq 4dtv_corrected Ka Ks Ka/Ks' >${abbr1}_${abbr2}.all-results.txt
	#给结果文件添加标题
	cat ${abbr1}_${abbr2}.all-results.txt|sed 's/ /,/g'   >${abbr1}_${abbr2}.kaks4DTv.csv
}






