#!/bin/bash


#####
#此程序运行于河南农大农学院的大服务器上
#最终会输出表达矩阵 两个csv 文件
#####

##处理RNA-seq的数据
###仅限于处理杂交种和亲本的RNA-seq
:<<EOF
多行注释
使用时需要修改的地方：
bef数组中的元素需要修改
EOF

##获取样本的无重复的名称
#ls *.gz |awk 'NR%2'|awk -F "." '{print $1}'  >sample
#bef=($(awk ‘{print $1}’ sample))

##全局变量3个，基因组文件不可随意修改，需要同时修改多个注释文件
bg1='.R1.fq.gz'
bg2='.R2.fq.gz'
#bef=(CO1_2 CO1_3 STTM166_1 STTM166_3)
#group_name='166'

bef=(SSSL-IM1 SSSL-IM1 SSSL-IM2 SSSL-IM2 SSSL-SM1 SSSL-SM1 SSSL-SM2 SSSL-SM2 SSSL-SM3 SSSL-SM3 SSSL-SPM1 SSSL-SPM1 SSSL-SPM2 SSSL-SPM2 SSSL-SPM3 SSSL-SPM3 X178-IM1 X178-IM1 X178-IM2 X178-IM2 X178-SM1 X178-SM1 X178-SM2 X178-SM2 X178-SM3 X178-SM3 X178-SPM1 X178-SPM1 X178-SPM2 X178-SPM2 X178-SPM3 X178-SPM3)
###group_name 指的是这一组材料的组名，当有多个组时，方便识别输出文件是哪个组的。
group_name='wangll'
genome="/public/home/tang/chaim/maize/Zea_mays.B73_RefGen_v4.42.fa"
gtf="/public/home/tang/chaim/maize/Zea_mays.B73_RefGen_v4.42.gtf"
##自定义函数用于输出程序运行的时间点  用法：echo start softname
echotime (){
if [ $1 = start ]
then
	echo "Start $2 at $(date +%x_%X)"
elif [ $1 = end ]
then
	echo "End $2 at $(date +%x_%X)"
else
	echo "input 参数错误！"
fi
}
echotime start analysis_RNA-seq

#复制下机数据到新的文件夹data ,尽量避免操作原始文件
#find ./Cleandata -name '*fq.gz'|xargs -i cp {} ./data
 #在存放clean的文件夹运行配套的重命名脚本rename.bash
 #bash rename.bash 
 #1.质控 
##获取输入的数组的长度
sample_num=${#bef[@]}

mkdir QC
for ((i=0;i<$sample_num;i++));
do
inA1=${bef[$i]}$bg1;
inA2=${bef[$i]}$bg2;
out1=${bef[$i]}"paired-R1.fq.gz";
out2=${bef[$i]}"paired-R2.fq.gz";
unpaired1=${bef[$i]}"unpaired-R1.fq.gz";
unpaired2=${bef[$i]}"unpaired-R2.fq.gz";
fastp --thread 16 -i $inA1 -o $out1 -I $inA2 -O $out2 --html ./QC/${bef[$i]}".html" --json ./QC/${bef[$i]}".json" 2>./QC/${bef[$i]}fastp_out ;
echo success end  QC_${bef[$i]} fastp at `date`;
done


#2.比对 （此处默认已经生成hisat2的索引文件,不再运行生成索引这一步）
echotime start hisat2;

#/public/home/tang/chaim/soft/hisat2/hisat2-2.1.0/extract_exons.py $gtf >genome.exon
#/public/home/tang/chaim/soft/hisat2/hisat2-2.1.0/extract_splice_sites.py $gtf >genome.ss
#/disks/backup/chaim/soft/hisat2/hisat2_extract_snps_haplotypes_VCF.py Zea_mays.B73_RefGen_v4.42.fa zea_mays.vcf genome.snp
##2.1构建索引文件，非常耗时，至少1-2h,多则1-2天
#hisat2-build -p 48 ${genome} --ss genome.ss --exon genome.exon /public/home/tang/chaim/maize/genome_tran 
#内存>200G的服务器，可以增加snp信息，命令如下。
#hisat2-build -p 48 ${genome} --snp genome.snp --ss genome.ss --exon genome.exon /disks/backup/chaim/maize/genome_tran &

##2.2开始比对
for ((i=0;i<$sample_num;i++));
	do	
	out1=${bef[$i]}"paired-R1.fq.gz";
	out2=${bef[$i]}"paired-R2.fq.gz";
	
	hisat2 -x /public/home/tang/chaim/maize/genome_tran -p 16 -1 $out1 -2 $out2 -S ${bef[$i]}".map.sam" --dta-cufflinks --novel-splicesite-outfile ${bef[$i]}".nsplice" 2>${bef[$i]}"hisat2_out"
		echo "complete all hisat2 . 完成所有的比对"
	done
echo end hisat2;

#####等第2步执行完成之后，才可执行第3步。
#3.排序
for ((i=0;i<$sample_num;i++));
do
samtools sort -@ 48 -O bam -o ${bef[$i]}".map.bam" ${bef[$i]}".map.sam" 2>${bef[$i]}"samtool_out"
#构建索引文件，方便使用IGV查看拼接结果
samtools index -@ 48 ${bef[$i]}".map.bam" ${bef[$i]}".map.bai" &

stringtie ${bef[$i]}".map.bam" -G $gtf -p 24 -o ${bef[$i]}".gtf" 2>${bef[$i]}"stringtie_first"
done

##3.1合并先前所有的样本数据gtf为merge.gtf

##用于遍历原始的材料名数组，生成原始材料名+”.gtf“的所有元素
for ele in ${bef[@]};
do
 echo $ele.gtf 
 str="$str ${ele}.gtf"
done
echo ${str}

stringtie --merge -G $gtf -p 16 -o merged_$group_name.gtf ${str} 2>stringtie_merge_$group_name
##4. 组装定量
for ((i=0;i<$sample_num;i++));
do
	#判断输出文件夹是否存在，不存在则新建。
	if [-d ${bef[$i]}"_out"]
	then
		echo "文件已存在，不再新建"
	else
		mkdir ${bef[$i]}"_out"
	fi
	##第一轮组装定量
	stringtie ${bef[$i]}".map.bam" -G merged_$group_name.gtf -p 8 -b ${bef[$i]}"_out" -e -o ${bef[$i]}"-st.gtf"
done

##5.数据整合转换
for ((i=0;i<$sample_num;i++));
do
	echo "${bef[$i]}-st ./${bef[$i]}-st.gtf" >>gtf2_$group_name
done
##6.输出计算矩阵，gene_count_matrix.csv  基因表达矩阵    transcript_count_matrix.csv  转录本表达矩阵
python2.7 /public/home/tang/chaim/soft/stringtie-2.0.6/prepDE.py -i gtf2_$group_name
mv gene_count_matrix.csv gene_$group_name.gene.csv 
mv transcript_count_matrix.csv $group_name.trans.csv
echo "complete all.next step run in R-studio"
echo "output file is gene_$group_name.csv & trans_$group_name.csv "
echotime end ALL-analysis

