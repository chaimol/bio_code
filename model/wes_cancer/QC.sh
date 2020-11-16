#!/bin/sh
<<!
##Fastqc查看数据质量
cat config1 | while read id
do
	fastqc --outdir ./3.qc/raw_qc/ --threads 16 ./1.raw_fq/${id}*.fastq.gz >> ./3.qc/raw_qc/${id}_fastqc.log 2>&1 
done 
#multiqc讲fastqc生成的多个报告文件整合成为一个报告文件

!
multiqc  ./3.qc/raw_qc/*zip  -o ./3.qc/raw_qc/multiqc

## trim_galore.sh
##使用trim_galore质控数据，去除接头
cat config | while read id
do
	fq1=./1.raw_fq/${id}_1.fastq.gz
	fq2=./1.raw_fq/${id}_2.fastq.gz
	trim_galore  --paired -q 28 --phred33 --length 30 --stringency 3 --gzip --cores 8 -o ./2.clean_fq  $fq1  $fq2 >> ./2.clean_fq/${id}_trim.log 2>&1
done
#nohup bash trim_galore.sh &

##再次使用fastqc检测质控之后的数据质量
cat config | while read id
do
	fastqc --outdir ./3.qc/clean_qc/ --threads 16 ./2.clean_fq/${id}*.fq.gz >> ./3.qc/clean_qc/${id}_fastqc.log 2>&1 
done 

multiqc  ./3.qc/clean_qc/*zip  -o ./3.qc/clean_qc/multiqc


###############################质控完成之后，进行比对###########################################################
##先使用bwa bwtsw 模式构建索引文件，这个只用构建一次，以后只要不换参考基因组，就可以一直用，不需要再次构建。
##### conda activate wes
##### cd ~/wes_cancer/data
##### gunzip Homo_sapiens_assembly38.fasta.gz
##### time bwa index -a bwtsw -p gatk_hg38 ~/wes_cancer/data/Homo_sapiens_assembly38.fasta
##### INDEX=~/wes_cancer/data/gatk_hg38

#此处我使用的是ensembl下载的参考基因组，而且自己构建索引在指定地点，故修改了此处的代码。
#genome=/share/home/chaimao1/human/hg38/hg38.fa
#cd /share/home/chaimao1/human/hg38/index/bwa/
#bwa index -a bwtsw -p hg38 $genome 
<<!
cd ~/wes_cancer/project
## bwa.sh 开始比对
INDEX=/share/home/chaimao1/human/hg38/index/bwa/hg38
##直接结合samtools，把sam文件转换成bam文件
cat config | while read id
do
	echo "start bwa for ${id}" `date`
	fq1=./2.clean_fq/${id}_1_val_1.fq.gz
	fq2=./2.clean_fq/${id}_2_val_2.fq.gz
	bwa mem -M -t 16 -R "@RG\tID:${id}\tSM:${id}\tLB:WXS\tPL:Illumina" ${INDEX} ${fq1} ${fq2} | samtools sort -@ 10 -m 10G  -o  ./4.align/${id}.bam -
	echo "end bwa for ${id}" `date`
done


#####提取指定的文件的指定染色体

##samtools index case5_germline.bam
##samtools view -h case5_germline.bam chr12 | samtools view -Sb - > small.bam
##samtools index small.bam

######可以进一步提取指定区间内的bam文件
samtools index case5_germline.bam
samtools view -h case5_germline.bam chr1:25151000-30151000 | samtools view -Sb - > chr1_small.bam
samtools index chr1_small.bam
!
