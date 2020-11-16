#!/bin/sh
###肿瘤外显子实战项目
##########搭建开发环境############################
## 首先在用户的主目录下创建 wes_cancer 文件夹作为工作目录
mkdir ~/wes_cancer
cd ~/wes_cancer
## 在 ~/wes_cancer 中创建 biosoft project data 三个文件夹
## biosoft 存放软件安装包
## project 存放分析过程产生的文件
## data 存放数据库文件
mkdir biosoft project data
cd project
## 在 project 文件夹中创建若干个文件夹，分别存放每一步产生的文件
mkdir -p 0.sra 1.raw_fq 2.clean_fq 3.qc/{raw_qc,clean_qc} 4.align/{qualimap,flagstat,stats} 5.gatk/gvcf 6.mutect 7.annotation/{vep,annovar,funcotator,snpeff} 8.cnv/{gatk,cnvkit,gistic,facet} 9.pyclone 10.signature

####################安装miniconda,安装软件##################
cd ~/wes_cancer/biosoft
wget https://mirrors.tuna.tsinghua.edu.cn/anaconda/miniconda/Miniconda3-latest-Linux-x86_64.sh
##下载完成后，安装miniconda。一直按enter或者输入yes就可以。
bash Miniconda3-latest-Linux-x86_64.sh

##重新载入环境变量文件
source ~/.bashrc
##在国内的修改conda的源为清华源
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda/
conda config --set show_channel_urls yes

## 新建小环境 wes
conda create -n wes python=3
## 激活 wes 小环境
conda activate wes
## 安装必要的生信软件
conda install -y sra-tools fastqc trim-galore multiqc bwa samtools gnuplot qualimap subread vcftools bedtools cnvkit 
conda install -y -c hcc aspera-cli=3.7.7
##gatk软件需要自行安装

###############下载必要的文件，包括基因组、索引、变异信息、注释信息##########################################
##使用nohup wget url & 格式下载，让数据在后台自动下载。下载完成后，检查文件的完整性。如果服务器网络不好，或者无法访问外网，那就本地电脑用下载工具下载，再上传。
cd ~/wes_cancer/data/
nohup wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta.gz &
## gatk
nohup wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/dbsnp_146.hg38.vcf.gz & 
nohup wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/dbsnp_146.hg38.vcf.gz.tbi & 
nohup wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz & 
nohup wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi & 
nohup wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta.gz & 
nohup wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta.fai & 
nohup wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.dict & 
nohup wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz & 
nohup wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi & 
nohup wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/funcotator/funcotator_dataSources.v1.6.20190124s.tar.gz &

## bed
wget ftp://ftp.ncbi.nlm.nih.gov/pub/CCDS/current_human/CCDS.current.txt
cat CCDS.current.txt | grep  "Public" | perl -alne '{/\[(.*?)\]/;next unless $1;$gene=$F[2];$exons=$1;$exons=~s/\s//g;$exons=~s/-/\t/g;print "$F[0]\t$_\t$gene" foreach split/,/,$exons;}'|sort -u |bedtools sort -i |awk '{if($3>$2) print "chr"$0"\t0\t+"}'  > hg38.exon.bed

############安装Aspera,下载原始数据###################
cd  ~/wes_cancer/biosoft
wget -c  https://download.asperasoft.com/download/sw/connect/3.8.1/ibm-aspera-connect-3.8.1.161274-linux-g2.12-64.tar.gz
tar -zxvf ibm-aspera-connect-3.8.1.161274-linux-g2.12-64.tar.gz
bash ibm-aspera-connect-3.8.1.161274-linux-g2.12-64.sh
##添加环境变量，更新环境变量
echo 'export PATH=~/.aspera/connect/bin/:$PATH' >> ~/.bashrc
source ~/.bashrc


#######测序数据的下载##############################
#####从这个链接地址https://sra-explorer.info/获取下载测序数据的链接，搜索SRP070662,找到前面是exome的外显子数据，勾选要下载的数据的复选框，点击Add 30 datasets,选择FastQ Download,选择ASpera command for download FastQ files.最后选择下面自动生成的地址，复制或者下载到本地，拷贝到服务器上，写成sh脚本就可以下载。
##下载时候遇到的错误原因 1.未打开33001端口 2.我们的服务器集群的其他node节点未联网，只有登录节点联网，只能在登录节点下载。
###root权限下打开33001端口的命令：（必须是root，非root就放弃吧）
###iptables -I INPUT -p udp --dport 33001 -j ACCEPT
###iptables -I OUTPUT -p udp --dport 33001 -j ACCEPT

##下载之前，可以先使用awk和grep处理命令。
#下面的文件结尾是_1.fastq.gz和_2.fastq.gz的是我们需要的文件。
#将原始下载文献保存为name.txt
cat name.txt|grep -E '_1.fastq.gz|_2.fastq.gz' >name2.txt
#再把name2的内容粘贴回原始下载脚本即可。
#!/usr/bin/env bash
cd ~/wes_cancer/project/1.raw_fq
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/008/SRR3182418/SRR3182418.fastq.gz . && mv SRR3182418.fastq.gz SRR3182418_exome_sequencing_of_case5_germline.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/008/SRR3182418/SRR3182418_1.fastq.gz . && mv SRR3182418_1.fastq.gz SRR3182418_exome_sequencing_of_case5_germline_1.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/008/SRR3182418/SRR3182418_2.fastq.gz . && mv SRR3182418_2.fastq.gz SRR3182418_exome_sequencing_of_case5_germline_2.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/009/SRR3182419/SRR3182419.fastq.gz . && mv SRR3182419.fastq.gz SRR3182419_exome_sequencing_of_case2_germline.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/009/SRR3182419/SRR3182419_1.fastq.gz . && mv SRR3182419_1.fastq.gz SRR3182419_exome_sequencing_of_case2_germline_1.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/009/SRR3182419/SRR3182419_2.fastq.gz . && mv SRR3182419_2.fastq.gz SRR3182419_exome_sequencing_of_case2_germline_2.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/000/SRR3182420/SRR3182420.fastq.gz . && mv SRR3182420.fastq.gz SRR3182420_exome_sequencing_of_case4_germline.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/000/SRR3182420/SRR3182420_1.fastq.gz . && mv SRR3182420_1.fastq.gz SRR3182420_exome_sequencing_of_case4_germline_1.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/000/SRR3182420/SRR3182420_2.fastq.gz . && mv SRR3182420_2.fastq.gz SRR3182420_exome_sequencing_of_case4_germline_2.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/003/SRR3182423/SRR3182423.fastq.gz . && mv SRR3182423.fastq.gz SRR3182423_exome_sequencing_of_case1_germline.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/003/SRR3182423/SRR3182423_1.fastq.gz . && mv SRR3182423_1.fastq.gz SRR3182423_exome_sequencing_of_case1_germline_1.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/003/SRR3182423/SRR3182423_2.fastq.gz . && mv SRR3182423_2.fastq.gz SRR3182423_exome_sequencing_of_case1_germline_2.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/001/SRR3182421/SRR3182421.fastq.gz . && mv SRR3182421.fastq.gz SRR3182421_exome_sequencing_of_case3_germline.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/001/SRR3182421/SRR3182421_1.fastq.gz . && mv SRR3182421_1.fastq.gz SRR3182421_exome_sequencing_of_case3_germline_1.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/001/SRR3182421/SRR3182421_2.fastq.gz . && mv SRR3182421_2.fastq.gz SRR3182421_exome_sequencing_of_case3_germline_2.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/002/SRR3182422/SRR3182422.fastq.gz . && mv SRR3182422.fastq.gz SRR3182422_exome_sequencing_of_case6_germline.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/002/SRR3182422/SRR3182422_1.fastq.gz . && mv SRR3182422_1.fastq.gz SRR3182422_exome_sequencing_of_case6_germline_1.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/002/SRR3182422/SRR3182422_2.fastq.gz . && mv SRR3182422_2.fastq.gz SRR3182422_exome_sequencing_of_case6_germline_2.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/004/SRR3182424/SRR3182424.fastq.gz . && mv SRR3182424.fastq.gz SRR3182424_exome_sequencing_of_case4_biorep_A.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/004/SRR3182424/SRR3182424_1.fastq.gz . && mv SRR3182424_1.fastq.gz SRR3182424_exome_sequencing_of_case4_biorep_A_1.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/004/SRR3182424/SRR3182424_2.fastq.gz . && mv SRR3182424_2.fastq.gz SRR3182424_exome_sequencing_of_case4_biorep_A_2.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/005/SRR3182425/SRR3182425.fastq.gz . && mv SRR3182425.fastq.gz SRR3182425_exome_sequencing_of_case4_biorep_B_techrep_1.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/005/SRR3182425/SRR3182425_1.fastq.gz . && mv SRR3182425_1.fastq.gz SRR3182425_exome_sequencing_of_case4_biorep_B_techrep_1_1.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/005/SRR3182425/SRR3182425_2.fastq.gz . && mv SRR3182425_2.fastq.gz SRR3182425_exome_sequencing_of_case4_biorep_B_techrep_1_2.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/009/SRR3182429/SRR3182429.fastq.gz . && mv SRR3182429.fastq.gz SRR3182429_exome_sequencing_of_case3_biorep_C_techrep_1.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/009/SRR3182429/SRR3182429_1.fastq.gz . && mv SRR3182429_1.fastq.gz SRR3182429_exome_sequencing_of_case3_biorep_C_techrep_1_1.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/009/SRR3182429/SRR3182429_2.fastq.gz . && mv SRR3182429_2.fastq.gz SRR3182429_exome_sequencing_of_case3_biorep_C_techrep_1_2.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/000/SRR3182430/SRR3182430.fastq.gz . && mv SRR3182430.fastq.gz SRR3182430_exome_sequencing_of_case6_biorep_A_techrep_1.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/000/SRR3182430/SRR3182430_1.fastq.gz . && mv SRR3182430_1.fastq.gz SRR3182430_exome_sequencing_of_case6_biorep_A_techrep_1_1.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/000/SRR3182430/SRR3182430_2.fastq.gz . && mv SRR3182430_2.fastq.gz SRR3182430_exome_sequencing_of_case6_biorep_A_techrep_1_2.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/007/SRR3182427/SRR3182427.fastq.gz . && mv SRR3182427.fastq.gz SRR3182427_exome_sequencing_of_case3_biorep_A.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/007/SRR3182427/SRR3182427_1.fastq.gz . && mv SRR3182427_1.fastq.gz SRR3182427_exome_sequencing_of_case3_biorep_A_1.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/007/SRR3182427/SRR3182427_2.fastq.gz . && mv SRR3182427_2.fastq.gz SRR3182427_exome_sequencing_of_case3_biorep_A_2.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/008/SRR3182428/SRR3182428.fastq.gz . && mv SRR3182428.fastq.gz SRR3182428_exome_sequencing_of_case3_biorep_B.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/008/SRR3182428/SRR3182428_1.fastq.gz . && mv SRR3182428_1.fastq.gz SRR3182428_exome_sequencing_of_case3_biorep_B_1.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/008/SRR3182428/SRR3182428_2.fastq.gz . && mv SRR3182428_2.fastq.gz SRR3182428_exome_sequencing_of_case3_biorep_B_2.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/001/SRR3182431/SRR3182431.fastq.gz . && mv SRR3182431.fastq.gz SRR3182431_exome_sequencing_of_case6_biorep_B.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/001/SRR3182431/SRR3182431_1.fastq.gz . && mv SRR3182431_1.fastq.gz SRR3182431_exome_sequencing_of_case6_biorep_B_1.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/001/SRR3182431/SRR3182431_2.fastq.gz . && mv SRR3182431_2.fastq.gz SRR3182431_exome_sequencing_of_case6_biorep_B_2.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/006/SRR3182426/SRR3182426.fastq.gz . && mv SRR3182426.fastq.gz SRR3182426_exome_sequencing_of_case4_biorep_C.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/006/SRR3182426/SRR3182426_1.fastq.gz . && mv SRR3182426_1.fastq.gz SRR3182426_exome_sequencing_of_case4_biorep_C_1.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/006/SRR3182426/SRR3182426_2.fastq.gz . && mv SRR3182426_2.fastq.gz SRR3182426_exome_sequencing_of_case4_biorep_C_2.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/003/SRR3182433/SRR3182433.fastq.gz . && mv SRR3182433.fastq.gz SRR3182433_exome_sequencing_of_case1_biorep_A_techrep_1.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/003/SRR3182433/SRR3182433_1.fastq.gz . && mv SRR3182433_1.fastq.gz SRR3182433_exome_sequencing_of_case1_biorep_A_techrep_1_1.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/003/SRR3182433/SRR3182433_2.fastq.gz . && mv SRR3182433_2.fastq.gz SRR3182433_exome_sequencing_of_case1_biorep_A_techrep_1_2.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/004/SRR3182434/SRR3182434.fastq.gz . && mv SRR3182434.fastq.gz SRR3182434_exome_sequencing_of_case1_biorep_B.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/004/SRR3182434/SRR3182434_1.fastq.gz . && mv SRR3182434_1.fastq.gz SRR3182434_exome_sequencing_of_case1_biorep_B_1.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/004/SRR3182434/SRR3182434_2.fastq.gz . && mv SRR3182434_2.fastq.gz SRR3182434_exome_sequencing_of_case1_biorep_B_2.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/005/SRR3182435/SRR3182435.fastq.gz . && mv SRR3182435.fastq.gz SRR3182435_exome_sequencing_of_case1_biorep_C.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/005/SRR3182435/SRR3182435_1.fastq.gz . && mv SRR3182435_1.fastq.gz SRR3182435_exome_sequencing_of_case1_biorep_C_1.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/005/SRR3182435/SRR3182435_2.fastq.gz . && mv SRR3182435_2.fastq.gz SRR3182435_exome_sequencing_of_case1_biorep_C_2.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/006/SRR3182436/SRR3182436.fastq.gz . && mv SRR3182436.fastq.gz SRR3182436_exome_sequencing_of_case2_biorep_A.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/006/SRR3182436/SRR3182436_1.fastq.gz . && mv SRR3182436_1.fastq.gz SRR3182436_exome_sequencing_of_case2_biorep_A_1.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/006/SRR3182436/SRR3182436_2.fastq.gz . && mv SRR3182436_2.fastq.gz SRR3182436_exome_sequencing_of_case2_biorep_A_2.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/007/SRR3182437/SRR3182437.fastq.gz . && mv SRR3182437.fastq.gz SRR3182437_exome_sequencing_of_case2_biorep_B_techrep_1.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/007/SRR3182437/SRR3182437_1.fastq.gz . && mv SRR3182437_1.fastq.gz SRR3182437_exome_sequencing_of_case2_biorep_B_techrep_1_1.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/007/SRR3182437/SRR3182437_2.fastq.gz . && mv SRR3182437_2.fastq.gz SRR3182437_exome_sequencing_of_case2_biorep_B_techrep_1_2.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/008/SRR3182438/SRR3182438.fastq.gz . && mv SRR3182438.fastq.gz SRR3182438_exome_sequencing_of_case2_biorep_C.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/008/SRR3182438/SRR3182438_1.fastq.gz . && mv SRR3182438_1.fastq.gz SRR3182438_exome_sequencing_of_case2_biorep_C_1.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/008/SRR3182438/SRR3182438_2.fastq.gz . && mv SRR3182438_2.fastq.gz SRR3182438_exome_sequencing_of_case2_biorep_C_2.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/009/SRR3182439/SRR3182439.fastq.gz . && mv SRR3182439.fastq.gz SRR3182439_exome_sequencing_of_case5_biorep_A.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/009/SRR3182439/SRR3182439_1.fastq.gz . && mv SRR3182439_1.fastq.gz SRR3182439_exome_sequencing_of_case5_biorep_A_1.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/009/SRR3182439/SRR3182439_2.fastq.gz . && mv SRR3182439_2.fastq.gz SRR3182439_exome_sequencing_of_case5_biorep_A_2.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/002/SRR3182432/SRR3182432.fastq.gz . && mv SRR3182432.fastq.gz SRR3182432_exome_sequencing_of_case6_biorep_C.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/002/SRR3182432/SRR3182432_1.fastq.gz . && mv SRR3182432_1.fastq.gz SRR3182432_exome_sequencing_of_case6_biorep_C_1.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/002/SRR3182432/SRR3182432_2.fastq.gz . && mv SRR3182432_2.fastq.gz SRR3182432_exome_sequencing_of_case6_biorep_C_2.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/000/SRR3182440/SRR3182440.fastq.gz . && mv SRR3182440.fastq.gz SRR3182440_exome_sequencing_of_case5_biorep_B_techrep_1.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/000/SRR3182440/SRR3182440_1.fastq.gz . && mv SRR3182440_1.fastq.gz SRR3182440_exome_sequencing_of_case5_biorep_B_techrep_1_1.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/000/SRR3182440/SRR3182440_2.fastq.gz . && mv SRR3182440_2.fastq.gz SRR3182440_exome_sequencing_of_case5_biorep_B_techrep_1_2.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/001/SRR3182441/SRR3182441.fastq.gz . && mv SRR3182441.fastq.gz SRR3182441_exome_sequencing_of_case4_techrep_2.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/001/SRR3182441/SRR3182441_1.fastq.gz . && mv SRR3182441_1.fastq.gz SRR3182441_exome_sequencing_of_case4_techrep_2_1.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/001/SRR3182441/SRR3182441_2.fastq.gz . && mv SRR3182441_2.fastq.gz SRR3182441_exome_sequencing_of_case4_techrep_2_2.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/003/SRR3182443/SRR3182443.fastq.gz . && mv SRR3182443.fastq.gz SRR3182443_exome_sequencing_of_case6_techrep_2.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/003/SRR3182443/SRR3182443_1.fastq.gz . && mv SRR3182443_1.fastq.gz SRR3182443_exome_sequencing_of_case6_techrep_2_1.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/003/SRR3182443/SRR3182443_2.fastq.gz . && mv SRR3182443_2.fastq.gz SRR3182443_exome_sequencing_of_case6_techrep_2_2.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/002/SRR3182442/SRR3182442.fastq.gz . && mv SRR3182442.fastq.gz SRR3182442_exome_sequencing_of_case3_techrep_2.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/002/SRR3182442/SRR3182442_1.fastq.gz . && mv SRR3182442_1.fastq.gz SRR3182442_exome_sequencing_of_case3_techrep_2_1.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/002/SRR3182442/SRR3182442_2.fastq.gz . && mv SRR3182442_2.fastq.gz SRR3182442_exome_sequencing_of_case3_techrep_2_2.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/004/SRR3182444/SRR3182444.fastq.gz . && mv SRR3182444.fastq.gz SRR3182444_exome_sequencing_of_case1_techrep_2.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/004/SRR3182444/SRR3182444_1.fastq.gz . && mv SRR3182444_1.fastq.gz SRR3182444_exome_sequencing_of_case1_techrep_2_1.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/004/SRR3182444/SRR3182444_2.fastq.gz . && mv SRR3182444_2.fastq.gz SRR3182444_exome_sequencing_of_case1_techrep_2_2.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/005/SRR3182445/SRR3182445.fastq.gz . && mv SRR3182445.fastq.gz SRR3182445_exome_sequencing_of_case2_techrep_2.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/005/SRR3182445/SRR3182445_1.fastq.gz . && mv SRR3182445_1.fastq.gz SRR3182445_exome_sequencing_of_case2_techrep_2_1.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/005/SRR3182445/SRR3182445_2.fastq.gz . && mv SRR3182445_2.fastq.gz SRR3182445_exome_sequencing_of_case2_techrep_2_2.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/007/SRR3182447/SRR3182447.fastq.gz . && mv SRR3182447.fastq.gz SRR3182447_exome_sequencing_of_case5_biorep_C.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/007/SRR3182447/SRR3182447_1.fastq.gz . && mv SRR3182447_1.fastq.gz SRR3182447_exome_sequencing_of_case5_biorep_C_1.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/007/SRR3182447/SRR3182447_2.fastq.gz . && mv SRR3182447_2.fastq.gz SRR3182447_exome_sequencing_of_case5_biorep_C_2.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/006/SRR3182446/SRR3182446.fastq.gz . && mv SRR3182446.fastq.gz SRR3182446_exome_sequencing_of_case5_techrep_2.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/006/SRR3182446/SRR3182446_1.fastq.gz . && mv SRR3182446_1.fastq.gz SRR3182446_exome_sequencing_of_case5_techrep_2_1.fastq.gz
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR318/006/SRR3182446/SRR3182446_2.fastq.gz . && mv SRR3182446_2.fastq.gz SRR3182446_exome_sequencing_of_case5_techrep_2_2.fastq.gz

##########################数据下载完成后，根据规则，重命名。#######################

#####################使用sra-tools的prefetch下载sra数据#########
#SRR_Acc_List.txt里文件格式
#SRR3182446
#SRR3182447
#...
#SRR3182418

cat SRR_Acc_List.txt | while read id
do
	prefetch ${id} -O  ./0.sra
done


####SRA转fq格式
cd ./0.sra
## fasterq-dump
cat config | while read id
do
	time fasterq-dump -3 -e 16 ${id}.sra -O ../1.raw_fq --outfile ${id}.fastq  >> ../1.raw_fq/${id}_sra2fq.log 2>&1
	time pigz -p 10 -f ../1.raw_fq/${id}_1.fastq >../1.raw_fq/${id}_1.fastq.gz
	time pigz -p 10 -f ../1.raw_fq/${id}_2.fastq >../1.raw_fq/${id}_2.fastq.gz
done 



cd ~/wes_cancer/project/1.raw_fq
##删除大小<1M的文件
find ./*gz -size 1M | xargs rm
##重命名文件名(截取字符串31-100位的字符)
ls *.gz | while read id; do mv ${id} ${id:31:100}; done

##Fastqc查看数据质量
cat config | while read id
do
	fastqc --outdir ./3.qc/raw_qc/ --threads 16 ./1.raw_fq/${id}*.fastq.gz >> ./3.qc/raw_qc/${id}_fastqc.log 2>&1 
done 
#multiqc讲fastqc生成的多个报告文件整合成为一个报告文件
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
genome=/share/home/chaimao1/human/hg38/hg38.fa
cd /share/home/chaimao1/wes_cancer/data
bwa index -a bwtsw -p hg38 $genome 

cd ~/wes_cancer/project
## bwa.sh 开始比对
INDEX= /share/home/chaimao1/wes_cancer/data/hg38
##直接结合samtools，把sam文件转换成bam文件
cat config | while read id
do
	echo "start bwa for ${id}" `date`
	fq1=./2.clean_fq/${id}_1_val_1.fq.gz
	fq2=./2.clean_fq/${id}_2_val_2.fq.gz
	bwa mem -M -t 16 -R "@RG\tID:${id}\tSM:${id}\tLB:WXS\tPL:Illumina" ${INDEX} ${fq1} ${fq2} | samtools view --threads 24 -bS - > ./4.align/${id}.bam
	echo "end bwa for ${id}" `date`
done


#####提取指定的文件的指定染色体

##samtools index case5_germline.bam
##samtools view -h case5_germline.bam chr12 | samtools view -Sb - > small.bam
##samtools index small.bam

######可以进一步提取指定区间内的bam文件
samtools index ./4.align/case5_germline.bam
samtools view -h ./4.align/case5_germline.bam chr1:25151000-30151000 | samtools view -Sb - > ./4.align/chr1_small.bam
samtools index ./4.align/chr1_small.bam


##############比对结果的bam文件质控####################################
## stats.sh
####方法1：使用samtools 的stats来检测
cat config | while read id
do
	bam=./4.align/${id}.bam
	samtools stats -@ 16 --reference $genome ${bam} > ./4.align/stats/${id}.stat
	plot-bamstats -p ./4.align/stats/${id} ./4.align/stats/${id}.stat
done

####方法2：使用qualimap来检测
cat config | while read id
do
	qualimap bamqc --java-mem-size=10G -gff ~/wes_cancer/data/hg38.exon.bed -nr 100000 -nw 500 -nt 16 -bam ./4.align/${id}.bam -outdir ./qualimap/4.align/${id}
done


#plot-bamstats -p ./4.align/stats/case1_germline ./4.align/stats/case1_germline.stat



###########################GATK4寻找变异#######################################################
##数据准备
# dbsnp_146.hg38.vcf.gz 
# dbsnp_146.hg38.vcf.gz.tbi 
# Mills_and_1000G_gold_standard.indels.hg38.vcf.gz 
# Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi 
# Homo_sapiens_assembly38.fasta
# Homo_sapiens_assembly38.fasta.gz 
# Homo_sapiens_assembly38.fasta.fai 
# Homo_sapiens_assembly38.dict 
# 1000G_phase1.snps.high_confidence.hg38.vcf.gz 
# 1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi

####GATK去除重复序列片段
GATK=~/soft/gatk/gatk-4.1.7.0/gatk

cat config  | while read id
do
	BAM=./4.align/${id}.bam
	if [ ! -f ./5.gatk/ok.${id}_marked.status ]
	then
		echo "start MarkDuplicates for ${id}" `date`
		$GATK --java-options "-Xmx20G -Djava.io.tmpdir=./" MarkDuplicates \
		-I ${BAM} \
		--REMOVE_DUPLICATES=true \
		-O ./5.gatk/${id}_marked.bam \
		-M ./5.gatk/${id}.metrics \
		1>./5.gatk/${id}_log.mark 2>&1 
		
		if [ $? -eq 0 ]
		then
			touch ./5.gatk/ok.${id}_marked.status
		fi
		echo "end MarkDuplicates for ${id}" `date`
		samtools index -@ 16 -m 4G -b ./5.gatk/${id}_marked.bam ./5.gatk/${id}_marked.bai
	fi
done

##碱基质量重新校正BQSR
snp=~/wes_cancer/data/dbsnp_146.hg38.vcf.gz
indel=~/wes_cancer/data/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
ref=/share/home/chaimao1/human/hg38/hg38.fa
cat config  | while read id
do
	if [ ! -f ./5.gatk/${id}_bqsr.bam ]
	then
		echo "start BQSR for ${id}" `date`
		$GATK --java-options "-Xmx20G -Djava.io.tmpdir=./"  BaseRecalibrator \
		-R $ref  \
		-I ./5.gatk/${id}_marked.bam  \
		--known-sites ${snp} \
		--known-sites ${indel} \
		-O ./5.gatk/${id}_recal.table \
		1>./5.gatk/${id}_log.recal 2>&1 
		
		$GATK --java-options "-Xmx20G -Djava.io.tmpdir=./"  ApplyBQSR \
		-R $ref  \
		-I ./5.gatk/${id}_marked.bam  \
		-bqsr ./5.gatk/${id}_recal.table \
		-O ./5.gatk/${id}_bqsr.bam \
		1>./5.gatk/${id}_log.ApplyBQSR  2>&1 
		
		echo "end BQSR for ${id}" `date`
	fi
done

####找germline SNPs+INDELs  生殖细胞的单碱基突变
GATK=~/soft/gatk/gatk-4.1.7.0/gatk
snp=~/wes_cancer/data/dbsnp_146.hg38.vcf.gz
indel=~/wes_cancer/data/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
bed=~/wes_cancer/data/hg38.exon.bed

cat config  | while read id
do
	echo "start HC for ${id}" `date`
	$GATK --java-options "-Xmx20G -Djava.io.tmpdir=./" HaplotypeCaller -ERC GVCF \
	-R ${ref} \
	-I ./5.gatk/${id}_bqsr.bam \
	--dbsnp ${snp} \
	-L ${bed} \
	-O ./5.gatk/${id}_raw.vcf \
	1>./5.gatk/${id}_log.HC 2>&1
	echo "end HC for ${id}" `date`

done

cd ./5.gatk/gvcf
for chr in chr{1..22} chrX chrY chrM
do

time $GATK --java-options "-Xmx20G -Djava.io.tmpdir=./" GenomicsDBImport \
-R ${ref} \
$(ls ./*raw.vcf | awk '{print "-V "$0" "}') \
-L ${chr} \
--genomicsdb-workspace-path gvcfs_${chr}.db

time $GATK --java-options "-Xmx20G -Djava.io.tmpdir=./" GenotypeGVCFs \
-R ${ref} \
-V gendb://gvcfs_${chr}.db \
-O gvcfs_${chr}.vcf

done

$GATK --java-options "-Xmx20G -Djava.io.tmpdir=./" GatherVcfs \
$(for i in {1..22} X Y M;do echo "-I gvcfs_chr${i}.vcf" ;done) \
-O merge.vcf


## mutect.sh  使用Mutect2检测体细胞的单碱基突变
#config2的文件格式：
#case5_germline case5_techrep_2

cat config2 | while read id
do
	arr=(${id})
	sample=${arr[1]}
	T=./5.gatk/${arr[1]}_bqsr.bam
	N=./5.gatk/${arr[0]}_bqsr.bam
	echo "start Mutect2 for ${id}" `date`
	$GATK  --java-options "-Xmx20G -Djava.io.tmpdir=./"  Mutect2 -R ${ref} \
	-I ${T} -tumor  $(basename "$T" _bqsr.bam) \
	-I ${N} -normal $(basename "$N" _bqsr.bam) \
	-L ${bed}  \
	-O ./6.mutect/${sample}_mutect2.vcf

	$GATK  FilterMutectCalls \
  -R ${ref} \
	-V ./6.mutect/${sample}_mutect2.vcf \
	-O ./6.mutect/${sample}_somatic.vcf
	echo "end Mutect2 for ${id}" `date`

	cat ./6.mutect/${sample}_somatic.vcf | perl -alne '{if(/^#/){print}else{next unless $F[6] eq "PASS";next if $F[0] =~/_/;print } }' > ./6.mutect/${sample}_filter.vcf
done


















##下载安装annovar
cd ~/wes_cancer/biosoft
# wget 下载地址
tar -zxvf annovar.latest.tar.gz
cd annovar

##在服务器端下载此文件失败，可以直接在本地下载后再上传，windows下载速度比较快。
#nohup ./annotate_variation.pl -downdb -webfrom annovar gnomad_genome --buildver hg38 humandb/ >down.log 2>&1 &
#本地下载文件地址
#http://www.openbioinformatics.org/annovar/download/hg38_gnomad_genome.txt.gz
#http://www.openbioinformatics.org/annovar/download/hg38_gnomad_genome.txt.idx.gz
































