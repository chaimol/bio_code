#用于学习使用单细胞RNA-seq测序分析技术

#数据下载地址 
# https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz




setwd("E:/bioinformation_center/model/scRNA-seq/filtered_gene_bc_matrices/")
# install.packages('BiocManager') # 没有的话就安装一个
# BiocManager::install('multtest')
# install.packages('Seurat')

install.packages("patchwork")
library(dplyr)
library(Seurat)
library(patchwork)

#10x的数据读取
# Load the PBMC dataset 读取数据
pbmc.data <- Read10X(data.dir ="./hg19/")

dim(pbmc.data)



