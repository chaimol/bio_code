### 计算KaKs，4DTv或进化树等中间步骤的小脚本。此目录所有程序都来源于其他人的代码，非本人编写。
- axt2one-line.py #把多行的fasta转为1行
- calc_4DTv_eff_vcf.py  #从snpEff的输出结果vcf文件，计算4DTv，分析进化关系
- calculate_4DTV_correction.pl #计算4DTv
- genome_statistic.pl #统计基因组的所有信息（长度，contig10-90等）
- removeRedundantProteins.py #从每个基因组的蛋白文件中找出最长的转录本的蛋白序列。产出文件用于比对进化关系。
- vcf2phylip.py #转换vcf格式的snp文件为PHYLIP, FASTA, NEXUS, or binary NEXUS 用于构建进化树