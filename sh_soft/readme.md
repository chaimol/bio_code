# 使用SimpleM计算最佳P的阈值
### 如果你有格式为0,1,2的数字矩阵，格式和snpSample.txt一致，可以直接使用`Rscript SimpleM.R snpSample.txt`,输出的第3行是P值
### 如果你只有原始的vcf文件，则需要运行脚本`bash getPvalue.bash vcffile abbr`.自动完成变异的缺失值填充，输出在abbr.SimpleM.pvalue.第3行是P值.
- plink 需要手动给予执行权限 `chmod 757 plink`
- SimpleM.R 是一个分析变异的矫正的P值的脚本，由我修改后实现的
- snpSample.txt 是一个示例数字矩阵文件，要求所有的变异必须是数字，而且不能有缺失值。
