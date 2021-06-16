if(!require(dplyr))install.packages("scholar")
if(!require(dplyr))install.packages("dplyr")
library(scholar)
library(dplyr)

#scholar这个包实质是爬取了谷歌学术，整理对应的内容。简单粗暴好用。

#设定期刊的名称的向量
# jn = c("bioinformatics","methods is ecology and evolution","molecular biosystems","molecular biology and evolution")
# 
# #获取对应期刊的影响因子
# get_impactfactor(jn)
# 
# #获取cns的最新影响因子
# cns=c("cell","nature","science")
# get_impactfactor(cns)
# 
# #获取有bioinformatics相关的所有期刊
# bioinfor <- agrep("bioinformatics",scholar:::impactfactor$Journal,ignore.case = T,value = T,max.distance = 0.05)
# get_impactfactor(bioinfor)


journal_name <- read.csv("JournalName.csv",header = TRUE)$Journal
output <- get_impactfactor(journal_name)
write.csv(output,"journal_IF.csv")