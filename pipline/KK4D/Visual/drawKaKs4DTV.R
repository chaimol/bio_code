setwd("e:/bioinformation_center/pipline/KK4D/Visual")
### 载入R包
library(ggplot2)

#读取wgd-ksd产出的文件
data_kaks <- read.csv("Ath_Spo.kaks4DTv.csv",header = TRUE,na.strings = "")
#删除包含NA的行
#data_kaks <- na.omit(data_kaks)
##Ks
p0 <- ggplot(data=data_kaks, aes(x=Ks,fill=species)) +
  geom_density(alpha=0.4,color="red")+ xlab('Ks value')+ylab('Density')+
  xlim(0,15)+theme_classic() #设置x的最大值
p0
ggsave('Ks.pdf',dpi=300)


##Ka
p1 <- ggplot(data=data_kaks, aes(x=Ka)) + 
  geom_density(alpha=0.4,color="green")+ xlab('Ka value')+ylab('Density')+xlim(0,2)+theme_classic()
#labs(title = "Distribution of Ka distances", size = 1.5)+guides()
p1
ggsave('Ka.pdf',dpi=300)

##Ka/Ks
data_kaks <- na.omit(data_kaks)
data_kaks$Ka.Ks <- data_kaks$Ka/data_kaks$Ks
p2 <- ggplot(data=data_kaks, aes(x=Ka.Ks)) + 
  geom_density(alpha=0.4,color="orange")+ xlab('Ka/Ks value')+ylab('Density')+xlim(0,2)+theme_classic()
p2
ggsave('KaKs.pdf',dpi=300)

p4 <- ggplot(data=data_kaks, aes(x=dtv)) +
  geom_density(alpha=0.4,color="pink")+ xlab('4Dtv value')+ylab('Density')+
  theme_classic() #设置x的最大值
p0
ggsave('Ks.pdf',dpi=300)




install.packages("ggalluvial")
remotes::install_github("corybrunson/ggalluvial", ref = "optimization")
library("ggplot2")
library("ggalluvial")
titanic_wide <- data.frame(Titanic)
head(titanic_wide)
#>   Class    Sex   Age Survived Freq
#> 1   1st   Male Child       No    0
#> 2   2nd   Male Child       No    0
#> 3   3rd   Male Child       No   35
#> 4  Crew   Male Child       No    0
#> 5   1st Female Child       No    0
#> 6   2nd Female Child       No    0
ggplot(as.data.frame(Titanic),
       aes(y = Freq,
           axis1 = Survived, axis2 = Sex, axis3 = Class)) +
  geom_alluvium(aes(fill = Class),
                width = 0, knot.pos = 0, reverse = FALSE) +
  guides(fill = FALSE) +
  geom_stratum(width = 1/8, reverse = FALSE) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            reverse = FALSE) +
  scale_x_continuous(breaks = 1:3, labels = c("Survived", "Sex", "Class")) +
  coord_flip() +
  theme_classic()+theme(axis.line = element_blank())+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+ #NULL
  theme(axis.ticks.y=element_blank(),axis.text.y = element_text(face="bold.italic"))

ggplot(as.data.frame(Titanic),
       aes(y = Freq,
           axis1 = Survived, axis2 = Sex)) +
  guides(fill = FALSE) +
  geom_stratum(width = 1/8, reverse = FALSE) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            reverse = FALSE) +
  scale_x_continuous(breaks = 1:2, labels = c("Survived", "Sex")) +
  coord_flip() +
  theme_classic()+theme(axis.line = element_blank())+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+ #NULL
  theme(axis.ticks.y=element_blank(),axis.text.y = element_text(face="bold.italic"))


ggplot(as.data.frame(UCBAdmissions),
       aes(y = Freq, axis1 = Gender, axis2 = Dept)) +
  geom_alluvium(aes(fill = Admit), width = 1/12) +
  geom_stratum(width = 1/12, fill = "#b5f6ca") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Gender", "Dept"), expand = c(.05, .05)) +
  scale_fill_brewer(type = "qual", palette = "Set1")




