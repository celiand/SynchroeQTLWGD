genepair<-read.table(file="genetabgenepair02.txt",header=TRUE)
length(genepair[genepair$regulation=="singreg","regulation"])


##get expression of gene 1 and gene 2 for each pair
library(tidyverse)

Expression <- read_table2("/mnt/project/MSLab/Marie/2022_Synchrosmolt/QTLTools/Synchro_mt_TPM.bed")    %>% select(-c(1,2,3,5,6))
Genotype <- read_table2("/mnt/project/MSLab/Marie/2022_Synchrosmolt/QTLTools/imputed_0.01_2730_smolts.traw")   %>% select(-c(1,3,4,5,6))
Covariate <-  read_table2("covariate_synchro.txt") 
## fix the name in genotype
names(Genotype)[2:ncol(Genotype)] <- str_remove(names(Genotype)[2:ncol(Genotype)], "0_")
# scale the expression
Expression[, 2:ncol(Expression)] <- t(apply(Expression[, 2:ncol(Expression)], 1, scale))


##extract only LL individuals



genepair<-read.table(file="genetabgenepair01.txt",header=TRUE)



genepair1<-Expression[Expression$pid%in%genepair$gene1,]
genepair2<-Expression[Expression$pid%in%genepair$gene2,]



# Pivot genotype table to long format
#geno_table_long_p1 <- genepair1 %>%
#pivot_longer(cols = -SNP, names_to = "id", values_to = "genotype")

#Pivot expression table to long format
expr_table_long1 <- genepair1 %>%
  pivot_longer(cols = -pid, names_to = "id", values_to = "expression")

expr_table_long2 <- genepair2 %>%
  pivot_longer(cols = -pid, names_to = "id", values_to = "expression")

covar_kept<-Covariate[c(3,5),]
id<-colnames(covar_kept)
tcovar<-t(covar_kept)
long_covar<-cbind.data.frame(id,tcovar)
long_covar<-long_covar[-1,]
colnames(long_covar)<-c("id","sex","LL")

long_covar_LL<-long_covar[,c(1,3)]

#longcovarll<-long_covar_LL[long_covar_LL$LL=="LL",] ## 907 individuals



mergedpair1 <- merge(expr_table_long1 , long_covar, by = "id")
mergedpair2 <- merge(expr_table_long2 , long_covar, by = "id")

##filter

mergedpair1_filter<- mergedpair1[mergedpair1$LL=="LL",]
mergedpair2_filter<- mergedpair2[mergedpair2$LL=="LL",]

##can keep more column if wanted
mergedpair1_filter<-mergedpair1_filter[,1:3]
mergedpair2_filter<-mergedpair2_filter[,1:3]


library(dplyr)

# Left join table 1 with table 2 based on gene1 and pid
joined_table1_2 <- left_join(genepair, mergedpair1_filter, by = c("gene1" = "pid"))
joined_table1_2<-joined_table1_2[!is.na(joined_table1_2$expression),]

# Left join result with table 3 based on gene2 and pid
joined_table1_2_3 <- left_join(joined_table1_2, mergedpair2_filter, by = c("id"="id","gene2" = "pid"))
joined_table1_2_3 <- joined_table1_2_3 [!is.na(joined_table1_2_3$expression.y),]

colnames(joined_table1_2_3)[9:10]<-c("Expression_gene1","Expression_gene2")
write.table(joined_table1_2_3,file="expression_gene_pair02.txt",row.names=FALSE)




##analysis
library(ggplot2)

gene_pair_tab<-read.table(file="expression_gene_pair02.txt",header=TRUE)

ggplot(gene_pair_tab[gene_pair_tab$ID=="ENSSSAG00000034990_Ssal_ENSSSAG00000071166_Ssal",],aes(x=Expression_gene1,y=Expression_gene2))+geom_point()+theme_bw(20)

##coefficient

# calculate spearman coefficient for each condition
result <- gene_pair_tab %>% 
  group_by(ID) %>% 
  summarize(spearman_coeff = cor(Expression_gene1, Expression_gene2, method = "spearman"))

##if want p_val:
data<-c()
for (i in unique(gene_pair_tab$ID)){
  subtab<-gene_pair_tab[gene_pair_tab$ID==i,]
  res<-cor.test(subtab$Expression_gene1, subtab$Expression_gene2, method = "spearman")
  estimate<-res$estimate
  pval<-res$p.value
  data<-rbind.data.frame(data,c(i,estimate,pval))
}
colnames(data)<-c("ID","coeff","pval")

genepair<-read.table(file="genetabgenepair02.txt",header=TRUE)

results2<-left_join(result, genepair, by = c("ID"="ID"))


##stats calculation
results3<-as.data.frame(results2)

ks.test(results3[results3$regulation=="bothsamereg","spearman_coeff"],results3[results3$regulation=="bothdiffreg","spearman_coeff"])
ks.test(results3[results3$regulation=="bothsamereg","spearman_coeff"],results3[results3$regulation=="singreg","spearman_coeff"])
ks.test(results3[results3$regulation=="bothsamereg","spearman_coeff"],results3[results3$regulation=="noreg","spearman_coeff"])
ks.test(results3[results3$regulation=="bothdiffreg","spearman_coeff"],results3[results3$regulation=="singreg","spearman_coeff"])
ks.test(results3[results3$regulation=="bothdiffreg","spearman_coeff"],results3[results3$regulation=="noreg","spearman_coeff"])
ks.test(results3[results3$regulation=="singreg","spearman_coeff"],results3[results3$regulation=="noreg","spearman_coeff"])


x<-abs(results3[results3$regulation=="bothsamereg","spearman_coeff"])
#x<-abs(results3[results3$regulation=="bothdiffreg","spearman_coeff"])
#x<-abs(results3[results3$regulation=="singreg","spearman_coeff"])
#x<-abs(results3[results3$regulation=="noreg","spearman_coeff"])

mean(x)
sd(x)/sqrt(length(x))

##plots

ggplot(results2,aes(x=spearman_coeff,color=regulation))+geom_density()+theme_bw(20)


ggplot(results2,aes(x=spearman_coeff,fill=regulation))+geom_histogram(color="black")+theme_bw(20)+facet_wrap(~regulation,scales = 'free_y')+xlab("Sperman coefficient")


results2[results2$regulation=="noreg","regulation"]<-"B: No eQTL"
results2[results2$regulation=="singreg","regulation"]<-"D: Single copy eQTL"
results2[results2$regulation=="bothdiffreg","regulation"]<-"A: Distinct eQTL"
results2[results2$regulation=="bothsamereg","regulation"]<-"C: Shared eQTL"

ggplot(results2,aes(x=spearman_coeff,fill=regulation))+geom_histogram(color="black")+theme_bw(20)+facet_wrap(~regulation,scales = 'free_y')+xlab("Spearman coefficient")+scale_fill_manual(values=c("#F9B27C","#99B8EA","#8CE86D","#D186D8"))

df2 <- results2 %>%
  group_by(regulation) %>%
  summarise(Mean = mean(spearman_coeff))

ggplot(results2,aes(x=spearman_coeff,fill=regulation))+geom_histogram(color="black")+theme_bw(20)+facet_wrap(~regulation,scales = 'free_y')+xlab("Spearman coefficient")+scale_fill_manual(values=c("#F9B27C","#99B8EA","#8CE86D","#D186D8"))+geom_vline(data = df2, mapping = aes(xintercept = Mean),color="firebrick3",size=1.5)


##median calculation
results2<-as.data.frame(results2)
sharedeqtlcoeff<-results2[results2$regulation=="bothsamereg","spearman_coeff"]

medianshared<-median(sharedeqtlcoeff) ##0.816

length(results2[results2$regulation=="bothdiffreg","spearman_coeff"])
length(results2[results2$regulation=="bothdiffreg" & results2$spearman_coeff>=medianshared,"spearman_coeff"]) ##152/1834 --> 8%

length(results2[results2$regulation=="singreg","spearman_coeff"])
length(results2[results2$regulation=="singreg" & results2$spearman_coeff>=medianshared,"spearman_coeff"]) ## 108/1661 --> 6%


##

results2<-left_join(data, genepair, by = c("ID"="ID"))

## add A with high LD to C
## 2 - 5 - 7
results2[results2$ID%in%pairtokeep,"regulation"]<-"bothsamereg"

results2[results2$regulation=="noreg","regulation"]<-"B: No eQTL"
results2[results2$regulation=="singreg","regulation"]<-"D: Single copy eQTL"
results2[results2$regulation=="bothdiffreg","regulation"]<-"A: Distinct eQTL"
results2[results2$regulation=="bothsamereg","regulation"]<-"C: Shared eQTL"

ggplot(results2,aes(x=spearman_coeff,fill=regulation))+geom_histogram(color="black")+theme_bw(20)+facet_wrap(~regulation,scales = 'free_y')+xlab("Spearman coefficient")+scale_fill_manual(values=c("#F9B27C","#99B8EA","#8CE86D","#D186D8"))


##investigate gene with negative correlation

#negativedata<-results2[results2$spearman_coeff<0.25 & results2$regulation=="bothsamereg",]
negativedata<-results2[results2$spearman_coeff<0 & results2$regulation=="bothsamereg",]
negativedata<-results2[results2$spearman_coeff<0 & results2$regulation=="C: Shared eQTL",]

#negativedata<-results2[results2$spearman_coeff>0 & results2$regulation=="bothsamereg",]
#negativedata<-negativedata[negativedata$ID=="ENSSSAG00000034990_Ssal_ENSSSAG00000071166_Ssal",]

#ggplot(gene_pair_tab[gene_pair_tab$ID=="ENSSSAG00000011938_Ssal_ENSSSAG00000108018_Ssal",],aes(x=Expression_gene1,y=Expression_gene2))+geom_point()+theme_bw(20)+ggtitle("correlation for pair ENSSSAG00000011938_Ssal_ENSSSAG00000108018_Ssal")

#ggplot(gene_pair_tab[gene_pair_tab$ID=="ENSSSAG00000086178_Ssal_ENSSSAG00000087231_Ssal",],aes(x=Expression_gene1,y=Expression_gene2))+geom_point()+theme_bw(20)+ggtitle("correlation for pair ENSSSAG00000086178_Ssal_ENSSSAG00000087231_Ssal")

#ggplot(gene_pair_tab[gene_pair_tab$ID=="ENSSSAG00000113858_Ssal_ENSSSAG00000120807_Ssal",],aes(x=Expression_gene1,y=Expression_gene2))+geom_point()+theme_bw(20)+ggtitle("correlation for pair ENSSSAG00000113858_Ssal_ENSSSAG00000120807_Ssal")

negativedata

ggplot(gene_pair_tab[gene_pair_tab$ID=="ENSSSAG00000011938_Ssal_ENSSSAG00000108018_Ssal",],aes(x=Expression_gene1,y=Expression_gene2))+geom_point()+theme_bw(20)+ggtitle("correlation for pair ENSSSAG00000011938_Ssal_ENSSSAG00000108018_Ssal")

ggplot(gene_pair_tab[gene_pair_tab$ID=="ENSSSAG00000086178_Ssal_ENSSSAG00000087231_Ssal",],aes(x=Expression_gene1,y=Expression_gene2))+geom_point()+theme_bw(20)+ggtitle("correlation for pair ENSSSAG00000086178_Ssal_ENSSSAG00000087231_Ssal")

ggplot(gene_pair_tab[gene_pair_tab$ID=="ENSSSAG00000113858_Ssal_ENSSSAG00000120807_Ssal",],aes(x=Expression_gene1,y=Expression_gene2))+geom_point()+theme_bw(20)+ggtitle("correlation for pair ENSSSAG00000113858_Ssal_ENSSSAG00000120807_Ssal")


## what about others ?

ggplot(gene_pair_tab[gene_pair_tab$ID=="ENSSSAG00000009601_Ssal_ENSSSAG00000073923_Ssal",],aes(x=Expression_gene1,y=Expression_gene2))+geom_point()+theme_bw(20)+ggtitle("correlation for pair ENSSSAG00000009601_Ssal_ENSSSAG00000073923_Ssal")

ggplot(gene_pair_tab[gene_pair_tab$ID=="ENSSSAG00000053868_Ssal_ENSSSAG00000057392_Ssal",],aes(x=Expression_gene1,y=Expression_gene2))+geom_point()+theme_bw(20)+ggtitle("ENSSSAG00000053868_Ssal_ENSSSAG00000057392_Ssal")



negativedata<-results2[results2$spearman_coeff<=-0.1,]

library(googledrive)

for( i in negativedata$ID){
  pairstatus<-negativedata[negativedata$ID==i,"regulation"]
  title<-paste("correlation for pair ",i,sep="")
  p<-ggplot(gene_pair_tab[gene_pair_tab$ID==i,],aes(x=Expression_gene1,y=Expression_gene2))+geom_point()+theme_bw(20)+ggtitle(title)
  fname<-paste(pairstatus,"plotgp",i,".png",sep="")
  ggsave(filename=fname,plot=p,device="png",path="/net/fs-2/scale/OrionStore/Home/cedi/phDSalmon/rna_analysis/plotgenepair/",width=15,height=12,units="in")
  #pathlocal<-paste("/net/fs-2/scale/OrionStore/Home/cedi/phDSalmon/rna_analysis/plotgenepair/",fname,sep="")
  #drive_upload(media=pathlocal,name=fname,path="genepair",overwrite=TRUE)
}

###debug drive
fname<-"testpair"
p<-ggplot(gene_pair_tab[gene_pair_tab$ID=="ENSSSAG00000113858_Ssal_ENSSSAG00000120807_Ssal",],aes(x=Expression_gene1,y=Expression_gene2))+geom_point()+theme_bw(20)+ggtitle("correlation for pair ENSSSAG00000113858_Ssal_ENSSSAG00000120807_Ssal")
ggsave(filename=fname,plot=p,device="png",path="/net/fs-2/scale/OrionStore/Home/cedi/phDSalmon/rna_analysis/plotgenepair/",width=30,height=24,units="in")
pathlocal<-paste("/net/fs-2/scale/OrionStore/Home/cedi/phDSalmon/rna_analysis/plotgenepair/",fname,sep="")
drive_upload(media=pathlocal,name=fname,path="genepair",overwrite=TRUE)

library(googledrive)
fname<-"testpair"
pathlocal<-paste("/net/fs-2/scale/OrionStore/Home/cedi/phDSalmon/rna_analysis/plotgenepair/",fname,sep="")
drive_upload(media=pathlocal,name=fname,path="genepair",overwrite=TRUE)

##get their raw expression
library(tidyverse)
Expression <- read_table2("/mnt/project/MSLab/Marie/2022_Synchrosmolt/QTLTools/Synchro_mt_TPM.bed")    %>% select(-c(1,2,3,5,6))

neggeneexpp1<- Expression[Expression$pid%in%negativedata$gene1,]
neggeneexpp2<- Expression[Expression$pid%in%negativedata$gene2,]

neggene_long1<- neggeneexpp1 %>%
  pivot_longer(cols = -pid, names_to = "id", values_to = "expression")

neggene_long2<- neggeneexpp2 %>%
  pivot_longer(cols = -pid, names_to = "id", values_to = "expression")



covar_kept<-Covariate[c(3,5),]
id<-colnames(covar_kept)
tcovar<-t(covar_kept)
long_covar<-cbind.data.frame(id,tcovar)
long_covar<-long_covar[-1,]
colnames(long_covar)<-c("id","sex","LL")

long_covar_LL<-long_covar[,c(1,3)]

#longcovarll<-long_covar_LL[long_covar_LL$LL=="LL",] ## 907 individuals



mergedpair1neg <- merge(neggene_long1 , long_covar, by = "id")
mergedpair2neg <- merge(neggene_long2 , long_covar, by = "id")

##filter

mergedpair1neg_filter<- mergedpair1neg[mergedpair1neg$LL=="LL",]
mergedpair2neg_filter<- mergedpair2neg[mergedpair2neg$LL=="LL",]

##can keep more column if wanted
mergedpair1neg_filter<-mergedpair1neg_filter[,1:3]
mergedpair2neg_filter<-mergedpair2neg_filter[,1:3]


library(dplyr)

# Left join table 1 with table 2 based on gene1 and pid
joined_table1_2neg <- left_join(genepair, mergedpair1neg_filter, by = c("gene1" = "pid"))
joined_table1_2neg<-joined_table1_2neg[!is.na(joined_table1_2neg$expression),]

# Left join result with table 3 based on gene2 and pid
joined_table1_2_3neg <- left_join(joined_table1_2neg, mergedpair2neg_filter, by = c("id"="id","gene2" = "pid"))
joined_table1_2_3neg <- joined_table1_2_3neg [!is.na(joined_table1_2_3neg$expression.y),]

colnames(joined_table1_2_3neg)[6:7]<-c("Expression_gene1","Expression_gene2")
write.table(joined_table1_2_3neg,file="expression_gene_pair_pos02.txt",row.names=FALSE)


##get their genotype

rna_full_table<-read.table(file="Synchro_uniq_fulltable_withinfo_04.txt",header=TRUE)
rna_full_table<-rna_full_table[rna_full_table$condition=="LL",]
neggene<-rna_full_table[rna_full_table$gene%in%negativedata$gene1 | rna_full_table$gene%in%negativedata$gene2,]
neggene$snpchr2<-ifelse(neggene$snpschrom<10,paste("ssa0",neggene$snpschrom,sep=""),paste("ssa",neggene$snpschrom,sep=""))
neggene$snpid<-paste(neggene$snpchr2,neggene$snpstart,sep="_")

snpid<-unique(neggene$snpid)
snpid<-snpid[-1]

gtkept<-Genotype[Genotype$SNP%in%snpid,]

geno_table_neg <- gtkept %>% pivot_longer(cols = -SNP, names_to = "id", values_to = "genotype")

long_covar$id<-paste("0_",long_covar$id,sep="")
gene_tab_neg_merge <- merge(geno_table_neg  , long_covar, by = "id")  
gene_tab_neg_merge_filter<- gene_tab_neg_merge[gene_tab_neg_merge$LL=="LL",]

pairdata<-data.frame(ID=c("ENSSSAG00000011938_Ssal_ENSSSAG00000108018_Ssal","ENSSSAG00000086178_Ssal_ENSSSAG00000087231_Ssal","ENSSSAG00000113858_Ssal_ENSSSAG00000120807_Ssal"),SNP=c("ssa05_72178688","ssa07_56580022","ssa26_27976437"))

pairdata<-data.frame(ID=c("ENSSSAG00000011938_Ssal_ENSSSAG00000108018_Ssal","ENSSSAG00000086178_Ssal_ENSSSAG00000087231_Ssal","ENSSSAG00000113858_Ssal_ENSSSAG00000120807_Ssal","ENSSSAG00000043068_Ssal_ENSSSAG00000054986_Ssal","ENSSSAG00000043068_Ssal_ENSSSAG00000054986_Ssal","ENSSSAG00000001778_Ssal_ENSSSAG00000045167_Ssal","ENSSSAG00000001778_Ssal_ENSSSAG00000045167_Ssal","ENSSSAG00000047981_Ssal_ENSSSAG00000057888_Ssal","ENSSSAG00000047981_Ssal_ENSSSAG00000057888_Ssal","ENSSSAG00000092357_Ssal_ENSSSAG00000109189_Ssal","ENSSSAG00000092357_Ssal_ENSSSAG00000109189_Ssal"),SNP=c("ssa05_72178688","ssa07_56580022","ssa26_27976437","ssa14_50870510","ssa14_51354372","ssa14_51370830","ssa14_50870510","ssa14_80446627","ssa14_50816381","ssa17_62355402","ssa17_62354384"))

pairdata<-data.frame(ID=c("ENSSSAG00000034990_Ssal_ENSSSAG00000071166_Ssal"),SNP=c("ssa14_50871630"))


mergegt<-merge(gene_tab_neg_merge_filter,pairdata,by="SNP")

write.table(mergegt,file="genotype_gene_pair_pos02.txt",row.names=FALSE)

##plot
genetab<-read.table(file = "expression_gene_pair_neg01.txt",header=TRUE)
genotab<-read.table(file = "genotype_gene_pair_neg01.txt",header=TRUE)

fulltabneg<-merge(genetab,genotab,by=c("ID","id"))
fulltabneg$sumexpression<-fulltabneg$Expression_gene1+fulltabneg$Expression_gene2

## see if sum of expression is the same
data0<-fulltabneg[fulltabneg$genotype=="0" & fulltabneg$ID=="ENSSSAG00000011938_Ssal_ENSSSAG00000108018_Ssal","sumexpression"]
data1<-fulltabneg[fulltabneg$genotype=="1" & fulltabneg$ID=="ENSSSAG00000011938_Ssal_ENSSSAG00000108018_Ssal","sumexpression"]
data2<-fulltabneg[fulltabneg$genotype=="2" & fulltabneg$ID=="ENSSSAG00000011938_Ssal_ENSSSAG00000108018_Ssal","sumexpression"]

ks.test(data0,data1)
ks.test(data2,data1)
ks.test(data2,data0) ##no diff for this pair

data0<-fulltabneg[fulltabneg$genotype=="0" & fulltabneg$ID=="ENSSSAG00000086178_Ssal_ENSSSAG00000087231_Ssal","sumexpression"]
data1<-fulltabneg[fulltabneg$genotype=="1" & fulltabneg$ID=="ENSSSAG00000086178_Ssal_ENSSSAG00000087231_Ssal","sumexpression"]
data2<-fulltabneg[fulltabneg$genotype=="2" & fulltabneg$ID=="ENSSSAG00000086178_Ssal_ENSSSAG00000087231_Ssal","sumexpression"]

ks.test(data0,data1)
ks.test(data2,data1)
ks.test(data2,data0) ##diff for this pair

data0<-fulltabneg[fulltabneg$genotype=="0" & fulltabneg$ID=="ENSSSAG00000113858_Ssal_ENSSSAG00000120807_Ssal","sumexpression"]
data1<-fulltabneg[fulltabneg$genotype=="1" & fulltabneg$ID=="ENSSSAG00000113858_Ssal_ENSSSAG00000120807_Ssal","sumexpression"]
data2<-fulltabneg[fulltabneg$genotype=="2" & fulltabneg$ID=="ENSSSAG00000113858_Ssal_ENSSSAG00000120807_Ssal","sumexpression"]

ks.test(data0,data1)
ks.test(data2,data1)
ks.test(data2,data0) ##no diff for this pair



##


ggplot(fulltabneg[fulltabneg$ID=="ENSSSAG00000011938_Ssal_ENSSSAG00000108018_Ssal",],aes(x=Expression_gene1,y=Expression_gene2,color=as.factor(genotype)))+geom_point()+theme_bw(20)+ggtitle("correlation for pair ENSSSAG00000011938_Ssal_ENSSSAG00000108018_Ssal")

ggplot(fulltabneg[fulltabneg$ID=="ENSSSAG00000086178_Ssal_ENSSSAG00000087231_Ssal",],aes(x=Expression_gene1,y=Expression_gene2,color=as.factor(genotype)))+geom_point()+theme_bw(20)+ggtitle("correlation for pair ENSSSAG00000086178_Ssal_ENSSSAG00000087231_Ssal")

ggplot(fulltabneg[fulltabneg$ID=="ENSSSAG00000113858_Ssal_ENSSSAG00000120807_Ssal",],aes(x=Expression_gene1,y=Expression_gene2,color=as.factor(genotype)))+geom_point()+theme_bw(20)+ggtitle("correlation for pair ENSSSAG00000113858_Ssal_ENSSSAG00000120807_Ssal")


toplot<-ggplot(fulltabneg[fulltabneg$ID=="ENSSSAG00000011938_Ssal_ENSSSAG00000108018_Ssal",],aes(x=Expression_gene1,fill=as.factor(genotype)))+geom_density(alpha=0.5)+theme_bw(20)+theme(legend.position = "none")+ggtitle("correlation for collec10")+xlab("Expression copy 1")
rightplot<-ggplot(fulltabneg[fulltabneg$ID=="ENSSSAG00000011938_Ssal_ENSSSAG00000108018_Ssal",],aes(x=Expression_gene2,fill=as.factor(genotype)))+geom_density(alpha=0.5)+theme_bw(20)+coord_flip()+xlab("Expression copy 2")
scatter<-ggplot(fulltabneg[fulltabneg$ID=="ENSSSAG00000011938_Ssal_ENSSSAG00000108018_Ssal",],aes(x=Expression_gene1,y=Expression_gene2,color=as.factor(genotype)))+geom_point()+theme_bw(20)+theme(legend.position = "none")+xlab("Expression copy 1")+ylab("Expression copy 2")
empty <- ggplot() + 
  geom_point(aes(1,1), colour="white") +
  theme(                              
    plot.background = element_blank(), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.border = element_blank(), 
    panel.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  )

library(egg) 

ggarrange(toplot, empty, scatter, rightplot, 
          ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))

toplot<-ggplot(fulltabneg[fulltabneg$ID=="ENSSSAG00000086178_Ssal_ENSSSAG00000087231_Ssal",],aes(x=Expression_gene1,fill=as.factor(genotype)))+geom_density(alpha=0.5)+theme_bw(20)+theme(legend.position = "none")+ggtitle("correlation for tetraspanin-8")+xlab("Expression copy 1")
rightplot<-ggplot(fulltabneg[fulltabneg$ID=="ENSSSAG00000086178_Ssal_ENSSSAG00000087231_Ssal",],aes(x=Expression_gene2,fill=as.factor(genotype)))+geom_density(alpha=0.5)+theme_bw(20)+coord_flip()+xlab("Expression copy 2")
scatter<-ggplot(fulltabneg[fulltabneg$ID=="ENSSSAG00000086178_Ssal_ENSSSAG00000087231_Ssal",],aes(x=Expression_gene1,y=Expression_gene2,color=as.factor(genotype)))+geom_point()+theme_bw(20)+theme(legend.position = "none")+xlab("Expression copy 1")+ylab("Expression copy 2")

ggarrange(toplot, empty, scatter, rightplot, 
          ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))

toplot<-ggplot(fulltabneg[fulltabneg$ID=="ENSSSAG00000113858_Ssal_ENSSSAG00000120807_Ssal",],aes(x=Expression_gene1,fill=as.factor(genotype)))+geom_density(alpha=0.5)+theme_bw(20)+theme(legend.position = "none")+ggtitle("correlation for eif3")+xlab("Expression copy 1")
rightplot<-ggplot(fulltabneg[fulltabneg$ID=="ENSSSAG00000113858_Ssal_ENSSSAG00000120807_Ssal",],aes(x=Expression_gene2,fill=as.factor(genotype)))+geom_density(alpha=0.5)+theme_bw(20)+coord_flip()+xlab("Expression copy 2")+ylim(0,0.03)
scatter<-ggplot(fulltabneg[fulltabneg$ID=="ENSSSAG00000113858_Ssal_ENSSSAG00000120807_Ssal",],aes(x=Expression_gene1,y=Expression_gene2,color=as.factor(genotype)))+geom_point()+theme_bw(20)+theme(legend.position = "none")+xlab("Expression copy 1")+ylab("Expression copy 2")

ggarrange(toplot, empty, scatter, rightplot, 
          ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))

## make histogram
gene1<-fulltabneg[,c(1,2,3,6,8,9)]
gene2<-fulltabneg[,c(1,2,4,7,8,9)]

sumgene<-data.frame(ID=fulltabneg$ID,id=fulltabneg$id,gene=fulltabneg$ID,expression=fulltabneg$sumexpression,SNP=fulltabneg$SNP,genotype=fulltabneg$genotype)
colnames(sumgene)<-c("ID","id","gene","expression","SNP","genotype")
sumgene$pair<-"Sum of copies"

colnames(gene2)<-c("ID","id","gene","expression","SNP","genotype")
colnames(gene1)<-c("ID","id","gene","expression","SNP","genotype")

gene1$pair<-"copy 1"
gene2$pair<-"copy 2"


longgene<-rbind.data.frame(gene1,gene2,stringsAsFactors = FALSE)
longgene$name<-paste("Pair: ",sapply(strsplit(longgene$ID,"_Ssal"),FUN = `[[`, 1),sapply(strsplit(longgene$ID,"_Ssal"),FUN = `[[`, 2),"\nSNP: ",longgene$SNP,sep="")


ggplot(longgene,aes(x=as.factor(genotype),y=expression,fill=pair))+geom_boxplot()+theme_bw(20)+facet_wrap(~name,scales = "free_y",labeller = labeller(name = label_wrap_gen(width = 45)))+xlab("eQTL genotype")


longgene<-rbind.data.frame(gene1,gene2,sumgene,stringsAsFactors = FALSE)
longgene$name<-paste("Pair: ",sapply(strsplit(longgene$ID,"_Ssal"),FUN = `[[`, 1),sapply(strsplit(longgene$ID,"_Ssal"),FUN = `[[`, 2),"\nSNP: ",longgene$SNP,sep="")
longgene[longgene$name=="Pair: ENSSSAG00000011938_ENSSSAG00000108018\nSNP: ssa05_72178688","name"]<-"colec10"
longgene[longgene$name=="Pair: ENSSSAG00000086178_ENSSSAG00000087231\nSNP: ssa07_56580022","name"]<-"tetraspanin-8"
longgene[longgene$name=="Pair: ENSSSAG00000113858_ENSSSAG00000120807\nSNP: ssa26_27976437","name"]<-"eif3"

ggplot(longgene,aes(x=as.factor(genotype),y=expression,fill=pair))+geom_boxplot()+theme_bw(25)+facet_wrap(~name,scales = "free_y",labeller = labeller(name = label_wrap_gen(width = 45)))+xlab("eQTL genotype")+scale_fill_manual(values=c("firebrick1","dodgerblue1","darkorchid1"))

library(data.table)
dat <- data.table(longgene)
dat[name == "colec10",y_min := 0]
dat[name == "colec10",y_max := 70]
dat[name == "eif3",y_min := 0]
dat[name == "eif3",y_max := 500]
dat[name == "tetraspanin-8",y_min := 0]
dat[name == "tetraspanin-8",y_max := 2300]
ggplot(dat,aes(x=as.factor(genotype),y=expression,fill=pair))+geom_boxplot()+theme_bw(25)+facet_wrap(~name,scales = "free_y",labeller = labeller(name = label_wrap_gen(width = 45)))+xlab("eQTL genotype")+scale_fill_manual(values=c("firebrick1","dodgerblue1","darkorchid1"))+geom_blank(aes(y = y_min)) +geom_blank(aes(y = y_max))
ggplot(dat[dat$pair!="Sum of copies",],aes(x=as.factor(genotype),y=expression,fill=pair))+geom_boxplot()+theme_bw(25)+facet_wrap(~name,scales = "free_y",labeller = labeller(name = label_wrap_gen(width = 45)))+xlab("eQTL genotype")+scale_fill_manual(values=c("firebrick1","dodgerblue1","darkorchid1"))+geom_blank(aes(y = y_min)) +geom_blank(aes(y = y_max))


longgene[longgene$name=="Pair: ENSSSAG00000011938_ENSSSAG00000108018\nSNP: ssa05_72178688","name"]<-"colec10"
longgene[longgene$name=="Pair: ENSSSAG00000086178_ENSSSAG00000087231\nSNP: ssa07_56580022","name"]<-"tetraspanin-8"
longgene[longgene$name=="Pair: ENSSSAG00000113858_ENSSSAG00000120807\nSNP: ssa26_27976437","name"]<-"eif3"
#unique(longgene$name)

ggplot(longgene,aes(x=as.factor(genotype),y=expression,fill=pair))+geom_boxplot()+theme_bw(25)+facet_wrap(~name,scales = "free_y",labeller = labeller(name = label_wrap_gen(width = 45)))+xlab("eQTL genotype")


## for negative pairs with linkage
library(ggplot2)
genetab<-read.table(file = "expression_gene_pair_pos02.txt",header=TRUE)
genotab<-read.table(file = "genotype_gene_pair_pos02.txt",header=TRUE)

fulltabneg<-merge(genetab,genotab,by=c("ID","id"))
fulltabneg$sumexpression<-fulltabneg$Expression_gene1+fulltabneg$Expression_gene2



gene1<-fulltabneg[,c(1,2,3,6,8,9)]
gene2<-fulltabneg[,c(1,2,4,7,8,9)]

sumgene<-data.frame(ID=fulltabneg$ID,id=fulltabneg$id,gene=fulltabneg$ID,expression=fulltabneg$sumexpression,SNP=fulltabneg$SNP,genotype=fulltabneg$genotype)
colnames(sumgene)<-c("ID","id","gene","expression","SNP","genotype")
sumgene$pair<-"Sum of copies"

colnames(gene2)<-c("ID","id","gene","expression","SNP","genotype")
colnames(gene1)<-c("ID","id","gene","expression","SNP","genotype")

gene1$pair<-"copy 1"
gene2$pair<-"copy 2"


longgene<-rbind.data.frame(gene1,gene2,stringsAsFactors = FALSE)
longgene$name<-paste("Pair: ",sapply(strsplit(longgene$ID,"_Ssal"),FUN = `[[`, 1),sapply(strsplit(longgene$ID,"_Ssal"),FUN = `[[`, 2),"\nSNP: ",longgene$SNP,sep="")

ggplot(longgene,aes(x=as.factor(genotype),y=expression,fill=pair))+geom_boxplot()+theme_bw(20)+facet_wrap(~name,scales = "free_y",labeller = labeller(name = label_wrap_gen(width = 45)))+xlab("eQTL genotype")

longgene<-rbind.data.frame(gene1,gene2,sumgene,stringsAsFactors = FALSE)
ggplot(longgene,aes(x=as.factor(genotype),y=expression,fill=pair))+geom_boxplot()+theme_bw(25)+facet_wrap(~name,scales = "free_y",labeller = labeller(name = label_wrap_gen(width = 45)))+xlab("eQTL genotype")+scale_fill_manual(values=c("firebrick1","dodgerblue1","darkorchid1"))

truenegpair<-c("Pair: ENSSSAG00000011938_ENSSSAG00000108018\nSNP: ssa05_72178688","Pair: ENSSSAG00000086178_ENSSSAG00000087231\nSNP: ssa07_56580022","Pair: ENSSSAG00000113858_ENSSSAG00000120807\nSNP: ssa26_27976437")
newlonggene<-longgene[!(longgene$name%in%truenegpair),]

ggplot(newlonggene,aes(x=as.factor(genotype),y=expression,fill=pair))+geom_boxplot()+theme_bw(25)+facet_wrap(~name,scales = "free_y",labeller = labeller(name = label_wrap_gen(width = 45)))+xlab("eQTL genotype")+scale_fill_manual(values=c("firebrick1","dodgerblue1","darkorchid1"))

##for a positive pair

genetab<-read.table(file = "expression_gene_pair_pos01.txt",header=TRUE)
genetab$id=paste("0_",genetab$id,sep="")
genotab<-read.table(file = "genotype_gene_pair_pos01.txt",header=TRUE)

fulltabneg<-merge(genetab,genotab,by=c("ID","id"))


ggplot(fulltabneg[fulltabneg$ID=="ENSSSAG00000034990_Ssal_ENSSSAG00000071166_Ssal",],aes(x=Expression_gene1,y=Expression_gene2,color=as.factor(genotype)))+geom_point()+theme_bw(20)+ggtitle("correlation for pair ENSSSAG00000011938_Ssal_ENSSSAG00000108018_Ssal")


toplot<-ggplot(fulltabneg[fulltabneg$ID=="ENSSSAG00000034990_Ssal_ENSSSAG00000071166_Ssal",],aes(x=Expression_gene1,fill=as.factor(genotype)))+geom_density(alpha=0.5)+theme_bw(20)+theme(legend.position = "none")+ggtitle("correlation for proteasome 20S subunit alpha 3")+xlab("Expression copy 1")
rightplot<-ggplot(fulltabneg[fulltabneg$ID=="ENSSSAG00000034990_Ssal_ENSSSAG00000071166_Ssal",],aes(x=Expression_gene2,fill=as.factor(genotype)))+geom_density(alpha=0.5)+theme_bw(20)+coord_flip()+xlab("Expression copy 2")+ylim(0,0.03)
scatter<-ggplot(fulltabneg[fulltabneg$ID=="ENSSSAG00000034990_Ssal_ENSSSAG00000071166_Ssal",],aes(x=Expression_gene1,y=Expression_gene2,color=as.factor(genotype)))+geom_point()+theme_bw(20)+theme(legend.position = "none")+xlab("Expression copy 1")+ylab("Expression copy 2")
empty <- ggplot() + 
  geom_point(aes(1,1), colour="white") +
  theme(                              
    plot.background = element_blank(), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.border = element_blank(), 
    panel.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  )

library(egg) 

ggarrange(toplot, empty, scatter, rightplot, 
          ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))



##find genes position
gene1l<-unique(fulltabneg$gene1) ##position in order: 5: 87,782,557-87,817,196 / 17: 81,256,811-81,261,845 /  11: 38,167,566-38,174,827 
gene2l<-unique(fulltabneg$gene2) ##position in order : 2: 4,851,964-4,887,386 / 7: 62,703,781-62,708,799 / 26: 38,895,331-38,901,886 




## check snp
SNP_tab<-read.table(file="/mnt/SCRATCH/cedi/phDSalmon/sv_detection/filtervcf/Only_SNP_Farmed_Europe_filtered.recode.vcf.gz_1100000010.stat")

##sex of individuals
indsex<-read.table("synchro_parents.txt",header=TRUE,fill = TRUE)
males <- indsex[,1]
males <- males[males!="0"]
males <- paste(substr(males,1,9),substr(males,13,19),sep="")
females <- indsex[,2]
females <- females[females!=""]
females <- paste(substr(females,1,9),substr(females,13,19),sep="")
##make the name be the same that those in the sex table
##list of samples
#sampleeu<-read.table(file="/mnt/SCRATCH/cedi/phDSalmon/sv_detection/manta_sv_detect_11_11_2022/Farmed_European_sample.txt",header=FALSE)
sampleeu<-read.table(file="/mnt/SCRATCH/cedi/phDSalmon/sv_detection/filtervcf/sample_SNP_FE.txt",header=FALSE)
sampname<-sampleeu[,1]
sampname[1]<-"2014GNOMales_169_D03_RG"
sampname[2]<-"2014GNOMales_169_H12_RG"
sampname[3]<-"2014GNOMales_000_D02_RG"
samplename2<-paste(substr(sampname,1,4),"GNOS1",sapply(strsplit(sampname,"_"),FUN = `[[`, 2),sapply(strsplit(sampname,"_"),FUN = `[[`, 3),sep="")

males<-males[males%in%samplename2]
females<-females[females%in%samplename2]

colnames(SNP_tab)<-c("CHROM","POS",samplename2)

SNP_tab[SNP_tab$CHROM=="ssa26" & SNP_tab$POS==27976437,]



## xpehh from pauline
xpehhtab<-read.table(file="/mnt/SCRATCH/pabu/bedtools/input2/xp_ehh_euro.bed")

sharedeqtl<-results2[results2$regulation=="bothsamereg",]


rna_full_table<-read.table(file="Synchro_uniq_fulltable_withinfo_04.txt",header=TRUE)
rna_full_table<-rna_full_table[rna_full_table$condition=="LL",]
sharedgene<-rna_full_table[rna_full_table$gene%in%sharedeqtl$gene1 | rna_full_table$gene%in%sharedeqtl$gene2,]
sharedgene$snpid<-paste(sharedgene$snpschrom,sharedgene$snpstart,sep="_")

xpehhtab$snpid<-paste(xpehhtab$V1,xpehhtab$V2,sep="_")


#joinboth<-inner_join(sharedgene,xpehhtab,by=c("snpid"))

#joinboth<-merge(sharedgene,xpehhtab,by=c("snpid"))




##check long reads transcript quantity

quantitytable<-read.table(file="/mnt/SCRATCH/cedi/phDSalmon/rna_expression/WGD_geneexp/transcript_eif3_counts.tsv",header=TRUE)
gilltab<-quantitytable[,c(1,20:29)]

quantitytable<-read.table(file="/mnt/SCRATCH/cedi/phDSalmon/rna_expression/WGD_geneexp/transcript_tetraspanin_counts.tsv",header=TRUE)
gilltab<-quantitytable[,c(1,20:29)]

quantitytable<-read.table(file="/mnt/SCRATCH/cedi/phDSalmon/rna_expression/WGD_geneexp/transcript_collectin_counts.tsv",header=TRUE)
gilltab<-quantitytable[,c(1,20:29)]

quantitytable<-read.table(file="/mnt/SCRATCH/cedi/phDSalmon/rna_expression/WGD_geneexp/transcript_bothlof_counts.tsv",header=TRUE)
gilltab<-quantitytable[,c(1,20:29)]


gilltab[1,]



##check type of genes in each categorie
genepair<-read.table(file="genetabgenepair02.txt",header=TRUE)
genepairsamereg1<-genepair[genepair$regulation=="noreg","gene1"]
genepairsamereg2<-genepair[genepair$regulation=="noreg","gene2"]
genevector<-c(genepairsamereg1,genepairsamereg2)
genevector2<-as.data.frame(genevector)
write.table(genevector2,file="noeQTLgenes",row.names=FALSE,quote=FALSE,col.names = FALSE)



sharedeqtltab<-read.table(file="enrichment_sharedeQTL.csv",header=TRUE,sep=",")
sharedeqtltab$regulation<-"shared EQTL"
diffeqtltab<-read.table(file="enrichment_diff_eQTL.csv",header=TRUE,sep=",")
diffeqtltab$regulation<-"diff eQTL"
singeqtltab<-read.table(file="enrichment_single_eQTL.csv",header=TRUE,sep=",")
singeqtltab$regulation<-"single eQTL"
noeqtltab<-read.table(file="enrichment_noeQTL.csv",header=TRUE,sep=",")
noeqtltab$regulation<-"no eQTL"
wholeenrichment<-rbind.data.frame(sharedeqtltab,diffeqtltab,singeqtltab,noeqtltab,stringsAsFactors = FALSE)


library(ggplot2)
ggplot(wholeenrichment,aes(x=Pathway,y=-log10(Enrichment.FDR),color=regulation))+geom_point(size=wholeenrichment$Fold.Enrichment/5)+theme(axis.text.x = element_text(angle = 45, hjust=1))

wholeenrichment2<-wholeenrichment[wholeenrichment$Enrichment.FDR<0.0000000001,]

ggplot(wholeenrichment2,aes(x=Pathway,y=-log10(Enrichment.FDR),color=regulation))+geom_point(size=wholeenrichment2$Fold.Enrichment/5)+theme(axis.text.x = element_text(angle = 45, hjust=1))




## investigate lof pairs
singlelof<-c("ENSSSAG00000052119_Ssal_ENSSSAG00000059710_Ssal","ENSSSAG00000010241_Ssal_ENSSSAG00000042869_Ssal","ENSSSAG00000086243_Ssal_ENSSSAG00000109224_Ssal","ENSSSAG00000065265_Ssal_ENSSSAG00000087702_Ssal","ENSSSAG00000043807_Ssal_ENSSSAG00000096923_Ssal","ENSSSAG00000060532_Ssal_ENSSSAG00000067803_Ssal","ENSSSAG00000091200_Ssal_ENSSSAG00000121190_Ssal")
bothlof<-c("ENSSSAG00000058006_Ssal_ENSSSAG00000093636_Ssal","ENSSSAG00000084674_Ssal_ENSSSAG00000090558_Ssal","ENSSSAG00000089270_Ssal_ENSSSAG00000090739_Ssal")


#tablepair<-genepair[genepair$ID%in%bothlof,]
tablepair<-genepair[genepair$ID%in%singlelof,]
library(tidyverse)
##start of table making
Covariate <-  read_table2("covariate_synchro.txt") 
Expression <- read_table2("/mnt/project/MSLab/Marie/2022_Synchrosmolt/QTLTools/Synchro_mt_TPM.bed")    %>% select(-c(1,2,3,5,6))
Genotype <- read_table2("/mnt/project/MSLab/Marie/2022_Synchrosmolt/QTLTools/imputed_0.01_2730_smolts.traw")   %>% select(-c(1,3,4,5,6))

neggeneexpp1<- Expression[Expression$pid%in%tablepair$gene1,]
neggeneexpp2<- Expression[Expression$pid%in%tablepair$gene2,]

neggene_long1<- neggeneexpp1 %>%
  pivot_longer(cols = -pid, names_to = "id", values_to = "expression")

neggene_long2<- neggeneexpp2 %>%
  pivot_longer(cols = -pid, names_to = "id", values_to = "expression")



covar_kept<-Covariate[c(3,5),]
id<-colnames(covar_kept)
tcovar<-t(covar_kept)
long_covar<-cbind.data.frame(id,tcovar)
long_covar<-long_covar[-1,]
colnames(long_covar)<-c("id","sex","LL")

long_covar_LL<-long_covar[,c(1,3)]

#longcovarll<-long_covar_LL[long_covar_LL$LL=="LL",] ## 907 individuals



mergedpair1neg <- merge(neggene_long1 , long_covar, by = "id")
mergedpair2neg <- merge(neggene_long2 , long_covar, by = "id")

##filter

mergedpair1neg_filter<- mergedpair1neg[mergedpair1neg$LL=="LL",]
mergedpair2neg_filter<- mergedpair2neg[mergedpair2neg$LL=="LL",]

##can keep more column if wanted
mergedpair1neg_filter<-mergedpair1neg_filter[,1:3]
mergedpair2neg_filter<-mergedpair2neg_filter[,1:3]


library(dplyr)

# Left join table 1 with table 2 based on gene1 and pid
joined_table1_2neg <- left_join(genepair, mergedpair1neg_filter, by = c("gene1" = "pid"))
joined_table1_2neg<-joined_table1_2neg[!is.na(joined_table1_2neg$expression),]

# Left join result with table 3 based on gene2 and pid
joined_table1_2_3neg <- left_join(joined_table1_2neg, mergedpair2neg_filter, by = c("id"="id","gene2" = "pid"))
joined_table1_2_3neg <- joined_table1_2_3neg [!is.na(joined_table1_2_3neg$expression.y),]

colnames(joined_table1_2_3neg)[6:7]<-c("Expression_gene1","Expression_gene2")
write.table(joined_table1_2_3neg,file="expression_gene_pair_singlelof_01.txt",row.names=FALSE)


##get their genotype

rna_full_table<-read.table(file="Synchro_uniq_fulltable_withinfo_04.txt",header=TRUE)
rna_full_table<-rna_full_table[rna_full_table$condition=="LL",]
neggene<-rna_full_table[rna_full_table$gene%in%tablepair$gene1 | rna_full_table$gene%in%tablepair$gene2,]
neggene$snpchr2<-ifelse(neggene$snpschrom<10,paste("ssa0",neggene$snpschrom,sep=""),paste("ssa",neggene$snpschrom,sep=""))
neggene$snpid<-paste(neggene$snpchr2,neggene$snpstart,sep="_")

snpid<-unique(neggene$snpid)


gtkept<-Genotype[Genotype$SNP%in%snpid,]

geno_table_neg <- gtkept %>% pivot_longer(cols = -SNP, names_to = "id", values_to = "genotype")
geno_table_neg$id<-sapply(strsplit(geno_table_neg$id,"_"),FUN = `[[`, 2)

gene_tab_neg_merge <- merge(geno_table_neg  , long_covar, by = "id")  
gene_tab_neg_merge_filter<- gene_tab_neg_merge[gene_tab_neg_merge$LL=="LL",]

#pairdata<-data.frame(ID=c("ENSSSAG00000089270_Ssal_ENSSSAG00000090739_Ssal","ENSSSAG00000089270_Ssal_ENSSSAG00000090739_Ssal","ENSSSAG00000089270_Ssal_ENSSSAG00000090739_Ssal","ENSSSAG00000058006_Ssal_ENSSSAG00000093636_Ssal","ENSSSAG00000089270_Ssal_ENSSSAG00000090739_Ssal","ENSSSAG00000058006_Ssal_ENSSSAG00000093636_Ssal","ENSSSAG00000058006_Ssal_ENSSSAG00000093636_Ssal","ENSSSAG00000058006_Ssal_ENSSSAG00000093636_Ssal"),SNP=c("ssa09_18435847","ssa06_5855146","ssa06_16705377","ssa05_79024252","ssa03_97813453","ssa02_24479361","ssa02_18219589","ssa02_14165657"))
pairdata<-data.frame(ID=c("ENSSSAG00000052119_Ssal_ENSSSAG00000059710_Ssal","ENSSSAG00000010241_Ssal_ENSSSAG00000042869_Ssal","ENSSSAG00000010241_Ssal_ENSSSAG00000042869_Ssal","ENSSSAG00000086243_Ssal_ENSSSAG00000109224_Ssal","ENSSSAG00000086243_Ssal_ENSSSAG00000109224_Ssal","ENSSSAG00000065265_Ssal_ENSSSAG00000087702_Ssal","ENSSSAG00000065265_Ssal_ENSSSAG00000087702_Ssal","ENSSSAG00000065265_Ssal_ENSSSAG00000087702_Ssal","ENSSSAG00000060532_Ssal_ENSSSAG00000067803_Ssal","ENSSSAG00000060532_Ssal_ENSSSAG00000067803_Ssal","ENSSSAG00000060532_Ssal_ENSSSAG00000067803_Ssal","ENSSSAG00000091200_Ssal_ENSSSAG00000121190_Ssal","ENSSSAG00000091200_Ssal_ENSSSAG00000121190_Ssal","ENSSSAG00000091200_Ssal_ENSSSAG00000121190_Ssal","ENSSSAG00000043807_Ssal_ENSSSAG00000096923_Ssal","ENSSSAG00000043807_Ssal_ENSSSAG00000096923_Ssal","ENSSSAG00000043807_Ssal_ENSSSAG00000096923_Ssal"),SNP=c("ssa17_52896112","ssa10_43409869","ssa10_33543516","ssa17_53382503","ssa17_53426443","ssa06_4239457","ssa03_96416465","ssa03_72887980","ssa13_94213776","ssa12_99098559","ssa26_13562878","ssa05_80384971","ssa02_55830642","ssa02_23587078","ssa22_34446567","ssa14_48000689","ssa06_4543142"))

mergegt<-merge(gene_tab_neg_merge_filter,pairdata,by="SNP")

write.table(mergegt,file="genotype_gene_singlelof_01.txt",row.names=FALSE)


##plots

genetab<-read.table(file = "expression_gene_pair_bothlof_01.txt",header=TRUE)
genotab<-read.table(file = "genotype_gene_bothlof_01.txt",header=TRUE)

fulltabneg<-merge(genetab,genotab,by=c("ID","id"))


gene1<-fulltabneg[,c(1,2,3,6,8,9)]
gene2<-fulltabneg[,c(1,2,4,7,8,9)]

colnames(gene2)<-c("ID","id","gene","expression","SNP","genotype")
colnames(gene1)<-c("ID","id","gene","expression","SNP","genotype")

gene1$pair<-"gene1"
gene2$pair<-"gene2"

longgene<-rbind.data.frame(gene1,gene2,stringsAsFactors = FALSE)
longgene$name<-paste("Pair: ",sapply(strsplit(longgene$ID,"_Ssal"),FUN = `[[`, 1),sapply(strsplit(longgene$ID,"_Ssal"),FUN = `[[`, 2),"\nSNP: ",longgene$SNP,sep="")


ggplot(longgene,aes(x=as.factor(genotype),y=expression,fill=pair))+geom_boxplot()+theme_bw(25)+facet_wrap(~name,scales = "free_y",labeller = labeller(name = label_wrap_gen(width = 45)))




genetab<-read.table(file = "expression_gene_pair_singlelof_01.txt",header=TRUE)
genotab<-read.table(file = "genotype_gene_singlelof_01.txt",header=TRUE)

fulltabneg<-merge(genetab,genotab,by=c("ID","id"))


gene1<-fulltabneg[,c(1,2,3,6,8,9)]
gene2<-fulltabneg[,c(1,2,4,7,8,9)]

colnames(gene2)<-c("ID","id","gene","expression","SNP","genotype")
colnames(gene1)<-c("ID","id","gene","expression","SNP","genotype")

gene1$pair<-"gene1"
gene2$pair<-"gene2"

longgene<-rbind.data.frame(gene1,gene2,stringsAsFactors = FALSE)
longgene$name<-paste("Pair: ",sapply(strsplit(longgene$ID,"_Ssal"),FUN = `[[`, 1),sapply(strsplit(longgene$ID,"_Ssal"),FUN = `[[`, 2),"\nSNP: ",longgene$SNP,sep="")

##divide table in 2
names<-unique(longgene$name)
longene1<-longgene[longgene$name%in%names[1:9],]
longene2<-longgene[longgene$name%in%names[10:17],]
ggplot(longene1,aes(x=as.factor(genotype),y=expression,fill=pair))+geom_boxplot()+theme_bw(25)+facet_wrap(~name,scales = "free_y",labeller = labeller(name = label_wrap_gen(width = 45)))
ggplot(longene2,aes(x=as.factor(genotype),y=expression,fill=pair))+geom_boxplot()+theme_bw(25)+facet_wrap(~name,scales = "free_y",labeller = labeller(name = label_wrap_gen(width = 45)))

longgene$name<-paste("Pair",sapply(strsplit(longgene$ID,"_Ssal"),FUN = `[[`, 1),sapply(strsplit(longgene$ID,"_Ssal"),FUN = `[[`, 2),"SNP",longgene$SNP,sep="")

for (i in unique(longgene$name)){
  title<-i
  p<-ggplot(longgene[longgene$name==i,],aes(x=as.factor(genotype),y=expression,fill=pair))+geom_boxplot()+theme_bw(25)+ggtitle(title)
  fname<-paste("Pairplot",i,".png",sep="")
  ggsave(filename=fname,plot=p,device="png",path="/net/fs-2/scale/OrionStore/Home/cedi/phDSalmon/rna_analysis/",width=15,height=12,units="in")
}


## check similarity of the pair depending on their regulation status

## check the similarity 
library(fuzzyjoin)
library(dplyr)
library(IRanges)
library(ggplot2)


IDfile<-read.table(file="Homologous_blocks.csv",header=FALSE,sep=";")
IDfile<-IDfile[,1:4]
colnames(IDfile)<-c("number","identity","POS","CHROM")
IDfile$POS<-as.numeric(as.character(gsub("\\s+", "", IDfile$POS)))
IDfile$END<- ifelse(lead(IDfile$POS) > IDfile$POS, lead(IDfile$POS), IDfile$POS + 1000000)
IDfile[2582,"END"]<-42816660
IDfile$identity<- as.numeric(as.character(gsub(",", ".", IDfile$identity)))
IDfile$CHROM<-ifelse(substr(IDfile$CHROM,4,4)==0,substr(IDfile$CHROM,5,5),substr(IDfile$CHROM,4,5))

##remove start id duplicate
# Identify duplicate rows based on columns A and B
duplicated_rows <- duplicated(IDfile[, c("POS", "CHROM")])

# Subset the data frame to keep only non-duplicated rows
IDfile <- subset(IDfile, !duplicated_rows)




tablegene<-read.table(file="https://salmobase.org/datafiles/TSV/genes/AtlanticSalmon/Ssal_v3.1/Ensembl_genes.tsv",header=TRUE,sep="\t")
tablegene<-tablegene[!is.na(as.numeric(as.character(tablegene$seqname))),]

tablegene1<-tablegene
colnames(tablegene1)[5]<-"gene1"

tablegene2<-tablegene
colnames(tablegene2)[5]<-"gene2"

genepair<-read.table(file="genetabgenepair02.txt",header=TRUE)


tablegenwithpos1<-merge(genepair,tablegene1,by=c("gene1"))
tablegenwithpos2<-merge(tablegenwithpos1,tablegene2,by=c("gene2"))

##by doing this we goes from 11233 gene pair to 1771

tablegenwithpos<-tablegenwithpos2[,c(1:7,15:17)]
colnames(tablegenwithpos)[5:10]<-c("CHROM1","START1","END1","CHROM2","START2","END2")



for( i in 1:length(tablegenwithpos$gene2)){
  chrom1<-tablegenwithpos[i,"CHROM1"]
  pos1<-(tablegenwithpos[i,"START1"]+tablegenwithpos[i,"END1"])/2 ##take the middle of the gene
  similarity1<-IDfile[IDfile$CHROM==chrom1 & IDfile$POS<pos1 & IDfile$END>pos1,"identity"]
  if(length(similarity1)==0){
    similarity1<-NA
  }
  chrom2<-tablegenwithpos[i,"CHROM2"]
  pos2<-(tablegenwithpos[i,"START2"]+tablegenwithpos[i,"END2"])/2 ##take the middle of the gene
  similarity2<-IDfile[IDfile$CHROM==chrom2 & IDfile$POS<pos2 & IDfile$END>pos2,"identity"]
  if(length(similarity2)==0){
    similarity2<-NA
  }
  tablegenwithpos[i,"similarity1"]<-similarity1
  tablegenwithpos[i,"similarity2"]<-similarity2
}
tablegenwithposnona<-tablegenwithpos[!is.na(tablegenwithpos$similarity1),]
tablegenwithposnona<-tablegenwithposnona[!is.na(tablegenwithposnona$similarity2),]

tablegenwithposnona$meansimilarity<-(tablegenwithposnona$similarity1+tablegenwithposnona$similarity2)/2

ggplot(tablegenwithposnona,aes(x=similarity1,y=similarity2))+geom_point() ##some outliers where similarity 1 is very different from similarity 2

ggplot(tablegenwithposnona,aes(x=meansimilarity,fill=regulation))+geom_histogram(color="black")+theme_bw(20)+facet_wrap(~regulation,scales = 'free_y')




### plot of B3 (figure with detailled scenarios of pairs of shared eQTL)
genetab<-read.table(file="genetabgenepairwithtrans01.txt",header=TRUE)

genetabtest<-genetab[genetab$regulation=="bothsamereg",]

genetabb3<-genetab[genetab$regulation=="bothsamereg" & genetab$reegscenario=="inter_chr_non-ohnolog#inter_chr_non-ohnolog" & genetab$regtype=="both_both",]

Expression <- read_table2("/mnt/project/MSLab/Marie/2022_Synchrosmolt/QTLTools/Synchro_mt_TPM.bed")    %>% select(-c(1,2,3,5,6))
Genotype <- read_table2("/mnt/project/MSLab/Marie/2022_Synchrosmolt/QTLTools/imputed_0.01_2730_smolts.traw")   %>% select(-c(1,3,4,5,6))
Covariate <-  read_table2("covariate_synchro.txt") 
## fix the name in genotype
names(Genotype)[2:ncol(Genotype)] <- str_remove(names(Genotype)[2:ncol(Genotype)], "0_")


## Get expression

neggeneexpp1<- Expression[Expression$pid%in%genetabb3$gene1,]
neggeneexpp2<- Expression[Expression$pid%in%genetabb3$gene2,]

neggene_long1<- neggeneexpp1 %>%
  pivot_longer(cols = -pid, names_to = "id", values_to = "expression")

neggene_long2<- neggeneexpp2 %>%
  pivot_longer(cols = -pid, names_to = "id", values_to = "expression")



covar_kept<-Covariate[c(3,5),]
id<-colnames(covar_kept)
tcovar<-t(covar_kept)
long_covar<-cbind.data.frame(id,tcovar)
long_covar<-long_covar[-1,]
colnames(long_covar)<-c("id","sex","LL")

long_covar_LL<-long_covar[,c(1,3)]

#longcovarll<-long_covar_LL[long_covar_LL$LL=="LL",] ## 907 individuals



mergedpair1neg <- merge(neggene_long1 , long_covar, by = "id")
mergedpair2neg <- merge(neggene_long2 , long_covar, by = "id")

##filter

mergedpair1neg_filter<- mergedpair1neg[mergedpair1neg$LL=="LL",]
mergedpair2neg_filter<- mergedpair2neg[mergedpair2neg$LL=="LL",]

##can keep more column if wanted
mergedpair1neg_filter<-mergedpair1neg_filter[,1:3]
mergedpair2neg_filter<-mergedpair2neg_filter[,1:3]


library(dplyr)

# Left join table 1 with table 2 based on gene1 and pid
joined_table1_2neg <- left_join(genepair, mergedpair1neg_filter, by = c("gene1" = "pid"))
joined_table1_2neg<-joined_table1_2neg[!is.na(joined_table1_2neg$expression),]

# Left join result with table 3 based on gene2 and pid
joined_table1_2_3neg <- left_join(joined_table1_2neg, mergedpair2neg_filter, by = c("id"="id","gene2" = "pid"))
joined_table1_2_3neg <- joined_table1_2_3neg [!is.na(joined_table1_2_3neg$expression.y),]

colnames(joined_table1_2_3neg)[6:7]<-c("Expression_gene1","Expression_gene2")
write.table(joined_table1_2_3neg,file="expression_gene_pair_b301.txt",row.names=FALSE)



##get genotype

rna_full_table<-read.table(file="Synchro_uniq_fulltable_withinfo_04.txt",header=TRUE)
rna_full_table<-rna_full_table[rna_full_table$condition=="LL",]
neggene<-rna_full_table[rna_full_table$gene%in%genetabb3$gene1 | rna_full_table$gene%in%genetabb3$gene2,]
neggene$snpchr2<-ifelse(neggene$snpschrom<10,paste("ssa0",neggene$snpschrom,sep=""),paste("ssa",neggene$snpschrom,sep=""))
neggene$snpid<-paste(neggene$snpchr2,neggene$snpstart,sep="_")

## trans part

neggenetrans<-neggene[neggene$cistrans=="trans",]

snpid<-unique(neggenetrans$snpid)

gtkept<-Genotype[Genotype$SNP%in%snpid,]

geno_table_neg <- gtkept %>% pivot_longer(cols = -SNP, names_to = "id", values_to = "genotype")

gene_tab_neg_merge <- merge(geno_table_neg  , long_covar, by = "id")  
gene_tab_neg_merge_filter<- gene_tab_neg_merge[gene_tab_neg_merge$LL=="LL",]

pairdata<-data.frame(ID=c("ENSSSAG00000008620_Ssal_ENSSSAG00000044419_Ssal","ENSSSAG00000041645_Ssal_ENSSSAG00000070444_Ssal","ENSSSAG00000068682_Ssal_ENSSSAG00000075896_Ssal","ENSSSAG00000067127_Ssal_ENSSSAG00000072816_Ssal","ENSSSAG00000043064_Ssal_ENSSSAG00000079737_Ssal","ENSSSAG00000055063_Ssal_ENSSSAG00000096396_Ssal","ENSSSAG00000008406_Ssal_ENSSSAG00000118613_Ssal"),SNP=c("ssa14_50870510","ssa14_19961522","ssa10_113545511","ssa14_48007951","ssa14_49553986","ssa14_50815335","ssa14_32499888"))


mergegt<-merge(gene_tab_neg_merge_filter,pairdata,by="SNP")

write.table(mergegt,file="genotype_gene_pair_b3_trans_01.txt",row.names=FALSE)


##cis part
neggenecis<-neggene[neggene$cistrans=="cis",]

snpid<-unique(neggenecis$snpid)

gtkept<-Genotype[Genotype$SNP%in%snpid,]

geno_table_neg <- gtkept %>% pivot_longer(cols = -SNP, names_to = "id", values_to = "genotype")

gene_tab_neg_merge <- merge(geno_table_neg  , long_covar, by = "id")  
gene_tab_neg_merge_filter<- gene_tab_neg_merge[gene_tab_neg_merge$LL=="LL",]

genedata<-neggenecis[,c("gene","snpid")]
colnames(genedata)<-c("gene","SNP")

mergegt<-merge(gene_tab_neg_merge_filter,genedata,by="SNP")

write.table(mergegt,file="genotype_gene_pair_b3_cis_01.txt",row.names=FALSE)


## make the plot

expressiontab<-read.table(file="expression_gene_pair_b301.txt",header=TRUE)
transgenotype<-read.table(file="genotype_gene_pair_b3_trans_01.txt",header=TRUE)
cisgenotype<-read.table(file="genotype_gene_pair_b3_cis_01.txt",header=TRUE)


mergedexpr_trans<-merge(expressiontab,transgenotype,by=c("id","ID"))

mergedcisgene1<-merge(mergedexpr_trans,cisgenotype,by.x=c("id","gene1"),by.y=c("id","gene"))
mergedcisgene1<-mergedcisgene1[,-c(14,15)]
colnames(mergedcisgene1)[8:13]<-c("trans_eQTL","trans_genotype","sex","LL","cis_eQTL","cis_genotype")

mergedcisgene2<-merge(mergedexpr_trans,cisgenotype,by.x=c("id","gene2"),by.y=c("id","gene"))
mergedcisgene2<-mergedcisgene2[,-c(14,15)]
colnames(mergedcisgene2)[8:13]<-c("trans_eQTL","trans_genotype","sex","LL","cis_eQTL","cis_genotype")


mergedcisgene1$copy<-"copy1"
mergedcisgene2$copy<-"copy2"

fulltab<-rbind.data.frame(mergedcisgene1,mergedcisgene2,stringsAsFactors = FALSE)

library(ggplot2)

#ggplot(fulltab[fulltab$ID=="ENSSSAG00000008406_Ssal_ENSSSAG00000118613_Ssal",],aes(x=as.factor(trans_genotype),y=Expression_gene1,fill=as.factor(cis_genotype)))+geom_boxplot()+facet_wrap(~copy)


mergedcisgene1<-mergedcisgene1[,-c(7)]
mergedcisgene2<-mergedcisgene2[,-c(6)]

colnames(mergedcisgene1)[6]<-"Expression"
colnames(mergedcisgene2)[6]<-"Expression"
fulltab<-rbind.data.frame(mergedcisgene1,mergedcisgene2,stringsAsFactors = FALSE)

#ggplot(fulltab[fulltab$ID=="ENSSSAG00000008406_Ssal_ENSSSAG00000118613_Ssal",],aes(x=as.factor(trans_genotype),y=Expression,fill=as.factor(cis_genotype)))+geom_boxplot()+facet_wrap(~copy)

for (i in unique(fulltab$ID)){
  title<-i
  p<-ggplot(fulltab[fulltab$ID==i,],aes(x=as.factor(trans_genotype),y=Expression,fill=as.factor(cis_genotype)))+geom_boxplot()+facet_wrap(~copy)
  fname<-paste("B3_pair",i,".png",sep="")
  ggsave(filename=fname,plot=p,device="png",path="/net/fs-2/scale/OrionStore/Home/cedi/phDSalmon/rna_analysis/",width=15,height=12,units="in")
}

#unique(fulltab$ID)
#ggplot(fulltab[fulltab$ID=="ENSSSAG00000008620_Ssal_ENSSSAG00000044419_Ssal",],aes(x=as.factor(trans_genotype),y=Expression,fill=as.factor(cis_genotype)))+geom_boxplot()+facet_wrap(~copy)




## refine C scenario
genepair<-read.table(file="genetabgenepair02.txt",header=TRUE)

Cscenario<-genepair[genepair$regulation=="bothdiffreg",]


##find out which pairs have eQTL on the same chromosome
rna_full_table<-read.table(file="Synchro_uniq_fulltable_withinfo_04.txt",header=TRUE)
rna_full_table<-rna_full_table[rna_full_table$condition=="LL",]

toberefined<-c()
for(i in 1:length(Cscenario$gene1)){
  gene1<-Cscenario[i,"gene1"]
  gene2<-Cscenario[i,"gene2"]
  tabpair2<-rna_full_table[rna_full_table$gene==gene2,]
  tabpair1<-rna_full_table[rna_full_table$gene==gene1,]
  for(y in unique(tabpair1$snpschrom)){
    if(y%in%tabpair2$snpschrom){
      pair<-Cscenario[i,"ID"]
      subsettab1<-tabpair1[tabpair1$snpschrom==y,]
      subsettab2<-tabpair2[tabpair2$snpschrom==y,]
      subtabtemp<-rbind.data.frame(subsettab1,subsettab2,stringsAsFactors = FALSE)
      subtabtemp$pair<-pair
      toberefined<-rbind.data.frame(toberefined,subtabtemp,stringsAsFactors = FALSE)
    }
  }
}

length(unique(toberefined$pair)) ##417 pair out of 1834 with snps on same chromosome
hist(toberefined$snpschrom,breaks = 29)

#write.table(toberefined,file="candidatepairinCscnenariotoberefined.txt",row.names = FALSE)

toberefined<-read.table("candidatepairinCscnenariotoberefined.txt",header=TRUE)
toberefined$chrom<-ifelse(toberefined$snpschrom<10,paste("ssa0",toberefined$snpschrom,sep=""),paste("ssa",toberefined$snpschrom,sep=""))

completetab<-c()
for( i in unique(toberefined$pair)){
  subset<-toberefined[toberefined$pair==i,]
  subset$coord<-paste(subset$chrom,subset$snpstart,sep="_")
  combinations <- combn(subset$coord, 2)
  new_table <- data.frame(Column1 = combinations[1, ], Column2 = combinations[2, ],pair=i)
  completetab<-rbind.data.frame(completetab,new_table)
}

write.table(completetab,file="/mnt/SCRATCH/cedi/phDSalmon/rna_expression/WGD_geneexp/refinedpair.txt",row.names = FALSE,col.names = FALSE,quote = FALSE)



haptab<-data.frame(toberefined$chrom,toberefined$snpstart)
write.table(haptab,file="/mnt/SCRATCH/cedi/phDSalmon/rna_expression/WGD_geneexp/refinedcscenario.txt",row.names = FALSE,col.names = FALSE,quote = FALSE)

haptab2<-data.frame(toberefined$chrom,toberefined$snpstart,toberefined$pair)
haptab2$posid<-paste(haptab2$toberefined.chrom,haptab2$toberefined.snpstart,sep="_")

## import ld results

ldfile<-read.table(file="/mnt/SCRATCH/cedi/phDSalmon/rna_expression/WGD_geneexp/ld_cscenario.list.hap.ld",header=TRUE)
ldfile$queryid<-paste(ldfile$CHR1,ldfile$POS1,sep="_")
ldfile$resultid<-paste(ldfile$CHR2,ldfile$POS2,sep="_")

trimmed_candidate<-c()
for(i in unique(haptab2$toberefined.pair)){
  subtab<-haptab2[haptab2$toberefined.pair==i,]
  ldres<-ldfile[ldfile$queryid%in%subtab$posid,]
  for(y in subtab$posid){
    if(y%in%ldres$resultid){
      ressave<-ldres[ldres$resultid==y,]
      ressave$pair<-i
      trimmed_candidate<-rbind.data.frame(trimmed_candidate,ressave)
    }
  }
}

toberefined[toberefined$pair=="ENSSSAG00000049938_Ssal_ENSSSAG00000092379_Ssal",]


##ld result version plink

ldres<-read.table(file="/mnt/SCRATCH/cedi/phDSalmon/rna_expression/WGD_geneexp/LD_plink.txt",header=TRUE)
#summary(ldres$LD)

library(ggplot2)
ggplot(ldres,aes(x=LD))+geom_histogram()
library(tidyverse)

##get back gene pair info
#completetab$mergedpos<-paste(completetab$Column1,completetab$Column2,sep="_")
#ldres$mergedpos<-paste(ldres$POS1,ldres$POS2,sep="_")

ldres2<-merge(completetab,ldres,by.x = c("Column1","Column2"),by.y = c("POS1","POS2"))

ldres2<-ldres2[!duplicated(ldres2),]

pairtokeep<-ldres2[ldres2$LD>0.8,"pair"]
pairtokeep<-unique(pairtokeep)



## identify the shared eQTL for the 89 pairs
genepair<-read.table(file="genetabgenepair02.txt",header=TRUE)
sharedtab<-genepair[genepair$regulation=="bothsamereg",]

rna_full_table<-read.table(file="Synchro_uniq_fulltable_withinfo_04.txt",header=TRUE)
rna_full_table<-rna_full_table[rna_full_table$condition=="LL",]
rna_full_table$snpid<-paste(rna_full_table$snpschrom,rna_full_table$snpstart,sep="_")


for (i in 1:length(sharedtab$gene2)){
  gene1<-sharedtab[i,"gene1"]
  gene2<-sharedtab[i,"gene2"]
  subtab1<-rna_full_table[rna_full_table$gene==gene1,]
  subtab2<-rna_full_table[rna_full_table$gene==gene2,]
  sharedeQTL<-subtab1[subtab1$snpid%in%subtab2$snpid,"snpid"]
  sharedtab[i,"eQTL"]<-sharedeQTL
}

sharedeqtlrnatab<-rna_full_table[rna_full_table$snpid%in%sharedtab$eQTL,]

##some count and plots
library(ggplot2)
library(tidyverse)

dfcount <- sharedeqtlrnatab %>% group_by(snpid) %>% 
  summarise(total_count=n(),.groups = 'drop') %>%
  as.data.frame()

ggplot(dfcount,aes(x=total_count,fill=..x..))+geom_histogram(color="black",binwidth = 10)+scale_fill_gradient("Number of genes regulated",low = "dodgerblue1", high = "firebrick")+scale_y_continuous(breaks = seq(0, max(dfcount$total_count), by = 2))+theme_bw(20)


dfcount <- sharedtab %>% group_by(eQTL) %>% 
  summarise(total_count=n(),.groups = 'drop') %>%
  as.data.frame()

ggplot(dfcount,aes(x=total_count,fill=..x..))+geom_histogram(color="black",binwidth = 1)+scale_fill_gradient("Number of pair regulated",low = "dodgerblue1", high = "firebrick")+theme_bw(20)+scale_y_continuous(breaks = seq(0, 35, by = 5))

length(dfcount[dfcount$total_count==1,"total_count"])

##investigate if each eQTL is cis or trans or both (compute the number of cis gene regulated, the number of trans gene regulated)

eqtltab<-c()
for(i in unique(sharedeqtlrnatab$snpid)){
  cistab<-sharedeqtlrnatab[sharedeqtlrnatab$snpid==i & sharedeqtlrnatab$cistrans=="cis","cistrans"]
  transtab<-sharedeqtlrnatab[sharedeqtlrnatab$snpid==i & sharedeqtlrnatab$cistrans=="trans","cistrans"]
  nbcis<-length(cistab)
  nbtrans<-length(transtab)
  if(nbcis==0 & nbtrans>0){
    status<-"trans"
  }else if(nbcis>0 & nbtrans==0){
    status<-"cis"
  }else{
    status<-"both"
  }
  row<-c(i,nbcis,nbtrans,status)
  eqtltab<-rbind.data.frame(eqtltab,row,stringsAsFactors = FALSE)
}
colnames(eqtltab)<-c("snpid","ciscount","transcount","status")

eqtltab$snpchrom<-sapply(strsplit(eqtltab$snpid,"_"),FUN = `[[`, 1)
eqtltab$snppos<-sapply(strsplit(eqtltab$snpid,"_"),FUN = `[[`, 2)

##investigate their similarity score

library(fuzzyjoin)
library(dplyr)
library(IRanges)



IDfile<-read.table(file="Homologous_blocks.csv",header=FALSE,sep=";")
IDfile<-IDfile[,1:4]
colnames(IDfile)<-c("number","identity","POS","CHROM")
IDfile$POS<-as.numeric(as.character(gsub("\\s+", "", IDfile$POS)))
IDfile$END<- ifelse(lead(IDfile$POS) > IDfile$POS, lead(IDfile$POS), IDfile$POS + 1000000)
IDfile[2582,"END"]<-42816660
IDfile$identity<- as.numeric(as.character(gsub(",", ".", IDfile$identity)))
IDfile$CHROM<-ifelse(substr(IDfile$CHROM,4,4)==0,substr(IDfile$CHROM,5,5),substr(IDfile$CHROM,4,5))

for(i in 1:length(eqtltab$snpid)){
  poseqtl<-eqtltab[i,"snppos"]
  chromeqtl<-eqtltab[i,"snpchrom"]
  identity<-IDfile[IDfile$CHROM==chromeqtl & IDfile$POS<=poseqtl & IDfile$END>poseqtl,"identity"]
  eqtltab[i,"identity"]<-identity[1]
}

##investigate similarity score of gene pair regulated

genepair<-read.table(file="genetabgenepair02.txt",header=TRUE)
sharedtab<-genepair[genepair$regulation=="bothsamereg",]

rna_full_table<-read.table(file="Synchro_uniq_fulltable_withinfo_04.txt",header=TRUE)
rna_full_table<-rna_full_table[rna_full_table$condition=="LL",]
rna_full_table$snpid<-paste(rna_full_table$snpschrom,rna_full_table$snpstart,sep="_")

rna_full_table_gene1_tot<-rna_full_table[rna_full_table$gene%in%sharedtab$gene1,]
rna_full_table_gene2_tot<-rna_full_table[rna_full_table$gene%in%sharedtab$gene2,]

rna_full_table_gene1<-rna_full_table_gene1_tot[,2:4]
rna_full_table_gene1_nopos<-rna_full_table_gene1_tot[,2:3]
rna_full_table_gene1<-rna_full_table_gene1[!duplicated(rna_full_table_gene1_nopos),]

rna_full_table_gene2<-rna_full_table_gene2_tot[,2:4]
rna_full_table_gene2_nopos<-rna_full_table_gene2_tot[,2:3]
rna_full_table_gene2<-rna_full_table_gene2[!duplicated(rna_full_table_gene2_nopos),]


for(i in 1:length(rna_full_table_gene1$gene)){
  posgene<-rna_full_table_gene1[i,"gene_start"]
  chromeqtl<-rna_full_table_gene1[i,"gene_chr"]
  identity<-IDfile[IDfile$CHROM==chromeqtl & IDfile$POS<=posgene & IDfile$END>posgene,"identity"]
  rna_full_table_gene1[i,"identity"]<-identity[1]
}

for(i in 1:length(rna_full_table_gene2$gene)){
  posgene<-rna_full_table_gene2[i,"gene_start"]
  chromeqtl<-rna_full_table_gene2[i,"gene_chr"]
  identity<-IDfile[IDfile$CHROM==chromeqtl & IDfile$POS<=posgene & IDfile$END>posgene,"identity"]
  rna_full_table_gene2[i,"identity"]<-identity[1]
}

gene1id<-data.frame(status=rep("copy1",89),identity=rna_full_table_gene1$identity)
gene2id<-data.frame(status=rep("copy2",89),identity=rna_full_table_gene2$identity)
mergeddf<-rbind.data.frame(nonsharedsubset,gene1id,gene2id,stringsAsFactors = FALSE)

ggplot(mergeddf,aes(y=identity,x=status,fill=status))+geom_boxplot()+theme_bw(20)


##investigate if they are in LOR or AOR

duplication_table<-read.table(file="https://salmobase.org/datafiles/TSV/synteny/2021-11/AtlanticSalmon/synteny.tsv",header=TRUE,sep="\t")
duplication_table$region<-ifelse(duplication_table$lore_count>duplication_table$aore_count,"LORe","AORe")

for(i in 1:length(eqtltab$snpid)){
  poseqtl<-as.numeric(as.character(eqtltab[i,"snppos"]))
  chromeqtl<-as.numeric(as.character(eqtltab[i,"snpchrom"]))
  
  subsetx<-duplication_table[duplication_table$chromosome_x==chromeqtl & duplication_table$begin_x<poseqtl & duplication_table$end_x>poseqtl,"region"]
  subsety<-duplication_table[duplication_table$chromosome_y==chromeqtl & duplication_table$begin_y<poseqtl & duplication_table$end_y>poseqtl,"region"]
  
  if(length(subsetx)!=0){
    region<-subsetx
  }else if(length(subsety)!=0){
    region<-subsety
  }else{
    region<-"Unknown"
  }
  eqtltab[i,"region"]<-region
}

#write.table(eqtltab,"eqtlshared_01.txt",row.names = FALSE)


##some plots
eqtltab<-read.table("eqtlshared_01.txt",header=TRUE)

dfcount <- eqtltab %>% group_by(region) %>% 
  summarise(total_count=n(),.groups = 'drop') %>%
  as.data.frame()

ggplot(dfcount,aes(x=region,y=total_count,fill=region))+geom_bar(stat="identity")+theme_bw(20)

#ggplot(eqtltab,aes(x=identity,fill=..x..))+geom_histogram(color="black")+theme_bw(20)+scale_fill_gradient("Identity of the region compared to its duplicate",low = "dodgerblue1", high = "firebrick")


##compute identity for all others snps
nonsharedeQTL<-rna_full_table[!(rna_full_table$snpid%in%sharedtab$eQTL),]

library(fuzzyjoin)
library(dplyr)
library(IRanges)


IDfile<-read.table(file="Homologous_blocks.csv",header=FALSE,sep=";")
IDfile<-IDfile[,1:4]
colnames(IDfile)<-c("number","identity","POS","CHROM")
IDfile$POS<-as.numeric(as.character(gsub("\\s+", "", IDfile$POS)))
IDfile$END<- ifelse(lead(IDfile$POS) > IDfile$POS, lead(IDfile$POS), IDfile$POS + 1000000)
IDfile[2582,"END"]<-42816660
IDfile$identity<- as.numeric(as.character(gsub(",", ".", IDfile$identity)))
IDfile$CHROM<-ifelse(substr(IDfile$CHROM,4,4)==0,substr(IDfile$CHROM,5,5),substr(IDfile$CHROM,4,5))

for(i in 1:length(nonsharedeQTL$snpid)){
  poseqtl<-nonsharedeQTL[i,"snpstart"]
  chromeqtl<-nonsharedeQTL[i,"snpschrom"]
  identity<-IDfile[IDfile$CHROM==chromeqtl & IDfile$POS<=poseqtl & IDfile$END>poseqtl,"identity"]
  nonsharedeQTL[i,"identity"]<-identity[1]
}

nonsharedsubset<-data.frame(status="non shared",identity=nonsharedeQTL$identity)
sharedsubset<-data.frame(status="shared",identity=eqtltab$identity)

mergeddf<-rbind.data.frame(nonsharedsubset,sharedsubset,stringsAsFactors = FALSE)

ggplot(mergeddf,aes(y=identity,x=status,fill=status))+geom_boxplot()+theme_bw(20)

mean(mergeddf[mergeddf$status=="shared","identity"])
mean(mergeddf[mergeddf$status=="non shared","identity"],na.rm=TRUE)


## cis and trans plots

#dfcount <- eqtltab %>% group_by(status) %>% 
#summarise(total_count=n(),.groups = 'drop') %>%
#as.data.frame()
#dfcount[3,]<-c("cis",0)
#dfcount$total_count<-as.numeric(as.character(dfcount$total_count))

#ggplot(dfcount,aes(x=status,y=total_count,fill=status))+geom_bar(stat="identity")+theme_bw(20)+scale_y_continuous(breaks = seq(0, 25, by = 5))

eqtltab$transminucscis<-eqtltab$transcount-eqtltab$ciscount

subtab<-eqtltab[eqtltab$status=="both",]
ggplot(subtab, aes(x = as.factor(snppos), y =transminucscis))+geom_point(alpha = 0.75,size=subtab$ciscount)+theme_bw(20)


#sharedeqtlrnatab$newgene_chr<-ifelse(sharedeqtlrnatab$gene_chr<10,paste("ssa0",sharedeqtlrnatab$gene_chr,sep=""),paste("ssa",sharedeqtlrnatab$gene_chr,sep=""))
#sharedeqtlrnatab$neweqtl_chr<-ifelse(sharedeqtlrnatab$snpschrom<10,paste("ssa0",sharedeqtlrnatab$snpschrom,sep=""),paste("ssa",sharedeqtlrnatab$snpschrom,sep=""))

##same as Dom plot

lengthvector<-c(174498729,95481959,105780080,90536438,92788608,96060288,68862998,28860523,161282225,125877811,111868677,101677876,
                114417674,101980477,110670232,96486271,87489397,84084598,88107222,96847506,59819933,63823863,52460201,49354470,54385492,
                55994222,45305548,41468476,43051128)
chromadd<-c(0,cumsum(lengthvector+50000000)[1:28])
#chromname<-unique(sharedeqtlrnatab$gene_chr)
chromname<-c(1:29)
chrommid<-chromadd+lengthvector/2
listval<-cbind.data.frame(chromname,lengthvector,chromadd,chrommid)
sharedeqtlrnatab$newposgene<-sharedeqtlrnatab$gene_start+listval[match(sharedeqtlrnatab$gene_chr,listval$chromname),3]
sharedeqtlrnatab$newposeqtl<-sharedeqtlrnatab$snpstart+listval[match(sharedeqtlrnatab$snpschrom,listval$chromname),3]

ggplot(sharedeqtlrnatab,aes(x=newposeqtl,y=newposgene,color=as.factor(gene_chr)))+ scale_color_manual(values = rep(c("#276FBF", "#183059"), 29))+geom_point(alpha = 0.75,size=0.80)+scale_x_continuous(label =listval$chromname, breaks = listval$chrommid )+scale_y_continuous(label =listval$chromname, breaks = listval$chrommid )+theme_bw(20)


## see effect size
effectszie<-read.table(file="LL_cistrans_slope_Aug2023.txt",header=TRUE)

genepair<-read.table(file="genetabgenepair02.txt",header=TRUE)
sharedtab<-genepair[genepair$regulation=="bothsamereg",]

rna_full_table<-read.table(file="Synchro_uniq_fulltable_withinfo_04.txt",header=TRUE)
rna_full_table<-rna_full_table[rna_full_table$condition=="LL",]
rna_full_table$snpid<-paste(ifelse(rna_full_table$snpschrom<9,paste("ssa0",rna_full_table$snpschrom,sep=""),paste("ssa",rna_full_table$snpschrom,sep="")),rna_full_table$snpstart,sep="_")

for (i in unique(sharedtab$gene1)){
  gene2<-sharedtab[sharedtab$gene1==i,"gene2"]
  snpidgene1<-rna_full_table[rna_full_table$gene==i,"snpid"]
  snpidgene2<-rna_full_table[rna_full_table$gene==gene2,"snpid"]
  causalsnp<-intersect(snpidgene1, snpidgene2)
  sharedtab[sharedtab$gene1==i,"snpid"]<-causalsnp
}

effectsziegene1<-merge(sharedtab,effectszie,by.x="gene1",by.y="gene")
effectsziegene2<-merge(sharedtab,effectszie,by.x="gene2",by.y="gene")

effectsizegene1<-effectsziegene1[effectsziegene1$snpid==effectsziegene1$SNP,]
effectsizegene2<-effectsziegene2[effectsziegene2$snpid==effectsziegene2$SNP,]

effectsizegene1$copy<-"copy1"
effectsizegene2$copy<-"copy2"

effectszieshared<-merge(effectsizegene1,effectsizegene2,by="ID")
effectszieshared$snpchrom<-sapply(strsplit(effectszieshared$SNP.x,"_"),FUN = `[[`, 1)
effectszieshared$negativecorr<-"no"
effectszieshared[effectszieshared$ID=="ENSSSAG00000011938_Ssal_ENSSSAG00000108018_Ssal","negativecorr"]<-"yes"
effectszieshared[effectszieshared$ID=="ENSSSAG00000086178_Ssal_ENSSSAG00000087231_Ssal","negativecorr"]<-"yes"
effectszieshared[effectszieshared$ID=="ENSSSAG00000113858_Ssal_ENSSSAG00000120807_Ssal","negativecorr"]<-"yes"
effectszieshared<-effectszieshared[,-(9:13)]

library(ggplot2)

ggplot(effectszieshared,aes(x=slope.x,y=slope.y,color=negativecorr))+geom_point(size=2)+theme_bw(20)

ggplot(effectszieshared,aes(x=slope.x,y=slope.y,color=negativecorr))+geom_point(size=2)+theme_bw(20)+geom_smooth(method=lm, color='#2C3E50')+ geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")

ggplot(effectszieshared,aes(x=abs(slope.x),y=abs(slope.y),color=negativecorr))+geom_point(size=2)+theme_bw(20)

ggplot(effectszieshared,aes(x=abs(slope.x),y=abs(slope.y),color=negativecorr))+geom_point(size=2)+theme_bw(20)+geom_smooth(method=lm, color='#2C3E50')+ geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")


##investigate the others outliers...
genepairoutlier<-data.frame(ID=c("ENSSSAG00000087643_Ssal_ENSSSAG00000112530_Ssal","ENSSSAG00000053979_Ssal_ENSSSAG00000074975_Ssal"),gene1=c("ENSSSAG00000112530","ENSSSAG00000053979"),gene2=c("ENSSSAG00000087643","ENSSSAG00000074975"),snpid=c("ssa03_68202430","ssa17_58675294"))

#write.table(genepairoutlier,file="geneoutlier.tab",row.names = FALSE)

tableexpression<-Get_SNP_gene_info(genepairoutlier,feature="pair")

tableexpression$sumexpression<-tableexpression$Expression_gene1+tableexpression$Expression_gene2

gene1<-tableexpression[,c(1,2,3,6,8,9)]
gene2<-tableexpression[,c(1,2,4,7,8,9)]

sumgene<-data.frame(ID=tableexpression$ID,id=tableexpression$id,gene=tableexpression$ID,expression=tableexpression$sumexpression,SNP=tableexpression$SNP,genotype=tableexpression$genotype)
colnames(sumgene)<-c("ID","id","gene","expression","SNP","genotype")
sumgene$pair<-"Sum of copies"

colnames(gene2)<-c("ID","id","gene","expression","SNP","genotype")
colnames(gene1)<-c("ID","id","gene","expression","SNP","genotype")

gene1$pair<-"copy 1"
gene2$pair<-"copy 2"


longgene<-rbind.data.frame(gene1,gene2,stringsAsFactors = FALSE)

ggplot(longgene,aes(x=as.factor(genotype),y=expression,fill=pair))+geom_boxplot()+theme_bw(25)+facet_wrap(~ID,scales = "free_y",labeller = labeller(name = label_wrap_gen(width = 45)))

ggplot(longgene,aes(x=as.factor(genotype),y=expression,fill=pair))+geom_violin()+theme_bw(25)+facet_wrap(~ID,scales = "free_y",labeller = labeller(name = label_wrap_gen(width = 45)))



## see effect size for group in figure 1


## shuffled vs non shuffled
rna_full_table<-read.table(file="Synchro_uniq_fulltable_withinfo_04.txt",header=TRUE)
rna_full_table$snpid<-paste(ifelse(rna_full_table$snpschrom<9,paste("ssa0",rna_full_table$snpschrom,sep=""),paste("ssa",rna_full_table$snpschrom,sep="")),rna_full_table$snpstart,sep="_")


rna_full_table<-rna_full_table[rna_full_table$condition=="LL",]


dupstatus_list <- c("inter_chr_non-ohnolog","intra_chr_syntenic","intra_chr_non-syntenic","ohnolog")

rna_dupstatus<-rna_full_table[rna_full_table$dupstatus%in%dupstatus_list,] ##2210 eQTL-gene association not in dupstatus list


effectszie<-read.table(file="LL_cistrans_slope_Aug2023.txt",header=TRUE)

effectsize_rna<-merge(effectszie,rna_dupstatus,by.x = c("gene","SNP"),by.y=c("gene","snpid"))

library(ggplot2)
ggplot(effectsize_rna,aes(x=abs(slope),fill=dupstatus))+geom_histogram(color="black")+theme_bw(20)+facet_wrap(~dupstatus,scales = 'free_y')

ks.test(abs(effectsize_rna[effectsize_rna$dupstatus=="ohnolog","slope"]),abs(effectsize_rna[effectsize_rna$dupstatus=="intra_chr_syntenic","slope"]))
ks.test(abs(effectsize_rna[effectsize_rna$dupstatus=="ohnolog","slope"]),abs(effectsize_rna[effectsize_rna$dupstatus=="intra_chr_non-syntenic","slope"]))
ks.test(abs(effectsize_rna[effectsize_rna$dupstatus=="ohnolog","slope"]),abs(effectsize_rna[effectsize_rna$dupstatus=="inter_chr_non-ohnolog","slope"]))
ks.test(abs(effectsize_rna[effectsize_rna$dupstatus=="intra_chr_non-syntenic","slope"]),abs(effectsize_rna[effectsize_rna$dupstatus=="inter_chr_non-ohnolog","slope"]))
ks.test(abs(effectsize_rna[effectsize_rna$dupstatus=="intra_chr_syntenic","slope"]),abs(effectsize_rna[effectsize_rna$dupstatus=="inter_chr_non-ohnolog","slope"]))
ks.test(abs(effectsize_rna[effectsize_rna$dupstatus=="intra_chr_syntenic","slope"]),abs(effectsize_rna[effectsize_rna$dupstatus=="intra_chr_non-syntenic","slope"]))



### get both regulation

geneboth<-genepair[genepair$regulation=="bothdiffreg",]

genechr<-read.table(file="gene_SNPchr_num.txt",header=TRUE)

genechrsubset<-genechr[genechr$Var1%in%geneboth$gene1 | genechr$Var1%in%geneboth$gene2,]

write.table(genechrsubset,file="genediffeqtl.txt",row.names = FALSE)


##add chr information of genes
genepair<-read.table(file="genetabgenepair02.txt",header=TRUE)

tablegene<-read.table(file="https://salmobase.org/datafiles/TSV/genes/AtlanticSalmon/Ssal_v3.1/Ensembl_genes.tsv",header=TRUE,sep="\t")
tablegene<-tablegene[!is.na(as.numeric(as.character(tablegene$seqname))),]
tablegene$CHROM<-paste(ifelse(as.numeric(as.character(tablegene$seqname))>9,"ssa","ssa0"),tablegene$seqname,sep="")


for( i in 1:length(genepair$ID)){
  gene1<-genepair[i,"gene1"]
  gene2<-genepair[i,"gene2"]
  genepair[i,"gene1chr"]<-ifelse(length(tablegene[tablegene$gene_id==gene1,"CHROM"])==0,"NA",tablegene[tablegene$gene_id==gene1,"CHROM"])
  genepair[i,"gene2chr"]<-ifelse(length(tablegene[tablegene$gene_id==gene2,"CHROM"])==0,"NA",tablegene[tablegene$gene_id==gene2,"CHROM"])
}

write.table(genepair,file="genepair_withchr01.txt",row.names = FALSE)





## part 2 with new data
data<-read.table(file="synchroLL_finemap.txt",header=TRUE)

all_data<-read.table(file="rna_LL_all_august_version_01.txt",header=TRUE)

all_cis<-all_data[all_data$kind=="cis",c(1,2,5)]
all_cis$posterior_prob<-NA

colnames(all_cis)[2]<-"SNP"
data$kind<-"trans"

cistransdata<-rbind.data.frame(all_cis,data,stringsAsFactors = FALSE)


#cistransdata[cistransdata$gene=="ENSSSAG00000011938",]
#cistransdata[cistransdata$gene=="ENSSSAG00000108018",]

#cistransdata[cistransdata$gene=="ENSSSAG00000113858",]
#cistransdata[cistransdata$gene=="ENSSSAG00000120807",]


##function to find pair

genetab<-read.table(file="ss4r_dups_and_singletons_ENSrapid_convPipeline.tsv",header=TRUE)
genetab<-genetab[genetab$type=="ss4r",]
genetab<-genetab[,-c(1,2,6)]


getPairstatus<-function(threshold){
  data_subset<- cistransdata[is.na(cistransdata$posterior_prob) | cistransdata$posterior_prob >= threshold, ]
  for(y in 1:length(genetab$ID)){
    
    gene1row1<-genetab[y,"gene1"]
    tabpair1<-data_subset[data_subset$gene==gene1row1,]
    gene2row1<-genetab[y,"gene2"]
    tabpair2<-data_subset[data_subset$gene==gene2row1,]  
    
    ltab1<-length(tabpair1$kind)
    ltab2<-length(tabpair2$kind)
    
    reg<-"No eQTL"
    bothdiffreg<-"nobothdiffreg"
    bothsamereg<-"nobothsamereg"
    if(ltab1==0 && ltab2==0){
    }
    if(ltab1==0 && ltab2!=0){
      reg<-"Single copy eQTL"
    }
    if(ltab1!=0 && ltab2==0){
      reg<-"Single copy eQTL"
    }
    if(ltab1!=0 && ltab2!=0){
      issamesnp<-ifelse(tabpair1$SNP%in%tabpair2$SNP,"yes","no")
      if("no"%in%issamesnp){
        bothdiffreg<-"Distinct eQTL"
      }
      if("yes"%in%issamesnp){
        bothsamereg<-"Shared eQTL"
      }
    }
    
    
    
    ##say that if regulated by same snp then we keep only this result
    if(bothdiffreg=="Distinct eQTL" && bothsamereg=="Shared eQTL"){
      reg<-"Shared eQTL"
    }else if(bothdiffreg=="nobothdiffreg" && bothsamereg=="Shared eQTL"){
      reg<-"Shared eQTL"
    }else if(bothdiffreg=="Distinct eQTL" && bothsamereg=="nobothsamereg"){
      reg<-"Distinct eQTL"
    }
    
    genetab[y,"regulation"]<-reg
  }
  return(genetab$regulation)
}

for(y in seq(from=0,to=1,by=0.05)){
  result<-getPairstatus(y)
  namecol<-paste("regulation",y,sep="_")
  genetab[,namecol]<-result
}

write.table(genetab,file="genetabgenepair03.txt",row.names=FALSE)

##check how count change
library(dplyr)

# Extract regulation column names
regulation_columns <- grep("^regulation_", colnames(genetab), value = TRUE)

# Extract unique eQTL types
eqtl_types <- unique(genetab$regulation_0)

# Initialize an empty data frame to store results
summary_table <- data.frame(regulation_number = numeric(0))

# Add columns for each eQTL type
for (eqtl_type in eqtl_types) {
  summary_table[[eqtl_type]] <- integer(0)
}

# Loop through regulation columns
for (reg_col in regulation_columns) {
  regulation_number <- as.numeric(gsub("^regulation_", "", reg_col))
  row_data <- list(regulation_number)
  
  for (eqtl_type in eqtl_types) {
    nb_eqtl_type <- sum(genetab[[reg_col]] == eqtl_type, na.rm = TRUE)
    row_data[[eqtl_type]] <- nb_eqtl_type
  }
  
  # Add a row to the summary table
  summary_table <- rbind(summary_table, do.call(cbind, row_data))
}


library(dplyr)
library(tidyr)


# Reshape the summary_table to the desired format
summary_table_long <- summary_table %>%
  pivot_longer(cols = -V1, names_to = "eQTL_type", values_to = "count")

# Rename the eQTL types to be more human-readable
eqtl_type_names <- c("No eQTL", "Distinct eQTL", "Single copy eQTL", "Shared eQTL")
summary_table_long$eQTL_type <- factor(summary_table_long$eQTL_type, levels = eqtl_type_names)
colnames(summary_table_long)[1]<-"threshold"

summary_table_long$eQTL_type<-as.factor(summary_table_long$eQTL_type)
summary_table_long$eQTL_type<-factor(summary_table_long$eQTL_type,levels=c("Distinct eQTL","No eQTL","Shared eQTL", "Single copy eQTL"))


library(ggplot2)

ggplot(summary_table_long,aes(x=threshold,y=count,color=eQTL_type))+geom_line(size=2)+geom_point(size=4)+scale_color_manual(values=c("#F9B27C","#99B8EA","#8CE86D","#D186D8"))+theme_bw(20)



##analysis
library(ggplot2)

gene_pair_tab<-read.table(file="expression_gene_pair02.txt",header=TRUE)

##coefficient

# calculate spearman coefficient for each condition
result <- gene_pair_tab %>% 
  group_by(ID) %>% 
  summarize(spearman_coeff = cor(Expression_gene1, Expression_gene2, method = "spearman"))

genepair<-read.table(file="genetabgenepair03.txt",header=TRUE)

results2<-left_join(result, genepair, by = c("ID"="ID"))

df2 <- results2 %>%
  group_by(regulation_0) %>%
  summarise(Mean = mean(spearman_coeff))

ggplot(results2,aes(x=spearman_coeff,fill=regulation_0))+geom_histogram(color="black")+theme_bw(20)+facet_wrap(~regulation_0,scales = 'free_y')+xlab("Spearman coefficient")+scale_fill_manual(values=c("#F9B27C","#99B8EA","#8CE86D","#D186D8"))+geom_vline(data = df2, mapping = aes(xintercept = Mean),color="firebrick3",size=1.5)

#results2<-as.data.frame(results2)
#results2[results2$ID=="ENSSSAG00000011938_Ssal_ENSSSAG00000108018_Ssal",]
#results2[results2$ID=="ENSSSAG00000113858_Ssal_ENSSSAG00000120807_Ssal",]

##make expression boxplot
library(tidyverse)
Genotype <- read_table2("/mnt/project/MSLab/Marie/2022_Synchrosmolt/QTLTools/imputed_0.01_2730_smolts.traw")   %>% select(-c(1,3,4,5,6))
Expression <- read_table2("Synchro_mt_TPM.bed")    %>% select(-c(1,2,3,5,6))
#Genotype <- read_table2("imputed_0.01_2730_smolts.traw")   %>% select(-c(1,3,4,5,6))
Covariate <-  read_table2("covariate_synchro.txt")


negativedata<-results2[results2$spearman_coeff<0 & results2$regulation_0=="Shared eQTL",]

negativedatasubset<-negativedata[,c(1,3,4)]
negativedatasubset[9,]<-negativedatasubset[1,]
negativedatasubset[10,]<-negativedatasubset[2,]
negativedatasubset[11,]<-negativedatasubset[3,]
negativedatasubset[12:13,]<-negativedatasubset[4,]
negativedatasubset[14:15,]<-negativedatasubset[5,]
negativedatasubset[16,]<-negativedatasubset[6,]
negativedatasubset[17:20,]<-negativedatasubset[7,]
negativedatasubset[21:22,]<-negativedatasubset[8,]

##add snp info
neggene<-cistransdata[cistransdata$gene%in%negativedata$gene1 | cistransdata$gene%in%negativedata$gene2,]
neggene <- neggene  %>% arrange(gene)


negativedatasubset$snpid<-c("ssa14_51354372","ssa14_51354372","ssa14_51354372","ssa17_62041220","ssa02_14694003","ssa14_51354372","ssa07_56591350","ssa17_62354384","ssa14_51370830","ssa14_51370830","ssa14_51370830","ssa17_62047452","ssa17_62047483","ssa02_14695921","ssa02_14696160","ssa14_51370830","ssa07_56596682","ssa07_56604908","ssa07_56592523","ssa07_56594397","ssa17_62354986","ssa17_62357348")

tableexpression<-Get_SNP_gene_info(negativedatasubset,feature="pair")

tableexpression$sumexpression<-tableexpression$Expression_gene1+tableexpression$Expression_gene2

gene1<-tableexpression[,c(1,2,3,6,8,9)]
gene2<-tableexpression[,c(1,2,4,7,8,9)]

sumgene<-data.frame(ID=tableexpression$ID,id=tableexpression$id,gene=tableexpression$ID,expression=tableexpression$sumexpression,SNP=tableexpression$SNP,genotype=tableexpression$genotype)
colnames(sumgene)<-c("ID","id","gene","expression","SNP","genotype")
sumgene$pair<-"Sum of copies"

colnames(gene2)<-c("ID","id","gene","expression","SNP","genotype")
colnames(gene1)<-c("ID","id","gene","expression","SNP","genotype")

gene1$pair<-"copy 1"
gene2$pair<-"copy 2"


longgene<-rbind.data.frame(gene1,gene2,stringsAsFactors = FALSE)

ggplot(longgene,aes(x=as.factor(genotype),y=expression,fill=pair))+geom_boxplot()+theme_bw(25)+facet_wrap(~ID,scales = "free_y",labeller = labeller(name = label_wrap_gen(width = 45)))


##see effect size




effectszie<-read.table(file="LL_cistrans_slope_Aug2023.txt",header=TRUE)

genepair<-read.table(file="genetabgenepair03.txt",header=TRUE)
sharedtab<-genepair[genepair$regulation_0=="Shared eQTL",]

data<-read.table(file="synchroLL_finemap.txt",header=TRUE)

all_data<-read.table(file="rna_LL_all_august_version_01.txt",header=TRUE)

all_cis<-all_data[all_data$kind=="cis",c(1,2,5)]
all_cis$posterior_prob<-NA

colnames(all_cis)[2]<-"SNP"
data$kind<-"trans"

cistransdata<-rbind.data.frame(all_cis,data,stringsAsFactors = FALSE)


#sharedtabsnp<-GetOneeQTL(sharedtab)
#sharedtabsnp<-sharedtabsnp[sharedtabsnp$snp_gene1!="No Data",]



effectsziegene1<-merge(sharedtab,effectszie,by.x="gene1",by.y="gene")
effectsziegene2<-merge(sharedtab,effectszie,by.x="gene2",by.y="gene")

#effectsizegene1<-effectsziegene1[effectsziegene1$snpid==effectsziegene1$SNP,]
#effectsizegene2<-effectsziegene2[effectsziegene2$snpid==effectsziegene2$SNP,]

effectsizegene1<-effectsziegene1[effectsziegene1$kind=="trans",]
effectsizegene2<-effectsziegene2[effectsziegene2$kind=="trans",]


effectsizegene1$copy<-"copy1"
effectsizegene2$copy<-"copy2"

effectszieshared<-merge(effectsizegene1,effectsizegene2,by="ID")
effectszieshared$snpchrom<-sapply(strsplit(effectszieshared$SNP.x,"_"),FUN = `[[`, 1)
effectszieshared$negativecorr<-"no"
effectszieshared[effectszieshared$ID%in%negativedata$ID,"negativecorr"]<-"yes"


library(ggplot2)

ggplot(effectszieshared,aes(x=slope.x,y=slope.y,color=negativecorr))+geom_point(size=2)+theme_bw(20)

ggplot(effectszieshared,aes(x=slope.x,y=slope.y,color=negativecorr))+geom_point(size=2)+theme_bw(20)+geom_smooth(method=lm, color='#2C3E50')+ geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")

ggplot(effectszieshared,aes(x=abs(slope.x),y=abs(slope.y),color=negativecorr))+geom_point(size=2)+theme_bw(20)

ggplot(effectszieshared,aes(x=abs(slope.x),y=abs(slope.y),color=negativecorr))+geom_point(size=2)+theme_bw(20)+geom_smooth(method=lm, color='#2C3E50')+ geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")




##function to get only 1 eQTL in case of mutiple same eQTL (shared eQTL)
##data : a table with column gene1 and gene2 with ID of genes


GetOneeQTL<-function(table){
  for(i in 1:length(table$gene1)){
    ##get gene ID
    gene1<-table[i,"gene1"]
    gene2<-table[i,"gene2"]
    
    linkgene1<-paste("/net/fs-2/scale/OrionStore/Projects/MSLab/Marie/2022_Synchrosmolt/QTLTools/finemap/result/zscore/",gene1,".zscore.txt.cs",sep="")
    linkgene2<-paste("/net/fs-2/scale/OrionStore/Projects/MSLab/Marie/2022_Synchrosmolt/QTLTools/finemap/result/zscore/",gene2,".zscore.txt.cs",sep="")
    
    ##import score if it exists
    if(file.exists(linkgene1) && file.exists(linkgene2)){
      
      pairscoregene1<-read.table(file=linkgene1,header=TRUE)
      pairscoregene2<-read.table(file=linkgene2,header=TRUE)
      
      if(length(pairscoregene1$cs)>0 && length(pairscoregene2$cs)>0){
        ##format snp id
        pairscore_gene1<-strsplit(pairscoregene1$cs,"/")
        pairscore_gene2<-strsplit(pairscoregene2$cs,"/")
        
        ##for each sublist of snp for gene1
        for ( y in 1:length(pairscore_gene1)){
          sublist<-pairscore_gene1[y]
          ##check if the same list is in snp for gene2 in each sublist
          for (z in 1:length(pairscore_gene2)){
            sublist2<-pairscore_gene2[z]
            if(identical(sublist,sublist2)){
              sampled_snp <- sample(sublist, size = 1, replace = TRUE)
              pairscore_gene1[y]<-sampled_snp
              pairscore_gene2[z]<-sampled_snp
              break
            }else{
              
            }
          }
        }
        
        ##register snp info
        table[i,"snp_gene1"]<- paste(sapply(pairscore_gene1, paste, collapse =";"), collapse = "/")
        table[i,"snp_gene2"]<- paste(sapply(pairscore_gene2, paste, collapse =";"), collapse = "/")
      }else{
        table[i,"snp_gene1"]<-"No Data"
        table[i,"snp_gene2"]<-"No Data"
      }
    }else{
      table[i,"snp_gene1"]<-"file not found"
      table[i,"snp_gene2"]<-"file not found"
      
    }
  }
  return(table)
}

