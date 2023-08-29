
## get genes snp impact
tabSNP<-read.table(file="/mnt/users/cedi/files/SNP_data/tablesnpeff_SNP_Farmed_European_02.txt",header=TRUE)
tabSNP<-tabSNP[!duplicated(tabSNP),]

#sexdiffgenes<-c("ENSSSAG00000006666","ENSSSAG00000069409","ENSSSAG00000075378","ENSSSAG00000079572","ENSSSAG00000065844","ENSSSAG00000005963",
#"ENSSSAG00000091484","ENSSSAG00000052374","ENSSSAG00000079572","ENSSSAG00000112370","ENSSSAG00000097001","ENSSSAG00000005963")



tablegene<-read.table(file="https://salmobase.org/datafiles/TSV/genes/AtlanticSalmon/Ssal_v3.1/Ensembl_genes.tsv",header=TRUE,sep="\t")
tablegene<-tablegene[!is.na(as.numeric(as.character(tablegene$seqname))),]
tablegene$CHROM<-paste(ifelse(as.numeric(as.character(tablegene$seqname))>9,"ssa","ssa0"),tablegene$seqname,sep="")

#write.table(tablegene,file="/mnt/SCRATCH/cedi/phDSalmon/rna_expression/fusion_transrcipt/tablegene_ssaname.tab",row.names=FALSE)

#sexdifgenetab<-tablegene[tablegene$gene_id%in%sexdiffgenes,]



##use Marie's table
#data<-read.table(file="geneLL_class.csv",sep=",",header=TRUE)
data<-read.csv(file="geneLL_class.csv",sep=",",header=TRUE)
data[is.na(data$SNP),"SNP"]<-"none"
data$SNPstatus<-ifelse(data$SNP=="none","no","yes")

#test<-data[data$SNP!="none",]
#ggplot(test[test$cistrans=="cis",],aes(x=as.factor(expression_bin)))+geom_histogram(stat="count")

colnames(tabSNP)[5]<-"gene"

#test<-data[data$SNP!="none",]
#test2<-test[test$cistrans=="cis",] #155 data
#test2<-test[test$cistrans=="both",] #4370 data

###merge and count the number of SNPS affecting each gene
merged<-merge(data,tabSNP,by="gene")
merged$gene_length<-merged$Gene.end..bp.-merged$Gene.start..bp.

write.table(merged,file="mergedSNPeFFxgeneLLclass.txt",row.names=FALSE)

merged<-read.table(file="mergedSNPeFFxgeneLLclass.txt",header=TRUE)

##if i want to get gene with no impact snp
data$snpimpact<-ifelse(data$gene%in%merged$gene,"yes","no")
extractnosnp<-data[data$snpimpact=="no",]

library(dplyr)

dfcount <- merged %>% group_by(gene,expression_bin,gene_length,SNPstatus,cistrans) %>% 
  summarise(total_count=n(),.groups = 'drop') %>%
  as.data.frame()
dfcount$normalizedcount<-dfcount$total_count/dfcount$gene_length

library(ggplot2)

ggplot(dfcount,aes(x=as.factor(expression_bin),y=normalizedcount))+geom_boxplot()+ylim(0,0.05)+theme_bw(base_size = 20)
ggplot(dfcount,aes(x=SNPstatus,y=normalizedcount))+geom_boxplot()+ylim(0,0.05)+theme_bw(base_size = 20)

ggplot(dfcount,aes(x=cistrans,y=normalizedcount))+geom_boxplot()+ylim(0,0.05)+theme_bw(base_size = 20)

ggplot(dfcount,aes(x=as.factor(expression_bin),y=normalizedcount,fill=SNPstatus))+geom_boxplot()+ylim(0,0.05)+theme_bw(base_size = 20)

ggplot(dfcount,aes(x=as.factor(expression_bin),y=normalizedcount,fill=cistrans))+geom_boxplot()+ylim(0,0.05)+theme_bw(base_size = 20)

##if we look at impact now
dfcount <- merged %>% group_by(gene,expression_bin,gene_length,SNPstatus,impact) %>% 
  summarise(total_count=n(),.groups = 'drop') %>%
  as.data.frame()

df <- dfcount %>%
  group_by(gene) %>%
  top_n(1, desc(factor(impact, levels=c("HIGH","MODERATE","LOW","MODIFIER")))) %>%
  ungroup()

df2 <- df %>% group_by(impact,expression_bin) %>% 
  summarise(total_count=n(),.groups = 'drop') %>%
  as.data.frame()

ggplot(df2,aes(x=as.factor(expression_bin),y=total_count,fill=impact))+geom_bar(stat="identity",position="fill")+theme_bw(base_size = 20)

df2 <- df %>% group_by(impact,SNPstatus) %>% 
  summarise(total_count=n(),.groups = 'drop') %>%
  as.data.frame()

ggplot(df2,aes(x=SNPstatus,y=total_count,fill=impact))+geom_bar(stat="identity",position="fill")+theme_bw(base_size = 20)


##facet with cistrans

dfcount <- merged %>% group_by(gene,expression_bin,gene_length,SNPstatus,impact,cistrans) %>% 
  summarise(total_count=n(),.groups = 'drop') %>%
  as.data.frame()

df <- dfcount %>%
  group_by(gene) %>%
  top_n(1, desc(factor(impact, levels=c("HIGH","MODERATE","LOW","MODIFIER")))) %>%
  ungroup()

df2 <- df %>% group_by(impact,cistrans) %>% 
  summarise(total_count=n(),.groups = 'drop') %>%
  as.data.frame()

ggplot(df2,aes(x=cistrans,y=total_count,fill=impact))+geom_bar(stat="identity",position="fill")+theme_bw(base_size = 20)


df2 <- df %>% group_by(impact,expression_bin,cistrans) %>% 
  summarise(total_count=n(),.groups = 'drop') %>%
  as.data.frame()

ggplot(df2,aes(x=as.factor(expression_bin),y=total_count,fill=impact))+geom_bar(stat="identity",position="fill")+theme_bw(base_size = 20)+facet_wrap(~cistrans)

##compute ratio


df2[is.na(df2$cistrans),"cistrans"]<-"none"
result <- aggregate(total_count ~ expression_bin + cistrans, data = df2, sum)
colnames(result)[3]<-"sum_total_count"

result2 <- df2 %>% 
  left_join(result , by = c("cistrans", "expression_bin")) %>% 
  mutate(proportion = total_count / sum_total_count)

ggplot(result2[result2$impact=="HIGH",],aes(x=as.factor(expression_bin),y=proportion,fill=cistrans))+geom_bar(stat="identity",position="dodge")+theme_bw(base_size = 20)


df2 <- df %>% group_by(impact,expression_bin,SNPstatus) %>% 
  summarise(total_count=n(),.groups = 'drop') %>%
  as.data.frame()

result <- aggregate(total_count ~ expression_bin + SNPstatus, data = df2, sum)
colnames(result)[3]<-"sum_total_count"

result2 <- df2 %>% 
  left_join(result , by = c("SNPstatus", "expression_bin")) %>% 
  mutate(proportion = total_count / sum_total_count)

ggplot(result2,aes(x=as.factor(expression_bin),y=proportion,fill=SNPstatus))+geom_bar(stat="identity",position="dodge")+theme_bw(base_size = 20)+facet_wrap(~impact)



##a try without the modifier....

merged2<-merged[merged$impact!="MODIFIER",]

dfcount <- merged2 %>% group_by(gene,expression_bin,gene_length,SNPstatus,impact,cistrans) %>% 
  summarise(total_count=n(),.groups = 'drop') %>%
  as.data.frame()

df <- dfcount %>%
  group_by(gene) %>%
  summarize(
    expression_bin = first(expression_bin),
    gene_length = first(gene_length),
    SNPstatus = first(SNPstatus),
    cistrans = first(cistrans),
    score = sum(case_when(
      impact == "HIGH" ~ 3,
      impact == "MODERATE" ~ 2,
      impact == "LOW" ~ 1,
      TRUE ~ 0
    ) * total_count)
  ) %>%
  ungroup()



ggplot(df,aes(x=as.factor(expression_bin),y=score))+geom_boxplot()+theme_bw(base_size = 20)+ylim(0,30)

ggplot(df,aes(x=SNPstatus,y=score))+geom_boxplot()+theme_bw(base_size = 20)+ylim(0,30)

ggplot(df,aes(x=cistrans,y=score))+geom_boxplot()+theme_bw(base_size = 20)+ylim(0,30)


ggplot(df,aes(x=as.factor(expression_bin),y=score,fill=cistrans))+geom_boxplot()+theme_bw(base_size = 20)+ylim(0,30)


## now with regions

tabSNP<-read.table(file="/mnt/users/cedi/files/SNP_data/tablesnpeff_Farmed_European_withregion_02.txt",header=TRUE)
tabSNP<-tabSNP[!duplicated(tabSNP),]

data<-read.csv(file="geneLL_class.csv",sep=",",header=TRUE)
data[is.na(data$SNP),"SNP"]<-"none"
data$SNPstatus<-ifelse(data$SNP=="none","no","yes")
colnames(tabSNP)[5]<-"gene"

###merge and count the number of SNPS affecting each gene
merged<-merge(data,tabSNP,by="gene")
merged$gene_length<-merged$Gene.end..bp.-merged$Gene.start..bp.

write.table(merged,file="mergedSNPeFFxgeneLLclass_withregions.txt",row.names=FALSE)

merged<-read.table(file="mergedSNPeFFxgeneLLclass_withregions.txt",header=TRUE)


library(dplyr)

dfcount <- merged %>% group_by(gene,expression_bin,gene_length,SNPstatus,cistrans,regions) %>% 
  summarise(total_count=n(),.groups = 'drop') %>%
  as.data.frame()
dfcount$normalizedcount<-dfcount$total_count/dfcount$gene_length

library(ggplot2)

df2<- dfcount %>% group_by(expression_bin,regions) %>% summarise(sum_total_count=sum(total_count))

ggplot(df2,aes(x=as.factor(expression_bin),y=sum_total_count,fill=regions))+geom_bar(stat="identity",position="fill")+theme_bw(base_size = 20)


df2<- dfcount %>% group_by(cistrans,regions) %>% summarise(sum_total_count=sum(total_count))

ggplot(df2,aes(x=cistrans,y=sum_total_count,fill=regions))+geom_bar(stat="identity",position="fill")+theme_bw(base_size = 20)

###
merged2<-merged %>% mutate(regions2 = case_when(
  regions == "Downstream" | regions == "Upstream" ~ "Up_Down_stream",
  regions == "UTR_3_Prime" | regions == "UTR_5_Prime" ~ "UTR",
  regions == "Start codon" | regions == "Stop codon" ~ "Start_Stop_codon",
  TRUE ~ regions # Retain the original value for all other rows
))

dfcount <- merged2 %>% group_by(gene,expression_bin,gene_length,SNPstatus,cistrans,regions2) %>% 
  summarise(total_count=n(),.groups = 'drop') %>%
  as.data.frame()
dfcount$normalizedcount<-dfcount$total_count/dfcount$gene_length

df2<- dfcount %>% group_by(expression_bin,regions2) %>% summarise(sum_total_count=sum(total_count))

ggplot(df2,aes(x=as.factor(expression_bin),y=sum_total_count,fill=regions2))+geom_bar(stat="identity",position="fill")+theme_bw(base_size = 20)

ggplot(df2,aes(x=as.factor(expression_bin),y=sum_total_count,fill=regions2))+geom_bar(stat="identity")+theme_bw(base_size = 20)+facet_wrap(~regions2,scales = 'free_y')

df3<-df2 %>% group_by(expression_bin) %>%
  mutate(proportion = sum_total_count / sum(sum_total_count))

ggplot(df3,aes(x=as.factor(expression_bin),y=proportion,fill=regions2))+geom_bar(stat="identity")+theme_bw(base_size = 20)+facet_wrap(~regions2,scales = 'free_y')

df2<- dfcount %>% group_by(expression_bin,regions2,SNPstatus) %>% summarise(sum_total_count=sum(total_count))

df3<-df2 %>% group_by(expression_bin) %>%
  mutate(proportion = sum_total_count / sum(sum_total_count))

ggplot(df3,aes(x=as.factor(expression_bin),y=proportion,fill=SNPstatus))+geom_bar(stat="identity")+theme_bw(base_size = 20)+facet_wrap(~regions2,scales = 'free_y')



##get SNP that cis affect gene. get those who only affect one gene, those who affect several cis gene and those who affect several cis and trans genes

##but first how many singleton vs LOR vs AOR we expect

genetab<-read.table(file="ss4r_dups_and_singletons_ENSrapid_convPipeline.tsv",header=TRUE)

genetabdup<-genetab[genetab$type=="ss4r" & !is.na(genetab$redip.class),]
genetabsing<-genetab[genetab$type=="singleton",]
#total: 9770 + 7879 = 17649
# 55% are singleton / 33% AORe / 11% LORe
#test<-genetabdup[genetabdup$redip.class=="LORe",]

rna_full_table<-read.table(file="Synchro_uniq_fulltable_withinfo_04.txt",header=TRUE)

rna_full_table_LL<-rna_full_table[rna_full_table$condition=="LL",]
rna_full_table_LL$snpcoor<-paste(rna_full_table_LL$snpstart,rna_full_table_LL$snpschrom,sep="_")
#length(unique(rna_full_table_LL$snpcoor))
##now get cis eQTL and trans eQTL
datasnp<-c()
for( i in unique(rna_full_table_LL$snpcoor)){
  subtab<-rna_full_table_LL[rna_full_table_LL$snpcoor==i,]
  nbgene<-length(subtab$cistrans)
  ##cis scenarios
  if("cis"%in%subtab$cistrans){
    if(nbgene==1){
      cisstatus<-"one cis"
    }else{
      if("trans"%in%subtab$cistrans){
        cisstatus<-"cis and trans"
      }else{
        cisstatus<-"multiple cis"
      }
    }
  }else{
    cisstatus<-"no cis"
  }
  ##trans scenarios
  if("trans"%in%subtab$cistrans){
    transcount<-nbgene
  }else{
    transcount<-0
  }
  vect<-c(i,cisstatus,transcount)
  datasnp<-rbind.data.frame(datasnp,vect,stringsAsFactors = FALSE)
}
colnames(datasnp)<-c("snpcoord","cisstatus","transcount")
datasnp$POS<-sapply(strsplit(datasnp$snpcoord,"_"),FUN = `[[`, 1)
vectchr<-as.numeric(as.character(sapply(strsplit(datasnp$snpcoord,"_"),FUN = `[[`, 2)))
datasnp$CHROM<-ifelse(vectchr<10,paste("ssa0",vectchr,sep=""),paste("ssa",vectchr,sep=""))

#length(datasnp[datasnp$cisstatus=="cis and trans","cisstatus"])##620 are both
#length(datasnp[datasnp$cisstatus!="cis and trans" & datasnp$transcount>0,"cisstatus"]) ##5549 trans
#length(datasnp[datasnp$cisstatus=="one cis" | datasnp$cisstatus=="multiple cis" ,"cisstatus"]) ##5887 cis

FE_SNP<-read.table(file="/mnt/SCRATCH/cedi/phDSalmon/sv_detection/filtervcf/SNP_info_table_Farmed_Europe.txt",header=TRUE)
colnames(FE_SNP)<-c("CHROM","POS","END","ID","REF","ALT","QUAL","FREQ","POPULATION","VARIANT_TYPE")

mergedsnp<-merge(datasnp,FE_SNP,by=c("CHROM","POS"))
mergedsnp$MAF<-ifelse(mergedsnp$FREQ<0.5,mergedsnp$FREQ,1-mergedsnp$FREQ)

write.table(mergedsnp,file="snpLL_popinfo01.txt",row.names=FALSE)

##pop analysis
popsnp<-read.table(file="snpLL_popinfo01.txt",header=TRUE)

library(ggplot2)

#ggplot(popsnp,aes(x=cisstatus,y=MAF,fill=cisstatus))+geom_boxplot()
ggplot(popsnp[popsnp$transcount>0,],aes(x=transcount,y=MAF))+geom_point(alpha=0.1)+scale_x_continuous(trans='log10')

popsnp$bin<-ifelse(popsnp$MAF<0.1,0.1,ifelse(popsnp$MAF<0.2 & popsnp$MAF>=0.1,0.2,ifelse(popsnp$MAF<0.3 & popsnp$MAF>=0.2,0.3,ifelse(popsnp$MAF<0.4 & popsnp$MAF>=0.3,0.4,0.5))))

popsnpno0<-popsnp[popsnp$transcount>0,]

popsnpno0$countbin<-ifelse(popsnpno0$transcount==1,"1",ifelse(popsnpno0$transcount==2,"2",">2"))
popsnpno0$countbin<-as.factor(popsnpno0$countbin)
popsnpno0$countbin<-factor(popsnpno0$countbin,levels=c("1","2",">2"))


ggplot(popsnp[popsnp$transcount>0,],aes(x=transcount,fill=as.factor(bin)))+geom_histogram(color="black")+facet_wrap(~bin)+scale_x_continuous(trans='log10')

ggplot(popsnp[popsnp$transcount>0,],aes(x=as.factor(bin),y=transcount,fill=as.factor(bin)))+geom_boxplot()+scale_y_continuous(trans='log10')

ggplot(popsnpno0,aes(x=as.factor(countbin),y=MAF,fill=as.factor(countbin)))+geom_boxplot()+xlab("Number of genes impacted")+theme_light()+ylim(0,0.7)+theme_bw(20)

ks.test(popsnpno0[popsnpno0$countbin=="1","MAF"],popsnpno0[popsnpno0$countbin=="2","MAF"])

ks.test(popsnpno0[popsnpno0$countbin=="more2","MAF"],popsnpno0[popsnpno0$countbin=="2","MAF"])

ks.test(popsnpno0[popsnpno0$countbin=="more2","MAF"],popsnpno0[popsnpno0$countbin=="1","MAF"])

library(tidyverse)
library(ggsignif)

ggplot(popsnpno0,aes(x=as.factor(countbin),y=MAF,fill=as.factor(countbin)))+geom_boxplot()+xlab("Number of genes impacted")+geom_signif(comparisons = list(c("2", "1"),
                                                                                                                                                           c("2", ">2"),
                                                                                                                                                           c(">2", "1")
),
test = "ks.test", step_increase = 0.075,
map_signif_level = TRUE, tip_length = 0)



library(tidyverse)

poptrans<-popsnp[popsnp$transcount>0,]

dfcount <- poptrans %>% group_by(transcount,bin) %>% 
  summarise(total_count=n(),.groups = 'drop') %>%
  as.data.frame()

ggplot(dfcount,aes(x=as.factor(transcount),y=total_count,fill=as.factor(bin)))+geom_bar(position="fill",stat="identity")

##

popsnponecis<-popsnp[popsnp$cisstatus=="one cis",]

rna_full_table<-read.table(file="Synchro_uniq_fulltable_withinfo_04.txt",header=TRUE)

rna_full_table_LL<-rna_full_table[rna_full_table$condition=="LL",]
rna_full_table_LL$snpcoord<-paste(rna_full_table_LL$snpstart,rna_full_table_LL$snpschrom,sep="_")


genetab<-read.table(file="ss4r_dups_and_singletons_ENSrapid_convPipeline.tsv",header=TRUE)

genetabdup<-genetab[genetab$type=="ss4r" & !is.na(genetab$redip.class),]
genetabsing<-genetab[genetab$type=="singleton",]

popsnponecismerged<-merge(popsnponecis,rna_full_table_LL,by=c("snpcoord"))

##find status of the gene
for( i in 1:length(popsnponecismerged$CHROM)){
  gene<-popsnponecismerged[i,"gene"]
  if(gene%in%genetabsing$gene1){
    popsnponecismerged[i,"gene status"]<-"singleton"
  }else if(gene%in%genetabdup$gene1){
    popsnponecismerged[i,"gene status"]<-genetabdup[genetabdup$gene1==gene,"redip.class"]
  }else if(gene%in%genetabdup$gene2){
    popsnponecismerged[i,"gene status"]<-genetabdup[genetabdup$gene2==gene,"redip.class"]
  }else{
    popsnponecismerged[i,"gene status"]<-"Unknown"
  }
}

#write.table(popsnponecismerged,file="cisgenestatus02.txt",row.names = FALSE)

popsnponecismerged<-read.table(file="cisgenestatus02.txt",header=TRUE)

ggplot(popsnponecismerged,aes(x=MAF,fill=`gene.status`))+geom_density(color="black")+facet_wrap(~`gene.status`)

ggplot(popsnponecismerged,aes(x=gene.status,y=MAF,fill=gene.status))+geom_boxplot()+xlab("gene status")+theme_light()+ylim(0,0.7)+theme_bw(20)

ggplot(popsnponecismerged,aes(x=gene.status,y=MAF,fill=gene.status))+geom_violin()

ks.test(popsnponecismerged[popsnponecismerged$gene.status=="AORe","MAF"],popsnponecismerged[popsnponecismerged$gene.status=="Unknown","MAF"])

ks.test(popsnponecismerged[popsnponecismerged$gene.status=="AORe","MAF"],popsnponecismerged[popsnponecismerged$gene.status=="singleton","MAF"])

ks.test(popsnponecismerged[popsnponecismerged$gene.status=="LORe","MAF"],popsnponecismerged[popsnponecismerged$gene.status=="singleton","MAF"])

ks.test(popsnponecismerged[popsnponecismerged$gene.status=="Unknown","MAF"],popsnponecismerged[popsnponecismerged$gene.status=="singleton","MAF"])

ks.test(popsnponecismerged[popsnponecismerged$gene.status=="LORe","MAF"],popsnponecismerged[popsnponecismerged$gene.status=="AORe","MAF"])

ks.test(popsnponecismerged[popsnponecismerged$gene.status=="LORe","MAF"],popsnponecismerged[popsnponecismerged$gene.status=="Unknown","MAF"])


ggplot(popsnponecismerged,aes(x=gene.status,y=MAF,fill=gene.status))+geom_boxplot()+xlab("gene status")+geom_signif(comparisons = list(c("AORe", "LORe"),
                                                                                                                                       c("AORe", "singleton"),
                                                                                                                                       c("AORe", "Unknown"), c("LORe","singleton"),c("singleton","Unknown"),c("LORe","Unknown")),test = "ks.test", step_increase = 0.075,map_signif_level = TRUE, tip_length = 0)

#tabcisexp<-popsnponecismerged[,c("CHROM","POS")]
#write.table(tabcisexp,file="tabcis.txt",row.names=FALSE,col.names =FALSE,quote=FALSE,sep="\t")

#tabtransexp<-popsnp[popsnp$transcount>0,c("CHROM","POS")]
#write.table(tabtransexp,file="tabtrans.txt",row.names=FALSE,col.names =FALSE,quote=FALSE,sep="\t")



##check fst
fstcis<-read.table(file="/mnt/SCRATCH/cedi/phDSalmon/rna_expression/WGD_geneexp/FE_vs_WE_singlecis.txt.weir.fst",header=TRUE)
colnames(fstcis)[3]<-"FST"
fstcismerged<-merge(popsnponecismerged,fstcis,by=c("CHROM","POS"))

ggplot(fstcismerged,aes(x=FST,fill=`gene.status`))+geom_histogram(color="black")+facet_wrap(~`gene.status`)

ggplot(fstcismerged,aes(x=gene.status,y=FST,fill=gene.status))+geom_boxplot()+scale_y_continuous(trans='log')+ylab("log(FST)")

ggplot(fstcismerged,aes(x=gene.status,y=FST,fill=gene.status))+geom_boxplot()+theme_light()+theme_bw(20)+xlab("gene status")

ks.test(fstcismerged[fstcismerged$gene.status=="AORe","FST"],fstcismerged[fstcismerged$gene.status=="Unknown","FST"])

ks.test(fstcismerged[fstcismerged$gene.status=="AORe","FST"],fstcismerged[fstcismerged$gene.status=="singleton","FST"])

ks.test(fstcismerged[fstcismerged$gene.status=="LORe","FST"],fstcismerged[fstcismerged$gene.status=="singleton","FST"])

ks.test(fstcismerged[fstcismerged$gene.status=="Unknown","FST"],fstcismerged[fstcismerged$gene.status=="singleton","FST"])

ks.test(fstcismerged[fstcismerged$gene.status=="LORe","FST"],fstcismerged[fstcismerged$gene.status=="Unknown","FST"])

ks.test(fstcismerged[fstcismerged$gene.status=="LORe","FST"],fstcismerged[fstcismerged$gene.status=="AORe","FST"])


fsttrans<-read.table(file="/mnt/SCRATCH/cedi/phDSalmon/rna_expression/WGD_geneexp/FE_vs_WE_trans.txt.weir.fst",header=TRUE)
colnames(fsttrans)[3]<-"FST"

transtable<-popsnp[popsnp$transcount>0,]
transtable<-popsnpno0
transtablemerged<-merge(transtable,fsttrans,by=c("CHROM","POS"))

ggplot(transtablemerged,aes(x=transcount,y=FST))+geom_point(alpha=0.5)+scale_x_continuous(trans='log10')

transtablemerged$color<-ifelse(transtablemerged$CHROM=="ssa14","orange","dodgerblue3")

ggplot(transtablemerged,aes(x=transcount,y=FST))+geom_point(alpha=0.5,color=transtablemerged$color)+scale_x_continuous(trans='log10')

ggplot(transtablemerged,aes(x=as.factor(countbin),y=FST,fill=as.factor(countbin)))+geom_boxplot()+theme_light()+xlab("Number of genes impacted")+theme_bw(20)

ks.test(transtablemerged[transtablemerged$countbin==1,"FST"],transtablemerged[transtablemerged$countbin==2,"FST"])

ks.test(transtablemerged[transtablemerged$countbin=="more2","FST"],transtablemerged[transtablemerged$countbin==2,"FST"])

ks.test(transtablemerged[transtablemerged$countbin=="more2","FST"],transtablemerged[transtablemerged$countbin==1,"FST"])

##add freq of wild

transtablemerged$mafFE<-transtablemerged$MAF
transtablemerged$freqFE<-transtablemerged$FREQ



WE_SNP<-read.table(file="/mnt/SCRATCH/cedi/phDSalmon/sv_detection/filtervcf/SNP_info_table_Wild_Europe.txt",header=TRUE)
colnames(WE_SNP)<-c("CHROM","POS","END","ID","REF","ALT","QUAL","FREQ","POPULATION","VARIANT_TYPE")

WE_SNP$MAF<-ifelse(WE_SNP$FREQ<0.5,WE_SNP$FREQ,1-WE_SNP$FREQ)

mergedFEWE<-merge(transtablemerged,WE_SNP,by=c("CHROM","POS"))

mergedFEWE$mafWE<-mergedFEWE$MAF.y
mergedFEWE$freqWE<-mergedFEWE$FREQ.y

transtabWEFE<-mergedFEWE[,-c(20:28)]

ggplot(transtabWEFE,aes(x=freqFE,y=freqWE))+geom_point(alpha=0.5,size=transtabWEFE$transcount/50,color=transtabWEFE$color)+theme_bw(20)

ggplot(transtabWEFE,aes(x=mafFE,y=mafWE))+geom_point(alpha=0.5,size=transtabWEFE$transcount/50,color=transtabWEFE$color)+theme_bw(20)

#write.table(transtabWEFE,file="Marie_poptransfile.txt",row.names=FALSE)

##add table to find if gene counterpart is in the table or not
genetab<-read.table(file="ss4r_dups_and_singletons_ENSrapid_convPipeline.tsv",header=TRUE)
genetabss4r<-genetab[genetab$type=="ss4r",]

data$genedup<-ifelse(data$gene%in%genetabss4r$gene1,genetabss4r$gene2,ifelse(data$gene%in%genetabss4r$gene2,genetabss4r$gene1,NA))
data$presdup<-ifelse(data$genedup%in%data$gene,"yes","no")

write.table(data,file="geneLL_class_genedup.csv",row.names=FALSE)



## get the genes that are TF
library(stringr)
tfgenes<-ifelse(str_detect(tablegene$description, "transcription "),tablegene$gene_id,NA)
tfgenes<-tfgenes[!is.na(tfgenes)]
genewithdesc<-tablegene[!is.na(tablegene$description),"gene_id"]

tfgenes2<-ifelse(str_detect(tablegene$description, "transcription"),tablegene$description,NA)
tfgenes2<-tfgenes2[!is.na(tfgenes2)]
tfgenes2<-as.data.frame(tfgenes2)

tfgenes<-ifelse(str_detect(tablegene$description, "transcription factor"),tablegene$gene_id,NA)
tfgenes<-tfgenes[!is.na(tfgenes)]
genewithdesc<-tablegene[!is.na(tablegene$description),"gene_id"]

##add a column to check if the SNP are in a gene, and if this gene is a transcription factor
rna_full_table<-read.table(file="Synchro_uniq_fulltable_withregulation_03.txt",header=TRUE)
rna_full_table$snpschrom<-ifelse(rna_full_table$snpschrom>9,paste("ssa",rna_full_table$snpschrom,sep=""),paste("ssa0",rna_full_table$snpschrom,sep=""))

for (i in 1:length(rna_full_table$gene)){
  snpchr<-rna_full_table[i,"snpschrom"]
  snpstart<-rna_full_table[i,"snpstart"]
  snpTF<-"nogene"
  for (y in 1:length(tablegene$seqname)){
    if(snpchr==tablegene[y,"CHROM"] && snpstart>=tablegene[y,"start"] && snpstart<=tablegene[y,"end"]){
      if(is.na(tablegene[y,"description"])){
        snpTF<-"gene"
      }else if(str_detect(tablegene[y,"description"], "transcription factor")){
        snpTF<-"TF"
        break
      }else{
        snpTF<-"gene"
      }
    }
  }
  rna_full_table[i,"snpTF"]<-snpTF
}

write.table(rna_full_table,file="Synchro_uniq_fulltable_withregulation_snpTF_01.txt")
#rna_full_table[rna_full_table$snpTT=="TF",]

#rna_full_table<-rna_full_table[!is.na(rna_full_table$snpTT),]




## do the same but considering that if the SNP is xkb around the TF then it is
library(stringr)
distancefromtf<-20000

tablegene$newstart<-tablegene$start-distancefromtf
tablegene$newend<-tablegene$end+distancefromtf

for (i in 1:length(rna_full_table$gene)){
  snpchr<-rna_full_table[i,"snpschrom"]
  snpstart<-rna_full_table[i,"snpstart"]
  tablegene$match<-ifelse(tablegene$newstart<=snpstart & tablegene$newend>=snpstart & tablegene$CHROM==snpchr,TRUE,FALSE)
  genestate<-ifelse(tablegene$match==TRUE & !is.na(tablegene$description) & str_detect(tablegene$description, "transcription factor"),"TF",ifelse(tablegene$match==FALSE,"nogene","gene"))
  if(length(genestate[genestate=="TF"])==0){
    if(length(genestate[genestate=="gene"])==0){
      rna_full_table[i,"snpTF"]<-"nogene"
    }else{
      rna_full_table[i,"snpTF"]<-"gene"
    }
  }else{
    rna_full_table[i,"snpTF"]<-"TF"
  }
}

write.table(rna_full_table,file="Synchro_uniq_fulltable_withregulation_snpTF_03.txt")

##non optimized
for (i in 1:length(rna_full_table$gene)){
  snpchr<-rna_full_table[i,"snpschrom"]
  snpstart<-rna_full_table[i,"snpstart"]
  snpTF<-"nogene"
  for (y in 1:length(tablegene$seqname)){
    if(snpchr==tablegene[y,"CHROM"] && snpstart>=tablegene[y,"newstart"] && snpstart<=tablegene[y,"newend"]){
      if(is.na(tablegene[y,"description"])){
        snpTF<-"gene"
      }else if(str_detect(tablegene[y,"description"], "transcription factor")){
        snpTF<-"TF"
        break
      }else{
        snpTF<-"gene"
      }
    }
  }
  rna_full_table[i,"snpTF"]<-snpTF
}

write.table(rna_full_table,file="Synchro_uniq_fulltable_withregulation_snpTF_02.txt")


##see if any xpehh from Pauline overlap with eQTL

xpehhtab<-read.table(file="/mnt/SCRATCH/pabu/bedtools/input2/xp_ehh_euro.bed")
colnames(xpehhtab)<-c("CHROM","POS","END")
xpehhtab$CHROM<-ifelse(xpehhtab$CHROM<10,paste("ssa0",xpehhtab$CHROM,sep=""),paste("ssa",xpehhtab$CHROM,sep=""))

rna_full_table<-read.table(file="Synchro_uniq_fulltable_withinfo_04.txt",header=TRUE)
rna_full_table<-rna_full_table[rna_full_table$condition=="LL",]
rna_full_table$CHROM<-ifelse(rna_full_table$snpschrom<10,paste("ssa0",rna_full_table$snpschrom,sep=""),paste("ssa",rna_full_table$snpschrom,sep=""))

options("scipen"=100, "digits"=4) ##disable scientific notation for later exportation
write.table(xpehhtab,file="xpehh_edit_pauline.bed",row.names = FALSE,col.names = FALSE,quote = FALSE,sep="\t")

rna_full_table_bed<-data.frame(rna_full_table$CHROM,rna_full_table$snpstart,rna_full_table$snpstart)
write.table(rna_full_table_bed,file="rna_eQTL.bed",row.names = FALSE,col.names = FALSE,quote = FALSE,sep="\t")


##results
tab_intersect<-read.table(file="/mnt/SCRATCH/cedi/phDSalmon/rna_expression/candidate_SNPs/intersect_eQTL_xpehh.txt")

tab_intersect$snpid<-paste(tab_intersect$V1,tab_intersect$V2,sep="_")
rna_full_table$snpid<-paste(rna_full_table$CHROM,rna_full_table$snpstart,sep="_")
subsetrna<-rna_full_table[rna_full_table$snpid%in%tab_intersect$snpid,]


##list genes
## cleavage and polyadenilation specific factor : https://www.ensembl.org/Salmo_salar/Gene/Summary?db=core;g=ENSSSAG00000041135;r=27:16940619-16963211
## choline phosphotransferase : https://www.ensembl.org/Salmo_salar/Gene/Summary?db=core;g=ENSSSAG00000089632;r=7:46103701-46155453
## makorin, ring finger protein : https://www.ensembl.org/Salmo_salar/Gene/Summary?db=core;g=ENSSSAG00000009919;r=17:77305030-77335674
## tuftelin interacting protein : https://www.ensembl.org/Salmo_salar/Gene/Summary?db=core;g=ENSSSAG00000074395;r=24:15578622-15587289;t=ENSSSAT00000134571
## myoneurin like gene: https://www.ensembl.org/Salmo_salar/Gene/Summary?db=core;g=ENSSSAG00000077916;r=29:27073639-27096456;t=ENSSSAT00000144957 / https://www.ncbi.nlm.nih.gov/gene/106590491
## unknown : https://www.ensembl.org/Salmo_salar/Gene/Summary?db=core;g=ENSSSAG00000003771;r=6:47321227-47327207 
## small subunit processome component: https://www.ensembl.org/Salmo_salar/Gene/Summary?db=core;g=ENSSSAG00000099529;r=27:6322783-6332626;t=ENSSSAT00000188575
## proteasome subunit : https://www.ensembl.org/Salmo_salar/Gene/Summary?db=core;g=ENSSSAG00000085574;r=27:10518241-10521084;t=ENSSSAT00000143990
## unknown : https://www.ensembl.org/Salmo_salar/Gene/Summary?db=core;g=ENSSSAG00000086391;r=27:12832006-12846200
## VPS52 subunit of GARP complex : https://www.ensembl.org/Salmo_salar/Gene/Summary?db=core;g=ENSSSAG00000078638;r=27:10789984-10840875 / https://www.ncbi.nlm.nih.gov/gene/6293
## unknown: https://www.ensembl.org/Salmo_salar/Gene/Summary?db=core;g=ENSSSAG00000009103;r=27:11417438-11456842
## creb regulated transcription coactivator : https://www.ensembl.org/Salmo_salar/Gene/Summary?db=core;g=ENSSSAG00000041719;r=27:12402729-12420179 / https://www.ncbi.nlm.nih.gov/gene/106588329
## unknown : https://www.ensembl.org/Salmo_salar/Gene/Summary?db=core;g=ENSSSAG00000077337;r=27:10386765-10391729
## lantibiotic synthetase component: https://www.ensembl.org/Salmo_salar/Gene/Summary?db=core;g=ENSSSAG00000078655;r=29:18443597-18510367
## dis3l2 exoribonuclease : https://www.ensembl.org/Salmo_salar/Gene/Summary?db=core;g=ENSSSAG00000004442;r=14:10851627-10867844 / https://www.ncbi.nlm.nih.gov/gene/129563
## sf3a1 : https://www.ensembl.org/Salmo_salar/Gene/Summary?db=core;g=ENSSSAG00000007482;r=20:21080066-21097383
## Unknown : https://www.ensembl.org/Salmo_salar/Gene/Summary?db=core;g=ENSSSAG00000102113;r=1:29134268-29137370;t=ENSSSAT00000250174
## mediator complex subunit 6 : https://www.ensembl.org/Salmo_salar/Gene/Summary?db=core;g=ENSSSAG00000000096;r=1:5108913-5125037



##make a bed for eQTL and SNPeff file

snpefftab<-read.table(file="/mnt/users/cedi/files/SNP_data/tablesnpeff_Farmed_European_withregion_02.txt",header=TRUE)
options("scipen"=100, "digits"=4) ##disable scientific notation for later exportation
bedfilesnpeff<-snpefftab[,1:3]

write.table(bedfilesnpeff,file="snpeff_FE.bed",row.names = FALSE,col.names = FALSE,quote = FALSE,sep="\t")





## male vs female dom plot
data<-read.table(file="both_effect_size.csv",header=TRUE,sep=",")

data$snppos<-as.numeric(as.character(sapply(strsplit(data$SNP,"_"),FUN = `[[`, 2)))
data$snpchromfull<-sapply(strsplit(data$SNP,"_"),FUN = `[[`, 1)
data$snpchrom<-ifelse(substr(data$snpchromfull,4,4)==0,substr(data$snpchromfull,5,5),substr(data$snpchromfull,4,5))

data<-data[data$kind!="",]  
data$sex<-sapply(strsplit(data$kind,"_"),FUN = `[[`, 2)

#data<-data[data$adjusted_p<0.001,]

##same as Dom plot

lengthvector<-c(174498729,95481959,105780080,90536438,92788608,96060288,68862998,28860523,161282225,125877811,111868677,101677876,
                114417674,101980477,110670232,96486271,87489397,84084598,88107222,96847506,59819933,63823863,52460201,49354470,54385492,
                55994222,45305548,41468476,43051128)
chromadd<-c(0,cumsum(lengthvector+50000000)[1:28])
#chromname<-unique(sharedeqtlrnatab$gene_chr)
chromname<-c(1:29)
chrommid<-chromadd+lengthvector/2
listval<-cbind.data.frame(chromname,lengthvector,chromadd,chrommid)
data$newposgene<-data$start+listval[match(data$chr,listval$chromname),3]
data$newpossnp<-data$snppos+listval[match(data$snpchrom,listval$chromname),3]

library(ggplot2)
ggplot(data,aes(x=newpossnp,y=newposgene,color=as.factor(chr)))+ scale_color_manual(values = rep(c("orange","dodgerblue3"), 29))+geom_point(alpha = 0.75,size=0.80)+scale_x_continuous(label =listval$chromname, breaks = listval$chrommid )+scale_y_continuous(label =listval$chromname, breaks = listval$chrommid )+theme_bw(20)+facet_wrap(~sex)



#sex specific relation
uniquesnp<-unique(data$SNP)
malespecific<-c() ##5996
femalespecific<-c() ##4407
shared<-c() ##1060
for( i in uniquesnp){
  male<-data[data$SNP==i & data$sex=="male","chr"]
  female<-data[data$SNP==i & data$sex=="female","chr"]
  
  if(length(female)==0){
    malespecific<-c(malespecific,i)
  }else if(length(male)==0){
    femalespecific<-c(femalespecific,i)
  }else{
    shared<-c(shared,i)
  }
}


datamale<-data[data$SNP%in%malespecific,]

#make table of random snp-gene pair
table<-data.frame(pid=c("ENSSSAG00000071363","ENSSSAG00000113161","ENSSSAG00000001220","ENSSSAG00000060568","ENSSSAG00000084166"),SNP=c("ssa24_14616345","ssa21_26540811","ssa03_89856471","ssa09_3632073","ssa11_89451052"))

tableexp<-Get_SNP_gene_info(data=table,feature="sex")

library(ggplot2)
ggplot(tableexp,aes(x=as.factor(Genotype),y=Expression,fill=sex))+theme_bw(20)+facet_wrap(~Gene,scales = 'free_y')+geom_boxplot()

##keep snps that regulate more than 2 genes
keep<-c()
for( i in uniquesnp){
  alldata<-data[data$SNP==i,"chr"]
  if(length(alldata)>2){
    keep<-c(keep,i)
  }
}

subdata<-data[data$SNP%in%keep,]





uniquesnp<-unique(subdata$SNP)
malespecific<-c() 
femalespecific<-c()
shared<-c() 
for( i in uniquesnp){
  male<-subdata[subdata$SNP==i & subdata$sex=="male","chr"]
  female<-subdata[subdata$SNP==i & subdata$sex=="female","chr"]
  
  if(length(female)==0){
    malespecific<-c(malespecific,i)
  }else if(length(male)==0){
    femalespecific<-c(femalespecific,i)
  }else{
    shared<-c(shared,i)
  }
}


##are most significant snp on chr2 sex specific or shared ?
datamale<-data[data$SNP%in%malespecific,]
datafemale<-data[data$SNP%in%femalespecific,]
datashared<-data[data$SNP%in%shared,]

datamale$group<-"malespe"
datafemale$group<-"femalespe"
datashared$group<-"shared"

totdata<-rbind.data.frame(datamale,datafemale,datashared,stringsAsFactors = FALSE)

ggplot(totdata,aes(x=newpossnp,y=newposgene,color=group))+geom_point(alpha = 0.75,size=1.1)+scale_x_continuous(label =listval$chromname, breaks = listval$chrommid )+scale_y_continuous(label =listval$chromname, breaks = listval$chrommid )+theme_bw(20)+facet_wrap(~sex)

totdata2<-totdata[totdata$snpchrom==2,]


table(totdata$group)
table(totdata2$group)
