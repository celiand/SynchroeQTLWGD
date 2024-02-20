#### Celian Diblasi
#### 20/02/2024
#### Investigate megative correlation ohnologs with a shared eQTL


negativedata<-results2[results2$spearman_coeff<0 & results2$regulation_0=="Shared eQTL",]

neggene<-cistransdata[cistransdata$gene%in%negativedata$gene1 | cistransdata$gene%in%negativedata$gene2,]
neggene <- neggene  %>% arrange(gene)

negativedatasubset<-negativedata[,c(1,3,4)]

negativedatasubset$snpid<-c("ssa06_12924981","ssa14_51354372","ssa14_51354372","ssa14_51354372","ssa17_62041220","ssa02_14694003","ssa14_51354372","ssa07_56591350","ssa17_62354384")

tableexpression<-Get_SNP_gene_info(negativedatasubset,feature="pair") ## function from function_Get_expression.R

tableexpression$sumexpression<-tableexpression$expression.x+tableexpression$expression.y

write.table(tableexpression,file="shared_negative_expression_01.txt",row.names = FALSE)



### make the plots




shareddata<-read.table(file="shared_negative_expression_01.txt",header=TRUE)
library(ggplot2)



gene1<-shareddata[,c("ID","id","gene1","expression.x","SNP","genotype")]
gene2<-shareddata[,c("ID","id","gene2","expression.y","SNP","genotype")]

sumgene<-data.frame(ID=shareddata$ID,id=shareddata$id,gene=shareddata$ID,expression=shareddata$sumexpression,SNP=shareddata$SNP,genotype=shareddata$genotype)
colnames(sumgene)<-c("ID","id","gene","expression","SNP","genotype")
sumgene$pair<-"Sum of copies"

colnames(gene2)<-c("ID","id","gene","expression","SNP","genotype")
colnames(gene1)<-c("ID","id","gene","expression","SNP","genotype")

gene1$pair<-"copy 1"
gene2$pair<-"copy 2"


longgene<-rbind.data.frame(gene1,gene2,sumgene,stringsAsFactors = FALSE)

### for all plots
ggplot(longgene,aes(x=as.factor(genotype),y=expression,fill=pair))+geom_boxplot()+theme_bw(25)+xlab("eQTL genotype")+scale_fill_manual(values=c("firebrick1","dodgerblue1","darkorchid1"))+facet_wrap(~ID,scales = "free_y")


### to plot a single pair
unique(longgene$ID)

longgene<-longgene[longgene$ID=="ENSSSAG00000058006_Ssal_ENSSSAG00000093636_Ssal",]
ggplot(longgene,aes(x=as.factor(genotype),y=expression,fill=pair))+geom_boxplot()+theme_bw(25)+xlab("eQTL genotype")+scale_fill_manual(values=c("firebrick1","dodgerblue1","darkorchid1"))


### check difference in expression between sum of copies for all genotypes





### Make plot of expression for both copies for one gene (dot plot and density)

## function to
## do the scatter plot
## do the boxplot
## do the statistical test of comparison on sum of expression between pairs
## table must have column ID, SNP, genotype, Expression_gene1, Expression_gene2, sumexpression, sex
Analyze_gene_pair<-function(ID,SNP,gene,table){
  library(ggplot2)
  library(egg) 
  
  ##first do the scatter plot
  subdata<-table[table$ID==ID & table$SNP==SNP,]
  titlescatter<-paste("correlation for ",gene,sep="")
  
  toplot<-ggplot(subdata,aes(x=Expression_gene1,fill=as.factor(genotype)))+geom_density(alpha=0.5)+theme_bw(20)+theme(legend.position = "none")+ggtitle(titlescatter)+xlab("Expression copy 1")
  rightplot<-ggplot(subdata,aes(x=Expression_gene2,fill=as.factor(genotype)))+geom_density(alpha=0.5)+theme_bw(20)+coord_flip()+xlab("Expression copy 2")
  scatter<-ggplot(subdata,aes(x=Expression_gene1,y=Expression_gene2,color=as.factor(genotype)))+geom_point()+theme_bw(20)+theme(legend.position = "none")+xlab("Expression copy 1")+ylab("Expression copy 2")
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
  
  
  
  scatterplot<-ggarrange(toplot, empty, scatter, rightplot, 
                         ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))
  filename_scat<-paste("Scatter_plot_",gene,".pdf",sep="")
  path="C:/Users/cedi/Downloads/Jupyte_R_not_working/"
  ggsave(filename = filename_scat,plot=scatterplot,path=path)
  
  
  ##then do the boxplot
  
  gene1<-table[,c("ID","id","gene1","Expression_gene1","SNP","genotype")]
  gene2<-table[,c("ID","id","gene2","Expression_gene2","SNP","genotype")]
  
  sumgene<-data.frame(ID=table$ID,id=table$id,gene=table$ID,expression=table$sumexpression,SNP=table$SNP,genotype=table$genotype)
  colnames(sumgene)<-c("ID","id","gene","expression","SNP","genotype")
  sumgene$pair<-"Sum of copies"
  
  colnames(gene2)<-c("ID","id","gene","expression","SNP","genotype")
  colnames(gene1)<-c("ID","id","gene","expression","SNP","genotype")
  
  gene1$pair<-"copy 1"
  gene2$pair<-"copy 2"
  
  
  longgene<-rbind.data.frame(gene1,gene2,sumgene,stringsAsFactors = FALSE)
  
  datasub<-longgene[longgene$ID==ID & longgene$SNP==SNP,]
  
  boxplot<-ggplot(datasub,aes(x=as.factor(genotype),y=expression,fill=pair))+geom_boxplot()+theme_bw(25)+xlab("eQTL genotype")+scale_fill_manual(values=c("firebrick1","dodgerblue1","darkorchid1"))
  filename_box<-paste("Boxplot_",gene,".pdf",sep="")
  path="C:/Users/cedi/Downloads/Jupyte_R_not_working/" ## indicate desired path
  ggsave(filename = filename_box,plot=boxplot,path=path)
  
  
  ##then do the test
  
  model <- aov(sumexpression ~ genotype, data = subdata)
  return(model)
}


### then run the function

fulltabneg<-read.table("expression_table_8negativepair.txt",header=TRUE) ### see end of this script on generation of this table
colnames(fulltabneg)[26:27]<-c("Expression_gene1","Expression_gene2")
fulltabneg<-fulltabneg[,-c(6:24)]

fulltabneg<-read.table("expression_table_9negativepairdistinctshared.txt",header=TRUE)
colnames(fulltabneg)[26:27]<-c("Expression_gene1","Expression_gene2")
fulltabneg<-fulltabneg[,-c(6:24)]



### with mfsd13al gene
mfsd13al<-Analyze_gene_pair(ID="ENSSSAG00000001778_Ssal_ENSSSAG00000045167_Ssal",SNP="ssa14_51370830",gene="major facilitator superfamily domain containing 13a-like",table=fulltabneg)
summary(mfsd13al)
sumary_mfsd13al<-summary(mfsd13al)


### adjust p-values
pval<-c(sumary_mfsd13al[[1]]$`Pr(>F)`[1],summary_idh3b[[1]]$`Pr(>F)`[1],summary_mterf3[[1]]$`Pr(>F)`[1],summary_napepld[[1]]$`Pr(>F)`[1],summary_hacl1[[1]]$`Pr(>F)`[1],summary_LOC106570013[[1]]$`Pr(>F)`[1],summary_LOC106597139[[1]]$`Pr(>F)`[1],summary_cpm[[1]]$`Pr(>F)`[1])

p.adjust(pval,method="bonferroni")

##### -- refine shared negative gene pair -- #####

transdata<-read.table(file="LLsmolt.trans.adjust.hits.txt",header=FALSE)
#head(transdata)
cisdata<-read.table(file="LLsmolt.cisnominal.txt",header=FALSE)



ENSSSAG00000120807<-transdata[transdata$V1=="ENSSSAG00000120807",]
ENSSSAG00000113858<-transdata[transdata$V1=="ENSSSAG00000113858",]

intersect(ENSSSAG00000113858$V4,ENSSSAG00000120807$V4)

ENSSSAG00000011938<-transdata[transdata$V1=="ENSSSAG00000011938",]
ENSSSAG00000108018<-transdata[transdata$V1=="ENSSSAG00000108018",]

obj<-intersect(ENSSSAG00000011938$V4,ENSSSAG00000108018$V4)


### identify gene in distinct eQTL that are present in the above files

negativedistinct<-results2[results2$regulation_0=="Distinct eQTL" & results2$spearman_coeff<=0,]
negativedistinct<-as.data.frame(negativedistinct)




for ( i in 1:length(negativedistinct$ID)){
  gene1<-negativedistinct[i,"gene1"]
  gene2<-negativedistinct[i,"gene2"]
  
  snpgene1<-cistransdata[cistransdata$gene==gene1,"SNP"]
  snpgene2<-cistransdata[cistransdata$gene==gene2,"SNP"]
  
  transgene1<-transdata[transdata$V1==gene1,]
  transgene2<-transdata[transdata$V1==gene2,]
  
  if(length(transgene1$V1)!=0 && length(transgene2$V1)!=0){
    intersecttrans<-intersect(transgene1$V4,snpgene2)
    if(length(intersecttrans)>=0){
      
    }else{
      intersecttrans<-intersect(transgene2$V4,snpgene1)
    }
  }else{
    intersecttrans<-NULL
  }
  
  cisgene1<-cisdata[cisdata$V1==gene1,]
  cisgene2<-cisdata[cisdata$V1==gene2,]
  
  if(length(cisgene1$V1)!=0 && length(cisgene2$V1)!=0){
    intersectcis<-intersect(cisgene1$V8,snpgene2)
    if(length(intersecttrans)>=0){
      
    }else{
      intersectcis<-intersect(cisgene2$V8,snpgene1)
    }
  }else{
    intersectcis<-NULL
  }
  
  if(length(intersectcis)!=0 & length(intersecttrans)!=0){
    val<-"shared"
    sampled_snp <- sample(intersecttrans, size = 1, replace = TRUE)
  }else if(length(intersectcis)!=0){
    val<-"shared"
    sampled_snp <- sample(intersectcis, size = 1, replace = TRUE)
  }else if(length(intersecttrans)!=0){
    val<-"shared"
    sampled_snp <- sample(intersecttrans, size = 1, replace = TRUE)
    
  }else{
    val<-"distinct"
    sampled_snp<-NA
  }
  
  negativedistinct[i,"indepth_test"]<-val
  negativedistinct[i,"commonsnp"]<-sampled_snp
}



negativedistinct$indepth_test

write.table(negativedistinct,file="negative_distinct_withindepth_02.txt",row.names = FALSE)

negativedistinct<-read.table(file="negative_distinct_withindepth_02.txt",header=TRUE)

negativedistinctshared<-negativedistinct[negativedistinct$indepth_test=="shared",]

negativedistinctsharedsub<-negativedistinctshared[,c(1,3,4,27)]
colnames(negativedistinctsharedsub)[4]<-"snpid"


negativedistinctexpression<-Get_SNP_gene_info(negativedistinctsharedsub,feature="pair")
negativedistinctexpression$sumexpression<-negativedistinctexpression$expression.x+negativedistinctexpression$expression.y

write.table(negativedistinctexpression,file="distinct_negative_shared_expression_01.txt",row.names = FALSE)




distinctdata<-read.table(file="distinct_negative_shared_expression_01.txt",header=TRUE)
library(ggplot2)



gene1<-distinctdata[,c("ID","id","gene1","expression.x","SNP","genotype")]
gene2<-distinctdata[,c("ID","id","gene2","expression.y","SNP","genotype")]

sumgene<-data.frame(ID=distinctdata$ID,id=distinctdata$id,gene=distinctdata$ID,expression=distinctdata$sumexpression,SNP=distinctdata$SNP,genotype=distinctdata$genotype)
colnames(sumgene)<-c("ID","id","gene","expression","SNP","genotype")
sumgene$pair<-"Sum of copies"

colnames(gene2)<-c("ID","id","gene","expression","SNP","genotype")
colnames(gene1)<-c("ID","id","gene","expression","SNP","genotype")

gene1$pair<-"copy 1"
gene2$pair<-"copy 2"


longgene<-rbind.data.frame(gene1,gene2,sumgene,stringsAsFactors = FALSE)


ggplot(longgene,aes(x=as.factor(genotype),y=expression,fill=pair))+geom_boxplot()+theme_bw(25)+xlab("eQTL genotype")+scale_fill_manual(values=c("firebrick1","dodgerblue1","darkorchid1"))+facet_wrap(~ID,scales = "free_y")

unique(longgene$ID)

longgene<-longgene[longgene$ID=="ENSSSAG00000009565_Ssal_ENSSSAG00000121915_Ssal",]
ggplot(longgene,aes(x=as.factor(genotype),y=expression,fill=pair))+geom_boxplot()+theme_bw(25)+xlab("eQTL genotype")+scale_fill_manual(values=c("firebrick1","dodgerblue1","darkorchid1"))


all_pval<-c()
for( i in unique(distinctdata$ID)){
  subdata<-distinctdata[distinctdata$ID==i,]
  model <- aov(sumexpression ~ genotype, data = subdata)
  sum_model<-summary(model)
  pval<-sum_model[[1]]$`Pr(>F)`[1]
  all_pval<-c(all_pval,pval)
}

p.adjust(all_pval,method="bonferroni")


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






### -- in depth analysis of the 18 negative correlated shared gene pair -- ###
negativedata<-results2[results2$spearman_coeff<0 & results2$regulation_0=="Shared eQTL",]
negativedata_subset<-negativedata[,c("ID","gene1","gene2")]
negativedata_subset$commonsnp<-c("ssa06_12924981","ssa14_51354372","ssa14_51354372","ssa14_51354372","ssa17_62041220","ssa02_14694003","ssa14_51354372","ssa07_56591350","ssa17_62354384")
negativedata_subset$status<-"trueshared"

negativedistinct<-read.table(file="negative_distinct_withindepth_02.txt",header=TRUE)
negativedistinctshared<-negativedistinct[negativedistinct$indepth_test=="shared",]
negativedistinctshared_subset<-negativedistinctshared[,c("ID","gene1","gene2","commonsnp")]
negativedistinctshared_subset$status<-"distinctshared"

all_shared<-rbind.data.frame(negativedata_subset,negativedistinctshared_subset,stringsAsFactors = FALSE)

write.table(all_shared,"allshared_info01.txt",row.names=FALSE)




#### WGD difference between LORe and AORe 



genetab<-read.table(file="ss4r_dups_and_singletons_ENSrapid_convPipeline.tsv",header=TRUE)

all_shared<-read.table("allshared_info01.txt",header=TRUE)

genepair<-read.table(file="genetabgenepair04.txt",header=TRUE)
region_results_exp<-merge(genetab,genepair,by="ID")
region_results_exp<-region_results_exp[!is.na(region_results_exp$redip.class),]

region_results_exp<-region_results_exp[,1:8]
region_results_exp$commonsnp<-NA
region_results_exp$status<-NA
region_results_exp$data<-"expected"

region_results<-merge(all_shared,genetab,by="ID")
region_results<-region_results[!is.na(region_results$redip.class),]
region_results$data<-"sharedneg"

alldata<-rbind.data.frame(region_results,region_results_exp,stringsAsFactors = FALSE)

library(ggplot2)
library(dplyr)

dfprop<- alldata %>% group_by(data, redip.class) %>% 
  summarise(total_count=n(),.groups = 'drop') %>%
  as.data.frame() %>% group_by(data) %>% 
  mutate(prop = total_count / sum(total_count))

ggplot(dfprop,aes(x=data,y=prop,fill=redip.class))+geom_bar(position = "fill", stat = "identity")+theme_bw(25)+ylab("Proportion")+scale_fill_manual(values=c("#5387b4","#db6d48"))


chidata<-data.frame(LORe=c(1984,12),AORe=c(5895,2))

chisq.test(chidata)



### check if there are any differences between the 4 categories and AORe and LORe regions

genetab<-read.table(file="ss4r_dups_and_singletons_ENSrapid_convPipeline.tsv",header=TRUE)


region_results<-merge(genetab,results3,by="ID")
region_results<-region_results[!is.na(region_results$redip.class),]

region_AORE<-region_results[region_results$redip.class=="AORe",]
ggplot(region_AORE,aes(x=spearman_coeff,fill=regulation_0))+geom_histogram(color="black")+theme_bw(20)+facet_wrap(~regulation_0,scales = 'free_y')+xlab("Spearman coefficient")+scale_fill_manual(values=c("#8CE86D","#F9B27C","#D186D8","#99B8EA"))+geom_vline(data = df2, mapping = aes(xintercept = Mean),color="firebrick3",size=1.5)

region_LORE<-region_results[region_results$redip.class=="LORe",]
ggplot(region_LORE,aes(x=spearman_coeff,fill=regulation_0))+geom_histogram(color="black")+theme_bw(20)+facet_wrap(~regulation_0,scales = 'free_y')+xlab("Spearman coefficient")+scale_fill_manual(values=c("#8CE86D","#F9B27C","#D186D8","#99B8EA"))+geom_vline(data = df2, mapping = aes(xintercept = Mean),color="firebrick3",size=1.5)



dfcount <- region_results%>% group_by(regulation_0,redip.class) %>% 
  summarise(total_count=n(),.groups = 'drop') %>%
  as.data.frame()  %>% group_by(redip.class)%>%
  mutate(prop = total_count / sum(total_count))


#ggplot(dfcount,aes(x="",y=prop,fill=regulation_0))+geom_bar(stat="identity",width=1)+
  #coord_polar("y",start=0)+theme_void()+scale_fill_manual(values=c("#F9B27C","#99B8EA","#8CE86D","#D186D8"))+
  #theme(legend.text = element_text(size=25))+facet_wrap(~redip.class)

ggplot(dfcount,aes(x=regulation_0,y=prop,fill=redip.class))+geom_bar(position = "fill", stat = "identity")+theme_bw(25)+ylab("Proportion")+scale_fill_manual(values=c("#5387b4","#db6d48"))+theme(axis.text.x = element_text(angle = 45, hjust=1))

chidata<-data.frame(LORe=c(59,241,342,194),AORe=c(97,603,936,548))

chisq.test(chidata) ## X-squared = 9.7479, df = 3, p-value = 0.02084



###### make expression table for 9 true negative gene pair ######



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



data=negativedatasubset
library(tidyverse)

Expression <- read_table2("Synchro_mt_TPM.bed")    %>% select(-c(1,2,3,5,6))
Genotype <- read_table2("imputed_0.01_2730_smolts.traw")   %>% select(-c(1,3,4,5,6))
Covariate <-  read_table2("covariate_synchro.txt")
## fix the name in genotype
names(Genotype)[2:ncol(Genotype)] <- str_remove(names(Genotype)[2:ncol(Genotype)], "0_")
covar_kept<-Covariate[c(3,5),]
id<-colnames(covar_kept)
tcovar<-t(covar_kept)
long_covar<-cbind.data.frame(id,tcovar)
long_covar<-long_covar[-1,]
colnames(long_covar)<-c("id","sex","LL")

neggeneexpp1<- Expression[Expression$pid%in%data$gene1,]
neggeneexpp2<- Expression[Expression$pid%in%data$gene2,]

neggene_long1<- neggeneexpp1 %>%
  pivot_longer(cols = -pid, names_to = "id", values_to = "expression")

neggene_long2<- neggeneexpp2 %>%
  pivot_longer(cols = -pid, names_to = "id", values_to = "expression")

long_covar_LL<-long_covar[,c(1,3)]

mergedpair1neg <- merge(neggene_long1 , long_covar, by = "id")
mergedpair2neg <- merge(neggene_long2 , long_covar, by = "id")

##filter

mergedpair1neg_filter<- mergedpair1neg[mergedpair1neg$LL=="LL",]
mergedpair2neg_filter<- mergedpair2neg[mergedpair2neg$LL=="LL",]

##can keep more column if wanted
mergedpair1neg_filter<-mergedpair1neg_filter[,1:3]
mergedpair2neg_filter<-mergedpair2neg_filter[,1:3]

# Left join table 1 with table 2 based on gene1 and pid
joined_table1_2neg <- left_join(genepair, mergedpair1neg_filter, by = c("gene1" = "pid"))
joined_table1_2neg<-joined_table1_2neg[!is.na(joined_table1_2neg$expression),]

# Left join result with table 3 based on gene2 and pid
joined_table1_2_3neg <- left_join(joined_table1_2neg, mergedpair2neg_filter, by = c("id"="id","gene2" = "pid"))
joined_table1_2_3neg <- joined_table1_2_3neg [!is.na(joined_table1_2_3neg$expression.y),]

colnames(joined_table1_2_3neg)[6:7]<-c("Expression_gene1","Expression_gene2")

gtkept<-Genotype[Genotype$SNP%in%data$snpid,]

geno_table_neg <- gtkept %>% pivot_longer(cols = -SNP, names_to = "id", values_to = "genotype")

gene_tab_neg_merge <- merge(geno_table_neg  , long_covar, by = "id")  
gene_tab_neg_merge_filter<- gene_tab_neg_merge[gene_tab_neg_merge$LL=="LL",]

mergegt<-merge(gene_tab_neg_merge_filter,data,by.x="SNP",by.y = "snpid")
mergegt<-mergegt[,-c(7,8)]

fulltabneg<-merge(joined_table1_2_3neg,mergegt,by=c("ID","id"))
fulltabneg$sumexpression<-fulltabneg$expression.x+fulltabneg$expression.y
final_table<-fulltabneg


write.table(final_table,file="expression_table_8negativepair.txt",row.names=FALSE)



### and for 9 refined negative gene pair

###same but get data from

negativedistinct<-read.table(file="negative_distinct_withindepth_02.txt",header=TRUE)

negativedistinctshared<-negativedistinct[negativedistinct$indepth_test=="shared",]
