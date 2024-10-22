#### Celian Diblasi
#### 15/03/2024
#### Investigate negative correlated ohnologs pair with a shared eQTL


library(ggplot2)
library(dplyr)
library(tidyr)

gene_pair_tab<-read.table(file="expression_gene_pair02.txt",header=TRUE)

##coefficient

## import data
result <- gene_pair_tab %>% 
  group_by(ID) %>% 
  summarize(spearman_coeff = cor(Expression_gene1, Expression_gene2, method = "spearman"))


genepair<-read.table(file="genetabgenepair05.txt",header=TRUE)

results2<-left_join(result, genepair, by = c("ID"="ID"))



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


## function to have genotype of a SNP and expression of a gene/ a pair of gene
## feature: sex (male vs female); pair
## data format: 
## if feature is sex: col pid with gene ID, col SNP with snp position in the following format: chrom_position
## if feature is pair: col ID with pair ID, col gene1 and gene2 with genes ID, snp position (snpid) in the following format: chrom_position
## expression: raw (default) or scaled
Get_SNP_gene_info<-function(data,feature,expression="raw"){
  library(tidyverse)
  Genotype <- read_table2("/mnt/project/MSLab/Marie/2022_Synchrosmolt/QTLTools/imputed_0.01_2730_smolts.traw")   %>% select(-c(1,3,4,5,6))
  Expression <- read_table2("Synchro_mt_TPM.bed")    %>% select(-c(1,2,3,5,6))
  Covariate <-  read_table2("covariate_synchro.txt")
  
  ## fix the name in genotype
  names(Genotype)[2:ncol(Genotype)] <- str_remove(names(Genotype)[2:ncol(Genotype)], "0_")
 
   ##covar table treatment is common in both analysis
  covar_kept<-Covariate[c(3,5),]
  id<-colnames(covar_kept)
  tcovar<-t(covar_kept)
  long_covar<-cbind.data.frame(id,tcovar)
  long_covar<-long_covar[-1,]
  colnames(long_covar)<-c("id","sex","LL")
  if(expression=="scaled"){
    # scale the expression
    Expression[, 2:ncol(Expression)] <- t(apply(Expression[, 2:ncol(Expression)], 1, scale))
  }else{
  }
  if(feature=="sex"){
    expr_kept<-Expression[Expression$pid%in%data$pid,]
    gt_kept<-Genotype[Genotype$SNP%in%data$SNP,]
    # Pivot genotype table to long format
    geno_table_long <- gt_kept %>%
      pivot_longer(cols = -SNP, names_to = "id", values_to = "genotype")
    #Pivot expression table to long format
    expr_table_long <- expr_kept %>%
      pivot_longer(cols = -pid, names_to = "id", values_to = "expression")
    # Join the three tables based on the 'id' column
    merged <- merge(geno_table_long, expr_table_long, by = "id")
    merged <- merge(merged, long_covar, by = "id")
    # Filter to only keep the rows where 'gene' and 'snp' match the values in 'infokeep'
    merged_filter <- subset(merged, pid %in% data$pid & SNP %in% data$SNP)
    merged_filter2<- merged_filter[merged_filter$LL=="LL",]
    merged_filter2$combination<-paste(merged_filter2$SNP,merged_filter2$pid,sep="_")
    data$combination<-paste(data$SNP,data$pid,sep="_")
    merged_filter3<-merged_filter2[merged_filter2$combination%in%data$combination,]
    colnames(merged_filter3)<-c("ID","SNP","Genotype","Gene","Expression","sex","condition","name")
    final_table<-merged_filter3
  }else if(feature=="pair"){
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
    if (exists("genepair")) {
      joined_table1_2neg <- left_join(genepair, mergedpair1neg_filter, by = c("gene1" = "pid"))
      joined_table1_2neg<-joined_table1_2neg[!is.na(joined_table1_2neg$expression),]
    } else {
      joined_table1_2neg <- left_join(data, mergedpair1neg_filter, by = c("gene1" = "pid"))
      joined_table1_2neg<-joined_table1_2neg[!is.na(joined_table1_2neg$expression),]
    }
    
    
    # Left join result with table 3 based on gene2 and pid
    joined_table1_2_3neg <- left_join(joined_table1_2neg, mergedpair2neg_filter, by = c("id"="id","gene2" = "pid"))
    joined_table1_2_3neg <- joined_table1_2_3neg [!is.na(joined_table1_2_3neg$expression.y),]
    
    colnames(joined_table1_2_3neg)[6:7]<-c("Expression_gene1","Expression_gene2")
    
    gtkept<-Genotype[Genotype$SNP%in%data$snpid,]
    
    geno_table_neg <- gtkept %>% pivot_longer(cols = -SNP, names_to = "id", values_to = "genotype")
    
    #long_covar$id<-paste("0_",long_covar$id,sep="")
    gene_tab_neg_merge <- merge(geno_table_neg  , long_covar, by = "id")  
    gene_tab_neg_merge_filter<- gene_tab_neg_merge[gene_tab_neg_merge$LL=="LL",]
    
    mergegt<-merge(gene_tab_neg_merge_filter,data,by.x="SNP",by.y = "snpid")
    mergegt<-mergegt[,-c(7,8)]
    
    fulltabneg<-merge(joined_table1_2_3neg,mergegt,by=c("ID","id"))
    #fulltabneg$sumexpression<-fulltabneg$Expression_gene1+fulltabneg$Expression_gene2
    final_table<-fulltabneg
    
  }else{
    stop("Wrong feature selected. Possibles values are 'pair' or 'sex'.")
  }
  return(final_table)
}





sharedpair<-negativedata[,c(1,3,4,6)]
colnames(sharedpair)[4]<-"snpid"

distinctpair<-negativedistinctshared[,c(1,3,4,27)]
colnames(distinctpair)[4]<-"snpid"


negativedatasubset<-rbind.data.frame(sharedpair,distinctpair,stringsAsFactors = FALSE)

### get expression values
tableexpression<-Get_SNP_gene_info(negativedatasubset,feature="pair",expression = "raw")


tableexpressionsub<-tableexpression[,c(1:6,47:51)]
tableexpressionsub <- tableexpressionsub %>%
  mutate(genotype = case_when(
    genotype == 0 ~ "AA",
    genotype == 1 ~ "Aa",
    genotype == 2 ~ "aa"
  ))

tableexpressionsub$sumexpression<-tableexpressionsub$expression.x+tableexpressionsub$expression.y

## save data
write.table(tableexpressionsub,"expression_table_18_negativesharedpair02.txt",row.names = FALSE)



##### -- Make basic plots and analysis -- #####
tableexpressionsub<-read.table(file="expression_table_18_negativesharedpair02.txt",header=TRUE)


gene1<-tableexpressionsub[,c(1,2,3,7,9,10)]
gene2<-tableexpressionsub[,c(1,2,4,8,9,10)]

sumgene<-data.frame(ID=tableexpressionsub$ID,id=tableexpressionsub$id,gene=tableexpressionsub$ID,expression=tableexpressionsub$sumexpression,SNP=tableexpressionsub$SNP,genotype=tableexpressionsub$genotype)
colnames(sumgene)<-c("ID","id","gene","expression","SNP","genotype")
sumgene$pair<-"Sum of copies"

colnames(gene2)<-c("ID","id","gene","expression","SNP","genotype")
colnames(gene1)<-c("ID","id","gene","expression","SNP","genotype")

gene1$pair<-"copy 1"
gene2$pair<-"copy 2"

longgene<-rbind.data.frame(gene1,gene2,sumgene,stringsAsFactors = FALSE)

ggplot(longgene,aes(x=as.factor(genotype),y=expression,fill=pair))+geom_boxplot()+theme_bw(25)+xlab("eQTL genotype")+scale_fill_manual(values=c("firebrick1","dodgerblue1","darkorchid1"))+facet_wrap(~ID,scales = "free_y")


#### plot a single pair
longgene<-longgene[longgene$ID=="ENSSSAG00000009565_Ssal_ENSSSAG00000121915_Ssal",]
ggplot(longgene,aes(x=as.factor(genotype),y=expression,fill=pair,color=pair))+geom_boxplot(color="black",width=0.6, position=position_dodge(width=0.9))+theme_bw(25)+xlab("eQTL genotype")+ylab("Expression (TPM)")+scale_fill_manual(values=c("firebrick1","dodgerblue1","darkorchid1"))+geom_point(alpha=0.3,position=position_jitterdodge(dodge.width=0.9, jitter.width=0.2),size=0.8)+scale_color_manual(values=c("firebrick1","dodgerblue1","darkorchid1"))




### function to plot for a single pair
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
  path="/mnt/users/cedi/"
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
  
  boxplot<-ggplot(datasub,aes(x=as.factor(genotype),y=expression,fill=pair,color=pair))+geom_boxplot(color="black",width=0.6, position=position_dodge(width=0.9))+theme_bw(20)+xlab("eQTL genotype")+ylab("Expression (TPM)")+scale_fill_manual(values=c("firebrick1","dodgerblue1","darkorchid1"))+geom_point(alpha=0.2,position=position_jitterdodge(dodge.width=0.9, jitter.width=0.2),size=0.6)+scale_color_manual(values=c("firebrick1","dodgerblue1","darkorchid1"))
  filename_box<-paste("Boxplot_",gene,".pdf",sep="")
  path="/mnt/users/cedi/"
  ggsave(filename = filename_box,plot=boxplot,path=path)
  
  
  ##then test if sum of expression is different between genotype
  
  model <- aov(sumexpression ~ genotype, data = subdata)
  return(model)
}

all_pairs<-read.table(file="expression_table_18_negativesharedpair02.txt",header=TRUE)


colnames(all_pairs)[7]<-"Expression_gene1"
colnames(all_pairs)[8]<-"Expression_gene2"

all_pairs[all_pairs$ID=="ENSSSAG00000001739_Ssal_ENSSSAG00000041942_Ssal","Gene_name"]<-"zgc_112148"
all_pairs[all_pairs$ID=="ENSSSAG00000001739_Ssal_ENSSSAG00000041942_Ssal","Gene_code"]<-"zgc"

all_pairs[all_pairs$ID=="ENSSSAG00000001778_Ssal_ENSSSAG00000045167_Ssal","Gene_name"]<-"major facilitator superfamily domain containing 13a-like"
all_pairs[all_pairs$ID=="ENSSSAG00000001778_Ssal_ENSSSAG00000045167_Ssal","Gene_code"]<-"mfsd13al"

all_pairs[all_pairs$ID=="ENSSSAG00000002173_Ssal_ENSSSAG00000114008_Ssal","Gene_name"]<-"STT3 oligosaccharyltransferase complex catalytic subunit B"
all_pairs[all_pairs$ID=="ENSSSAG00000002173_Ssal_ENSSSAG00000114008_Ssal","Gene_code"]<-"stt3b"

all_pairs[all_pairs$ID=="ENSSSAG00000009565_Ssal_ENSSSAG00000121915_Ssal","Gene_name"]<-"Src kinase associated phosphoprotein 2"
all_pairs[all_pairs$ID=="ENSSSAG00000009565_Ssal_ENSSSAG00000121915_Ssal","Gene_code"]<-"skap2"

all_pairs[all_pairs$ID=="ENSSSAG00000011938_Ssal_ENSSSAG00000108018_Ssal","Gene_name"]<-"collectin-10-like"
all_pairs[all_pairs$ID=="ENSSSAG00000011938_Ssal_ENSSSAG00000108018_Ssal","Gene_code"]<-"LOC106596782"

all_pairs[all_pairs$ID=="ENSSSAG00000036507_Ssal_ENSSSAG00000116949_Ssal","Gene_name"]<-"secretagogin"
all_pairs[all_pairs$ID=="ENSSSAG00000036507_Ssal_ENSSSAG00000116949_Ssal","Gene_code"]<-"segn"

all_pairs[all_pairs$ID=="ENSSSAG00000043068_Ssal_ENSSSAG00000054986_Ssal","Gene_name"]<-"isocitrate dehydrogenase 3 non-catalytic subunit beta"
all_pairs[all_pairs$ID=="ENSSSAG00000043068_Ssal_ENSSSAG00000054986_Ssal","Gene_code"]<-"idh3b"

all_pairs[all_pairs$ID=="ENSSSAG00000044129_Ssal_ENSSSAG00000117896_Ssal","Gene_name"]<-"ELOVL fatty acid elongase 8b"
all_pairs[all_pairs$ID=="ENSSSAG00000044129_Ssal_ENSSSAG00000117896_Ssal","Gene_code"]<-"elovl8a"

all_pairs[all_pairs$ID=="ENSSSAG00000047981_Ssal_ENSSSAG00000057888_Ssal","Gene_name"]<-"mitochondrial transcription termination factor 3"
all_pairs[all_pairs$ID=="ENSSSAG00000047981_Ssal_ENSSSAG00000057888_Ssal","Gene_code"]<-"mterf3"

all_pairs[all_pairs$ID=="ENSSSAG00000054393_Ssal_ENSSSAG00000117250_Ssal","Gene_name"]<-"N-acyl phosphatidylethanolamine phospholipase D"
all_pairs[all_pairs$ID=="ENSSSAG00000054393_Ssal_ENSSSAG00000117250_Ssal","Gene_code"]<-"napepld"

all_pairs[all_pairs$ID=="ENSSSAG00000058006_Ssal_ENSSSAG00000093636_Ssal","Gene_name"]<-"2-hydroxyacyl-CoA lyase 1"
all_pairs[all_pairs$ID=="ENSSSAG00000058006_Ssal_ENSSSAG00000093636_Ssal","Gene_code"]<-"hacl1"

all_pairs[all_pairs$ID=="ENSSSAG00000066297_Ssal_ENSSSAG00000078929_Ssal","Gene_name"]<-"E3 ubiquitin-protein ligase RNF19B-like"
all_pairs[all_pairs$ID=="ENSSSAG00000066297_Ssal_ENSSSAG00000078929_Ssal","Gene_code"]<-"LOC106570013"

all_pairs[all_pairs$ID=="ENSSSAG00000083156_Ssal_ENSSSAG00000116300_Ssal","Gene_name"]<-"Ribosome binding factor A"
all_pairs[all_pairs$ID=="ENSSSAG00000083156_Ssal_ENSSSAG00000116300_Ssal","Gene_code"]<-"rbfa"

all_pairs[all_pairs$ID=="ENSSSAG00000086178_Ssal_ENSSSAG00000087231_Ssal","Gene_name"]<-"tetraspanin-8-like"
all_pairs[all_pairs$ID=="ENSSSAG00000086178_Ssal_ENSSSAG00000087231_Ssal","Gene_code"]<-"LOC106597139"

all_pairs[all_pairs$ID=="ENSSSAG00000090152_Ssal_ENSSSAG00000110471_Ssal","Gene_name"]<-"CC068 protein"
all_pairs[all_pairs$ID=="ENSSSAG00000090152_Ssal_ENSSSAG00000110471_Ssal","Gene_code"]<-"cc068"

all_pairs[all_pairs$ID=="ENSSSAG00000091200_Ssal_ENSSSAG00000121190_Ssal","Gene_name"]<-"Protein associated with LIN7 2, MAGUK p55 family member"
all_pairs[all_pairs$ID=="ENSSSAG00000091200_Ssal_ENSSSAG00000121190_Ssal","Gene_code"]<-"PALS2"

all_pairs[all_pairs$ID=="ENSSSAG00000092357_Ssal_ENSSSAG00000109189_Ssal","Gene_name"]<-"carboxypeptidase M"
all_pairs[all_pairs$ID=="ENSSSAG00000092357_Ssal_ENSSSAG00000109189_Ssal","Gene_code"]<-"cpm"

all_pairs[all_pairs$ID=="ENSSSAG00000113858_Ssal_ENSSSAG00000120807_Ssal","Gene_name"]<-"Eukaryotic translation initiation factor 3, subunit M"
all_pairs[all_pairs$ID=="ENSSSAG00000113858_Ssal_ENSSSAG00000120807_Ssal","Gene_code"]<-"eif3m"


### for all pair, do the plots and save the test results
pval<-c()
for( i in unique(all_pairs$ID)){
  dataanalyzed<-Analyze_gene_pair(ID=i,SNP=unique(all_pairs[all_pairs$ID==i,"SNP"]),gene=unique(all_pairs[all_pairs$ID==i,"Gene_name"]),table=all_pairs)
  summary_data<-summary(dataanalyzed)
  pval<-c(pval,summary_data[[1]]$`Pr(>F)`[1])
}

### adjust p values for multiple testing






##### -- LORe and AORe analysis -- ####

### get duplicated information

genetab<-read.table(file="ss4r_dups_and_singletons_ENSrapid_convPipeline.tsv",header=TRUE)

all_shared<-read.table("allshared_info01.txt",header=TRUE)

genepair<-read.table(file="genetabgenepair05.txt",header=TRUE)

##first get AORe and LORe status in all ohnologs pairs
region_results_exp<-merge(genetab,genepair,by="ID")
region_results_exp<-region_results_exp[!is.na(region_results_exp$redip.class),]

region_results_exp<-region_results_exp[,1:8]
region_results_exp$commonsnp<-NA
region_results_exp$status<-NA
region_results_exp$data<-"expected"

## then same but only for the 18 pairs in dosage compensation
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



### check if there are any differences between the ohnolog pairs in AORe and LORe regions depending on their status

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



ggplot(dfcount,aes(x=regulation_0,y=prop,fill=redip.class))+geom_bar(position = "fill", stat = "identity")+theme_bw(25)+ylab("Proportion")+scale_fill_manual(values=c("#5387b4","#db6d48"))+theme(axis.text.x = element_text(angle = 45, hjust=1))

chidata<-data.frame(LORe=c(59,241,342,194),AORe=c(97,603,936,548))

chisq.test(chidata) ## X-squared = 9.7479, df = 3, p-value = 0.02084
p.adjust(pval,method="bonferroni")