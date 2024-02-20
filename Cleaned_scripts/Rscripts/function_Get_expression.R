## function to have genotype of a SNP and expression of a gene/ a pair of gene
## feature: sex (male vs female); pair
## data format: 
## if feature is sex: col pid with gene ID, col SNP with snp position in the following format: chrom_position
## if feature is pair: col ID with pair ID, col gene1 and gene2 with genes ID, snp position (snpid) in the following format: chrom_position
## expression: raw (default) or scaled
Get_SNP_gene_info<-function(data,feature,expression="raw"){
  library(tidyverse)
  #Expression <- read_table2("/mnt/project/MSLab/Marie/2022_Synchrosmolt/QTLTools/Synchro_mt_TPM.bed")    %>% select(-c(1,2,3,5,6))
  Genotype <- read_table2("/mnt/project/MSLab/Marie/2022_Synchrosmolt/QTLTools/imputed_0.01_2730_smolts.traw")   %>% select(-c(1,3,4,5,6))
  Expression <- read_table2("Synchro_mt_TPM.bed")    %>% select(-c(1,2,3,5,6))
  #Genotype <- read_table2("imputed_0.01_2730_smolts.traw")   %>% select(-c(1,3,4,5,6))
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
    joined_table1_2neg <- left_join(genepair, mergedpair1neg_filter, by = c("gene1" = "pid"))
    joined_table1_2neg<-joined_table1_2neg[!is.na(joined_table1_2neg$expression),]
    
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
    fulltabneg$sumexpression<-fulltabneg$Expression_gene1+fulltabneg$Expression_gene2
    final_table<-fulltabneg
    
  }else{
    stop("Wrong feature selected. Possibles values are 'pair' or 'sex'.")
  }
  return(final_table)
}