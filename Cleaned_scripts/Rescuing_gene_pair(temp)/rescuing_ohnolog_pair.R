#### Celian Diblasi
#### 15/03/2024
#### Rescuing ohnologs pair with distinct eQTL but actually sharing eQTL



##### -- refine shared negative gene pair -- #####

## import data of all eQTL


transdata<-read.table(file="LLsmolt.trans.adjust.hits.txt",header=FALSE)

cisdata<-read.table(file="LLsmolt.cisnominal.txt",header=FALSE)



### identify gene in distinct eQTL that are present in the above files

library(dplyr)
library(tidyr)

### import previous pair status
gene_pair_tab<-read.table(file="expression_gene_pair02.txt",header=TRUE)

##coefficient

# calculate spearman coefficient for each condition
result <- gene_pair_tab %>% 
  group_by(ID) %>% 
  summarize(spearman_coeff = cor(Expression_gene1, Expression_gene2, method = "spearman"))


genepair<-read.table(file="genetabgenepair05.txt",header=TRUE)

results2<-left_join(result, genepair, by = c("ID"="ID"))



negativedistinct<-results2[results2$regulation_0=="Distinct eQTL" & results2$spearman_coeff<=0,]
negativedistinct<-as.data.frame(negativedistinct)


### for each of the 53 pair with distinct eQTL and negative correlation, run the following loop

for ( i in 1:length(negativedistinct$ID)){
  
  ## get gene code
  gene1<-negativedistinct[i,"gene1"]
  gene2<-negativedistinct[i,"gene2"]
  
  
  ## get eQTL (cis)
  snpgene1<-cistransdata[cistransdata$gene==gene1,"SNP"]
  snpgene2<-cistransdata[cistransdata$gene==gene2,"SNP"]
  
  ## get eQTL (trans)
  transgene1<-transdata[transdata$V1==gene1,]
  transgene2<-transdata[transdata$V1==gene2,]
  
  ## intersect trans eQTL
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
  
  ## intersect cis eQTL
  if(length(cisgene1$V1)!=0 && length(cisgene2$V1)!=0){
    intersectcis<-intersect(cisgene1$V8,snpgene2)
    if(length(intersecttrans)>=0){
      
    }else{
      intersectcis<-intersect(cisgene2$V8,snpgene1)
    }
  }else{
    intersectcis<-NULL
  }
  
  ### if any eQTL is shared, then save shared status
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


## save data
write.table(negativedistinct,file="negative_distinct_withindepth_02.txt",row.names = FALSE)