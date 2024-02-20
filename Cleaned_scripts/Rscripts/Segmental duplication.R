### Celian Diblasi
### 19/02/2024
### investigation of segmental duplication


### -- Find segmental duplication -- ###

### get orthology data from Atlantic salmon 
ortho_tab<-read.table(file="https://salmobase.org/datafiles/orthology/2021-11/Ortho_pipeline/OGtbl.tsv",header=TRUE)

filteredtab<-ortho_tab[ortho_tab$spc=="Ssal",1:2]

filteredtab2<-filteredtab
for( i in unique(filteredtab2$OG)){
  lengthog<-length(filteredtab2[filteredtab2$OG==i,"OG"])
  ## removing OG with more than 10 genes
  if(lengthog>10){
    filteredtab2<-filteredtab2[!(filteredtab2$OG==i),]
  }
}

### get ensembl genes
tablegene<-read.table(file="https://salmobase.org/datafiles/TSV/genes/AtlanticSalmon/Ssal_v3.1/Ensembl_genes.tsv",header=TRUE,sep="\t")
tablegene<-tablegene[!is.na(as.numeric(as.character(tablegene$seqname))),]
tablegene$CHROM<-paste(ifelse(as.numeric(as.character(tablegene$seqname))>9,"ssa","ssa0"),tablegene$seqname,sep="")

mergeddata<-merge(filteredtab2,tablegene,by.x="geneID",by.y="gene_id")


## 14109 genes initially (## 16080 genes initially if not removing OG>10)

library(tidyverse)

mergeddata$OG_chrom<-paste(mergeddata$OG,mergeddata$seqname,sep="_")

mergeddata2<-mergeddata
## filter OG with only 1 gene on the same chromosome
for (i in mergeddata2$OG_chrom){
  lengthOG2<-length(mergeddata2[mergeddata2$OG_chrom==i,1])
  if(lengthOG2<2){
    mergeddata2<-mergeddata2[!(mergeddata2$OG_chrom)==i,]
  }
}

## 1089 that are part of OG group with at least 2 genes on the same chromosome




# Function to check for overlaps
check_overlaps <- function(data) {
  data <- data[order(data$start), ]
  overlaps <- logical(nrow(data))
  for (i in 1:(nrow(data) - 1)) {
    overlaps[i] <- any(data$start[i + 1] <= data$end[i] & data$end[i + 1] >= data$start[i])
  }
  return(overlaps)
}

##remove OG if any genes overlap


mergeddata3<-mergeddata2

for( i in unique(mergeddata3$OG_chrom)){
  subdata<-mergeddata3[mergeddata3$OG_chrom==i,]
  overlap<-check_overlaps(subdata)
  if(any(overlap)){
    mergeddata3<-mergeddata3[!(mergeddata3$OG_chrom==i),]
  }
}

## 906 genes within OG with no overlap


dist<-100000 ##define distance between genes to consider them as duplicates

for( i in unique(mergeddata3$OG_chrom)){
  subdata<-mergeddata3[mergeddata3$OG_chrom==i,]
  for(y in 1:length(subdata[,1])){
    vector<- seq(1, length(subdata[,1]))
    
    # Remove the current value form the vector
    vector <- vector[vector != y]
    
    startgene<-subdata[y,"start"]
    endgene<-subdata[y,"end"]
    
    boundariestartmax<-startgene+dist
    boundariestartmin<-startgene-dist
    boundarieendmax<-endgene+dist
    boundarieendmin<-endgene-dist
    
    startvector<-subdata[vector,"start"]
    endvector<-subdata[vector,"end"]
    
    if((any(startvector>boundariestartmin & startvector<boundariestartmax)) | (any(endvector>boundarieendmin & endvector<boundarieendmax))){
      
    }else{
      mergeddata3<-mergeddata3[!(mergeddata3$OG_chrom==i & mergeddata3$start==startgene & mergeddata3$end==endgene),]
      tests<-"oui"
    }
  }
}


## 597 genes left

length(unique(mergeddata3$OG_chrom))
length(unique(mergeddata3$OG))

###all genes that are part of a potential segmental duplicate group
write.table(mergeddata3,file="potentialSD.txt",row.names=FALSE)


#### -- Refine shared eQTL category (find if one of the eQTL regulation gene1 is also regulating gene2 but not lead eQTL) -- ####

transdata<-read.table(file="LLsmolt.trans.adjust.hits.txt",header=FALSE)

cisdata<-read.table(file="LLsmolt.cisnominal.txt",header=FALSE)

## for each OG term, check if two eQTls for two genes are on the same chromosome (likely to be LD with eachothers)
## then do pairs like that (pairID, gene1, gene2, eQTL1, eQTL2)
distincteQTl<-sd_eQTL[sd_eQTL$eQTL_status=="distinct_eQTL",]


pairdata<-c()
for( i in unique(distincteQTl$OG)){
  subdata<-distincteQTl[distincteQTl$OG==i,]
  for( y in 1:length(subdata$geneID)){
    subdata2<-subdata[subdata$snpchrom==subdata[y,"snpchrom"] & subdata$geneID!=subdata[y,"geneID"],]
    if(length(subdata2$geneID)>0){
      for(z in 1:length(subdata2$geneID)){
        ### test if the pair is not already registered
        existdata<-pairdata[pairdata$GeneID1==subdata2[z,"geneID"] & pairdata$GeneID2==subdata[y,"geneID"] & pairdata$Lead_SNP1==subdata2[z,"lead_SNP"] & pairdata$Lead_SNP2==subdata[y,"lead_SNP"],1]
        if(length(existdata)<1){
          rowpair<-data.frame(GeneID1=subdata[y,"geneID"],Lead_SNP1=subdata[y,"lead_SNP"],snpchrom1=subdata[y,"snpchrom"],snppos1=subdata[y,"snppos"],kind1=subdata[y,"kind"],GeneID2=subdata2[z,"geneID"],Lead_SNP2=subdata2[z,"lead_SNP"],snpchrom2=subdata2[z,"snpchrom"],snppos2=subdata2[z,"snppos"],kind2=subdata2[z,"kind"],OG=i)
          pairdata<-rbind.data.frame(pairdata,rowpair,stringsAsFactors = FALSE)
        }
      }
    }
  }
}





for ( i in 1:length(pairdata$GeneID1)){
  gene1<-pairdata[i,"GeneID1"]
  gene2<-pairdata[i,"GeneID2"]
  
  snpgene1<-pairdata[i,"Lead_SNP1"]
  snpgene2<-pairdata[i,"Lead_SNP2"]
  
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
  
  pairdata[i,"indepth_test"]<-val
  pairdata[i,"commonsnp"]<-sampled_snp
}


write.table(pairdata,file="refined_shared_eQTL_SD.txt",row.names = FALSE)





### number of genes and pair at each steps, and number of pairs with shared eQTL


## -- Summary -- ##
##some pairs have DC patterns
## 196 gene with distinct eQTL --> 87 pairs with eQTL on the same chromosome --> 50 with shared eQTL




### makes numbers clear:

##starting point:
length(unique(sd_eQTL$geneID))  
## 116 unique genes


### make all potential pairs
allpairs<-c()
for( i in unique(sd_eQTL$OG)){
  subgenes<-unique(sd_eQTL[sd_eQTL$OG==i,"geneID"])
  
  ##make all possible combinations
  if(length(subgenes)>1){
    
    all_combinations <- combn(subgenes, 2)
    for (y in 1:length(data.frame(all_combinations))){
      pairID<-paste(all_combinations[,y],collapse = "_")
      allpairs<-c(allpairs,pairID)
    }
  }
}

## 52 unique pairs of genes

length(unique(sharedsd$geneID))

## in shared eQTL, 4 pairs (8 genes)

length(unique(data$pair))

length(unique(c(data$gene1,data$gene2)))

## 23 refined pairs, 45 genes

## total: 27 pairs out of 52 possible with shared eQTL; 49 genes out of 116

### plot
library(ggplot2)

data_SD<-data.frame(Data=c("Initial","Nb Genes in OG < 10","At least 2 genes on the same chr","No overlap","Dist between genes<100 000bp","With eQTL regulation"),Nb_of_genes=c(16080,14109,1089,906,597,116))


data_SD$Data<-factor(data_SD$Data,levels=c("Initial","Nb Genes in OG < 10","At least 2 genes on the same chr","No overlap","Dist between genes<100 000bp","With eQTL regulation"))

ggplot(data_SD,aes(x=Data,y=Nb_of_genes,fill=Data))+geom_bar(stat="identity")+theme_bw(20)+theme(axis.text.x = element_text(angle = 45, hjust=1))+scale_y_log10()