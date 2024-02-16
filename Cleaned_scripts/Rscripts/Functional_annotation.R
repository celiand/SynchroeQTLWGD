### Celian Diblasi
### Functional annotation of eQTL in non coding RNA and transcription factors

### get TF and non coding RNA data


library(tidyverse)

## TFs
Eluc_TFs <- read_tsv("http://bioinfo.life.hust.edu.cn/AnimalTFDB4/static/download/TF_list_final/Esox_lucius_TF")
OGtbl <- read_tsv("https://salmobase.org/datafiles/TSV/og/2021-11.tsv")

Ssal_TFs <-
  Eluc_TFs %>% 
  left_join(select(filter(OGtbl,spc=="Eluc"),gene_id,teleost), by=c("Ensembl"="gene_id")) %>% 
  inner_join(select(filter(OGtbl,spc=="Ssal"),Ssal_geneID=gene_id,teleost), by="teleost") %>% 
  distinct(Ssal_geneID, Symbol, Family)


colnames(Ssal_TFs)[1]<-"gene_id"

tablegene<-read.table(file="https://salmobase.org/datafiles/TSV/genes/AtlanticSalmon/Ssal_v3.1/Ensembl_genes.tsv",header=TRUE,sep="\t")
tablegene<-tablegene[!is.na(as.numeric(as.character(tablegene$seqname))),]


mergetable<-left_join(tablegene,Ssal_TFs,by=c("gene_id"))



##non coding RNAs
ncrnatab<-read.table(file="non_coding_rna_noheader.gff",sep="\t") ##gff3 data only containing ncRNAs info
ncrnatab<-ncrnatab[,-c(6:9)]
colnames(ncrnatab)<-c("chrom","ensembl","type","start","end")

mirnatab<-ncrnatab[ncrnatab$type=="miRNA",]
lncRNAtab<-ncrnatab[ncrnatab$type=="lnc_RNA",]
ncRNAtab<-ncrnatab[ncrnatab$type=="ncRNA",]

#### Annotation of eQTLs

rna_full_table<-read.table(file="rna_LL_all_august_version_01.txt",header=TRUE) ### data come from XXX

## set up a distance around the location of the annotation to consider the eQTL as "nearby this annotation"
distTF<-20000
distncrna<-20000

##create an empty table
tabannotation<-c()
for(i in 1:length(rna_full_table$kind)){ ##run a loop for all eQTL
  annotations<-c()
  chreqtl<-rna_full_table[i,"snpchrom"]
  poseqtl<-rna_full_table[i,"snppos"]
  posminusTFeqtl<-poseqtl-distTF
  posplusTFeqtl<-poseqtl+distTF
  posminusRNA<-poseqtl-distncrna
  posplusRNA<-poseqtl+distncrna
  
  subtabmirna<-mirnatab[mirnatab$chrom==chreqtl & ((mirnatab$start>=posminusRNA & mirnatab$start<=posplusRNA) | (mirnatab$end>=posminusRNA & mirnatab$end<=posplusRNA)),]
  if(length(subtabmirna$end>0)){ ##if a mirna is found nearby
    annotations<-c(annotations,"miRNA")
  }
  
  sublncRNAtab<-lncRNAtab[lncRNAtab$chrom==chreqtl & ((lncRNAtab$start>=posminusRNA & lncRNAtab$start<=posplusRNA) | (lncRNAtab$end>=posminusRNA & lncRNAtab$end<=posplusRNA)),]
  if(length(sublncRNAtab$end>0)){ ##if a lncRNA is found nearby
    annotations<-c(annotations,"lncRNA")
  }
  
  subncRNAtab<-ncRNAtab[ncRNAtab$chrom==chreqtl & ((ncRNAtab$start>=posminusRNA & ncRNAtab$start<=posplusRNA) | (ncRNAtab$end>=posminusRNA & ncRNAtab$end<=posplusRNA)),]
  if(length(subncRNAtab$end>0)){ ##if a lncRNA is found nearby
    annotations<-c(annotations,"ncRNA")
  }
  
  subtabtf<-mergetable[mergetable$seqname==chreqtl & ((mergetable$start>=posminusTFeqtl & mergetable$start<=posplusTFeqtl) | (mergetable$end>=posminusTFeqtl & mergetable$end<=posplusTFeqtl)),]
  subtabtf2<-subtabtf[!is.na(subtabtf$Family),]
  if(length(subtabtf2$end>0)){ ##if a TF is found nearby
    annotations<-c(annotations,"TF")
  }
  
  if(is.null(annotations)){
    annotation<-"None"
  }else{
    annotation <- paste(annotations, collapse = "_")
  }
  
  ##get number of cis gene regulated
  snp<-rna_full_table[i,"lead_SNP"]
  cistab<-rna_full_table[rna_full_table$lead_SNP==snp & rna_full_table$kind=="cis","kind"]
  transtab<-rna_full_table[rna_full_table$lead_SNP==snp & rna_full_table$kind=="trans","kind"]
  nbcis<-length(cistab)
  nbtrans<-length(transtab)
  if(nbcis==0 & nbtrans>0){
    status<-"trans"
  }else if(nbcis>0 & nbtrans==0){
    status<-"cis"
  }else{
    status<-"both"
  }
  total<-nbcis+nbtrans
  
  row<-c(chreqtl,poseqtl,rna_full_table[i,"gene"],annotation,nbcis,nbtrans,total,status)
  tabannotation<-rbind.data.frame(tabannotation,row,stringsAsFactors = FALSE)
}
##format table and export
colnames(tabannotation)<-c("CHROM","POS","gene","annotation","nbcis","nbtrans","nbgene","status")
dups <- duplicated(subset(tabannotation, select = c("CHROM", "POS")))
tabannotation<-tabannotation[!dups,]

write.table(tabannotation,file="tabannotation05.txt")


#### Annotation of SNPs that are not eQTLs

## set up a distance around the location of the annotation to consider the SNP as "nearby this annotation"
distTF<-20000
distncrna<-20000

allsnps<-read.table(file="LL_0.01.recode.vcf.gz_1100000000.stat",header=TRUE) ## get all SNPs data
##fornat table
allsnps$chromnb<-ifelse(substr(allsnps$CHROM,4,4)==0,substr(allsnps$CHROM,5,5),substr(allsnps$CHROM,4,5))
allsnps$snpid<-paste(allsnps$chromnb,allsnps$POS,sep="_")

## remove SNPs that are eQTLs
anottable<-read.table(file="tabannotation05.txt",header=TRUE)
anottable$snpid<-paste(anottable$CHROM,anottable$POS,sep="_")

allsnpsnoeqtl<-allsnps[!(allsnps$snpid%in%anottable$snpid),]

##make a matrix to store the results
tabannotation<-matrix(nrow=354381,ncol=4)

for(i in 1:length(allsnpsnoeqtl$POS)){ ##run a loop for each SNPs
  annotations<-c()
  chreqtl<-allsnpsnoeqtl[i,"chromnb"]
  poseqtl<-allsnpsnoeqtl[i,"POS"]
  posminusTFeqtl<-poseqtl-distTF
  posplusTFeqtl<-poseqtl+distTF
  posminusRNA<-poseqtl-distncrna
  posplusRNA<-poseqtl+distncrna
  
  subtabmirna<-mirnatab[mirnatab$chrom==chreqtl & ((mirnatab$start>=posminusRNA & mirnatab$start<=posplusRNA) | (mirnatab$end>=posminusRNA & mirnatab$end<=posplusRNA)),]
  if(length(subtabmirna$end>0)){ ##if a mirna is found nearby
    annotations<-c(annotations,"miRNA")
  }
  
  sublncRNAtab<-lncRNAtab[lncRNAtab$chrom==chreqtl & ((lncRNAtab$start>=posminusRNA & lncRNAtab$start<=posplusRNA) | (lncRNAtab$end>=posminusRNA & lncRNAtab$end<=posplusRNA)),]
  if(length(sublncRNAtab$end>0)){ ##if a lncRNA is found nearby
    annotations<-c(annotations,"lncRNA")
  }
  
  subncRNAtab<-ncRNAtab[ncRNAtab$chrom==chreqtl & ((ncRNAtab$start>=posminusRNA & ncRNAtab$start<=posplusRNA) | (ncRNAtab$end>=posminusRNA & ncRNAtab$end<=posplusRNA)),]
  if(length(subncRNAtab$end>0)){ ##if a lncRNA is found nearby
    annotations<-c(annotations,"ncRNA")
  }
  
  subtabtf<-mergetable[mergetable$seqname==chreqtl & ((mergetable$start>=posminusTFeqtl & mergetable$start<=posplusTFeqtl) | (mergetable$end>=posminusTFeqtl & mergetable$end<=posplusTFeqtl)),]
  subtabtf2<-subtabtf[!is.na(subtabtf$Family),]
  if(length(subtabtf2$end>0)){ ##if a TF is found nearby
    annotations<-c(annotations,"TF")
  }
  
  if(is.null(annotations)){
    annotation<-"None"
  }else{
    annotation <- paste(annotations, collapse = "_")
  }
  row<-c(chreqtl,poseqtl,"none",annotation)
  #tabannotation<-rbind.data.frame(tabannotation,row,stringsAsFactors = FALSE)
  tabannotation[i,]<-row
}
##format table and export
tabannotation<-as.data.frame(tabannotation)
colnames(tabannotation)<-c("CHROM","POS","gene","annotation")
dups <- duplicated(subset(tabannotation, select = c("CHROM", "POS")))
tabannotation<-tabannotation[!dups,]

write.table(tabannotation,file="tabannotation_noeqtl02.txt")




### analyze the result (lead SNPs)

eqtltab<-read.table("tabannotation05.txt",header=TRUE) ### import eQTL overlapping elements
noneqtltab<-read.table("tabannotation_noeqtl02.txt",header=TRUE) ### import SNP with no connection overlapping elements
noneqtltab$nbcis<-0
noneqtltab$nbtrans<-0
noneqtltab$nbgene<-0
noneqtltab$status<-"none"

## merge everything together in one table
fulltab<-rbind.data.frame(eqtltab,noneqtltab,stringsAsFactors = FALSE)
fulltab<-fulltab[,-3]
write.table(fulltab,"full_eqtltab_02.txt",row.names = FALSE)


## import the table and get only annotation of interest
## because some might have several annotation
## do unique(tabeqtl$annotation) or table(tabeqtl$annotation) to check
tabeqtl<-read.table("full_eqtltab_02.txt",header=TRUE)
subtab<-tabeqtl[tabeqtl$annotation%in%c("None","lncRNA","miRNA","TF"),]


##rename annotation for later plots
subtab[subtab$status=="both","status"]<-"Cis&trans regulator"
subtab[subtab$status=="cis","status"]<-"Cis regulator"
subtab[subtab$status=="trans","status"]<-"Trans regulator"
subtab[subtab$status=="none","status"]<-"SNP with no regulatory connections"


subtab[subtab$annotation=="TF","annotation"]<-"Transcription factor"
subtab[subtab$annotation=="None","annotation"]<-"Other"

library(ggplot2)
library(tidyverse)

dfcount <- subtab %>% group_by(annotation,status) %>% 
  summarise(total_count=n(),.groups = 'drop') %>%
  as.data.frame()

#ggplot(dfcount,aes(x=status,y=total_count,fill=annotation))+geom_bar(stat="identity",position="fill")

#ggplot(subtab[subtab$status!="none",],aes(x=annotation,y=nbgene))+geom_boxplot()+ylim(0,15) ##majority regulate 1 gene whatever the annotation

dfprop<-dfcount %>% group_by(status) %>% mutate(proportion = total_count / sum(total_count))

dfprop <- dfprop %>%
  group_by(annotation) %>%
  mutate(prop_division = proportion / proportion[status == "SNP with no regulatory connections"])

colnames(dfprop)[2]<-c("Connection type")
#dfprop$`Connection type`<-factor(dfprop$`Connection type`,levels=c("SNP with no regulatory connections","Cis regulator","Cis&trans regulator","Trans regulator"))
dfprop<-dfprop %>% arrange(desc(total_count))

ggplot(dfprop,aes(x=annotation,y=prop_division,color=`Connection type`,shape=`Connection type`))+geom_hline(yintercept=1,color="lightgrey",linetype = "dashed",size=1.2)+geom_point(size=log(dfprop$total_count))+theme_bw(25)+scale_color_manual(values=c("firebrick1","dodgerblue3","darkorchid3","darkgoldenrod1"))+scale_shape_manual(values=c(15,16,17,18))+
  ylab("Enrichment")+xlab("Regulatory element type")+ guides(colour = guide_legend(override.aes = list(size=5)))+theme(axis.text.x = element_text(angle = 45, hjust=1)) 



##remove SNP with no connections from the plots
dfprop2<-dfprop[dfprop$`Connection type`!="SNP with no regulatory connections",]

ggplot(dfprop2,aes(x=annotation,y=prop_division,color=`Connection type`,shape=`Connection type`))+geom_hline(yintercept=1,color="lightgrey",linetype = "dashed",size=2.5)+geom_point(size=log(dfprop2$total_count))+theme_bw(25)+scale_color_manual(values=c("firebrick1","dodgerblue3","darkgoldenrod1"))+scale_shape_manual(values=c(15,16,18))+
  ylab("Enrichment")+xlab("Regulatory element type")+ guides(colour = guide_legend(override.aes = list(size=5)))+theme(axis.text.x = element_text(angle = 45, hjust=1)) 







#### make a chi test to see significance
##make a table
contingency_table<-table(subtab$annotation, subtab$status)
chisq.test(contingency_table)

chi_square_test <- chisq.test(contingency_table)
chi_square_test$residuals