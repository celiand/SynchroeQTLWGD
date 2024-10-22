### annotation with ncrna and TF

##Tf part
library(tidyverse)
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



## NCRNA part - use info from the gff file
ncrnatab<-read.table(file="/mnt/SCRATCH/cedi/phDSalmon/evolutionary_analysis/dnds/non_coding_rna_noheader.gff",sep="\t")
ncrnatab<-ncrnatab[,-c(6:9)]
colnames(ncrnatab)<-c("chrom","ensembl","type","start","end")

mirnatab<-ncrnatab[ncrnatab$type=="miRNA",]
lncRNAtab<-ncrnatab[ncrnatab$type=="lnc_RNA",]
ncRNAtab<-ncrnatab[ncrnatab$type=="ncRNA",]