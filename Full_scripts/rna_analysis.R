duplication_table<-read.table(file="https://salmobase.org/datafiles/TSV/synteny/2021-11/AtlanticSalmon/synteny.tsv",header=TRUE,sep="\t")

rna_table<-read.table(file="Synchro_eQTL_5Mwindow_uniq.csv",sep=",",header=TRUE)
rna_table$snpschrom<-sapply(strsplit(rna_table$SNP,"_"),FUN = `[[`, 1)
rna_table$snpschrom<-ifelse(substr(rna_table$snpschrom,4,4)==0,substr(rna_table$snpschrom,5,5),substr(rna_table$snpschrom,4,5))
rna_table$snpstart<-as.numeric(as.character(sapply(strsplit(rna_table$SNP,"_"),FUN = `[[`, 2)))




##debug
#rownm<-1
#i<-1
#(rna_table[rownm,"snpstart"]<=duplication_table[i,"end_x"] && rna_table[rownm,"snpstart"]>=duplication_table[i,"begin_x"] && rna_table[rownm,"snpschrom"]==duplication_table[i,"chromosome_x"])

TestDuplicate<- function(rownm){
  var<-"non"
  for (i in 1:length(duplication_table$end_x)){
    if((rna_table[rownm,"snpstart"]<=duplication_table[i,"end_x"] && rna_table[rownm,"snpstart"]>=duplication_table[i,"begin_x"] && rna_table[rownm,"snpschrom"]==duplication_table[i,"chromosome_x"]) | (rna_table[rownm,"snpstart"]<=duplication_table[i,"end_y"] && rna_table[rownm,"snpstart"]>=duplication_table[i,"begin_y"] && rna_table[rownm,"snpschrom"]==duplication_table[i,"chromosome_y"]) ){
      var<-"oui"
    }else{
    }
  }
  if(var=="oui"){
    return(TRUE)
  }else{
    return(FALSE)
  }
}

for(y in 1:length(rna_table$snpschrom)){
  if(TestDuplicate(y)==TRUE){
    rna_table[y,"indup"]<-"Indup"
  }else{
    rna_table[y,"indup"]<-"Notindup"
  }
  
}


write.table(rna_table,file="Synchro_uniq_withdupstatus.txt",row.names=FALSE)

barplot(prop.table(table(rna_table[rna_table$cistrans=="trans","indup"])))
barplot(prop.table(table(rna_table[rna_table$cistrans=="cis","indup"])))

#table(rna_table[rna_table$cistrans=="trans","indup"])
library(ggplot2)


#https://stackoverflow.com/questions/3695497/show-percent-instead-of-counts-in-charts-of-categorical-variables

ggplot(rna_table,aes(x=indup))+geom_bar()

##if i want to do with dots...

counttot<-c(table(rna_table$indup)[1],table(rna_table$indup)[2])
counttrans<-c(table(rna_table[rna_table$cistrans=="trans","indup"])[1],table(rna_table[rna_table$cistrans=="trans","indup"])[2])
countcis<-c(table(rna_table[rna_table$cistrans=="cis","indup"])[1],table(rna_table[rna_table$cistrans=="cis","indup"])[2])

totvec<-c(counttot/sum(counttot))
transvec<-c(counttrans/sum(counttrans))
cisvec<-c(countcis/sum(countcis))
vecglob<-c(totvec,transvec,cisvec)
dupstat<-c("Indup","Outdup","Indup","Outdup","Indup","Outdup")
location<-c("TOT","TOT","TRANS","TRANS","CIS","CIS")

tableglobal<-cbind.data.frame(vecglob,dupstat,location)


ggplot(tableglobal,aes(x=dupstat,y=vecglob,color=location))+geom_point()+scale_color_manual(values=c("firebrick1","black","dodgerblue2"))+ylim(0,1)+ggtitle("Proportion of SNPs in or out of a duplication region")


## if i want to have a chromosome by chromosome ratio
tabtot<-c()
for (i in 1:29){
  valcounttrans<-table(rna_table[rna_table$cistrans=="trans" & rna_table$snpschrom==i,"indup"])
  if(length(valcounttrans)==2){
    ratiotrans<-valcounttrans[1]/(valcounttrans[1]+valcounttrans[2])
    vecteurtrans<-c(i,ratiotrans,"trans")
  }else{
    if(names(valcounttrans[1])=="Indup"){
      vecteurtrans<-c(i,1,"trans")
    }else if(names(valcounttrans[1])=="Notindup"){
      vecteurtrans<-c(i,0,"trans")
    }
  }
  valcountcis<-table(rna_table[rna_table$cistrans=="cis" & rna_table$snpschrom==i,"indup"])
  if(length(valcountcis)==2){
    ratiocis<-valcountcis[1]/(valcountcis[1]+valcountcis[2])
    vecteurcis<-c(i,ratiocis,"cis")
  }else{
    if(names(valcountcis[1])=="Indup"){
      vecteurcis<-c(i,1,"cis")
    }else if(names(valcountcis[1])=="Notindup"){
      vecteurcis<-c(i,0,"cis")
    }
    
  }
  valcounttot<-table(allsnpdata[allsnpdata$CHROM==i,"indup"])
  if(length(valcounttot)==2){
    ratiotot<-valcounttot[1]/(valcounttot[1]+valcounttot[2])
    vecteurtot<-c(i,ratiotot,"tot")
  }else{
    if(names(valcounttot[1])=="Indup"){
      vecteurtot<-c(i,1,"tot")
    }else if(names(valcounttot[1])=="Notindup"){
      vecteurtot<-c(i,0,"tot")
    }
    
  }
  tabtot<-rbind.data.frame(tabtot,vecteurtrans,vecteurcis,vecteurtot)
}
colnames(tabtot)<-c("chrom","ratio","effect")
tabtot$chrom<-as.numeric(as.character(tabtot$chrom))
tabtot$ratio<-as.numeric(as.character(tabtot$ratio))


ggplot(tabtot,aes(x=chrom,y=ratio,color=effect))+geom_point()+ylim(0,1)+scale_color_manual(values=c("firebrick1","black","dodgerblue2"))

##if we had all snp data

allsnpdata<-read.table(file="imputed_filtered_smolts.vcf.gz_1110000000.stat",header=TRUE)
allsnpdata$CHROM<-ifelse(substr(allsnpdata$CHROM,4,4)==0,substr(allsnpdata$CHROM,5,5),substr(allsnpdata$CHROM,4,5))




TestDuplicate2<- function(rownm){
  var<-"non"
  for (i in 1:length(duplication_table$end_x)){
    if((allsnpdata[rownm,"POS"]<=duplication_table[i,"end_x"] && allsnpdata[rownm,"POS"]>=duplication_table[i,"begin_x"] && allsnpdata[rownm,"CHROM"]==duplication_table[i,"chromosome_x"]) | (allsnpdata[rownm,"POS"]<=duplication_table[i,"end_y"] && allsnpdata[rownm,"POS"]>=duplication_table[i,"begin_y"] && allsnpdata[rownm,"CHROM"]==duplication_table[i,"chromosome_y"]) ){
      var<-"oui"
    }else{
    }
  }
  if(var=="oui"){
    return(TRUE)
  }else{
    return(FALSE)
  }
}

for(y in 1:length(allsnpdata$CHROM)){
  if(TestDuplicate2(y)==TRUE){
    allsnpdata[y,"indup"]<-"Indup"
  }else{
    allsnpdata[y,"indup"]<-"Notindup"
  }
  
}

write.table(allsnpdata,file="Synchro_all_withdupstatus.txt",row.names=FALSE)

##test chi test
tabchi<-rbind.data.frame(countcis,counttrans)
colnames(tabchi)<-c("indup","outdup")
row.names(tabchi)<-c("cis","trans")

chisq.test(tabchi)

##ggplot(rna_table,aes(x=indup))+geom_bar()
##not working as expected... but i learned
##ggplot(rna_table,aes(x=indup,color=cistrans))+geom_point(aes(y=..count../sum(..count..)), stat="count")




## with the other data (not uniq)

rna_table<-read.table(file="Synchro_eQTL_5Mwindow.csv",sep=",",header=TRUE)
rna_table$snpschrom<-sapply(strsplit(rna_table$SNP,"_"),FUN = `[[`, 1)
rna_table$snpschrom<-ifelse(substr(rna_table$snpschrom,4,4)==0,substr(rna_table$snpschrom,5,5),substr(rna_table$snpschrom,4,5))
rna_table$snpstart<-as.numeric(as.character(sapply(strsplit(rna_table$SNP,"_"),FUN = `[[`, 2)))


rna_table_trans<-rna_table[rna_table$cistrans=="trans",]
## things to check for, for each SNP - gene pair (in trans only ?):
## if gene is in duplicated region but not SNP
## reverse also
## if both are in dupicated region
## for the above case, distinguish cases where they are in a WGD pair or not

##if a part only of the chromosome are written as ssaxx
#rna_table_trans$gene_chr<-ifelse(substr(rna_table_trans$gene_chr,1,1)=="s",ifelse(substr(rna_table_trans$gene_chr,4,4)==0,substr(rna_table_trans$gene_chr,5,5),substr(rna_table_trans$gene_chr,4,5)),rna_table_trans$gene_chr)


##add a test row
#rna_table_trans[32638,]<-c("test","ENSSSAG00000116091",11,4901297,"ssa26_5654562",0.01,NA,"trans",26,5654562)
#rna_table_trans<-rna_table_trans[30000:32638,]


rna_table_trans$snpstart<-as.numeric(as.character(rna_table_trans$snpstart))
rna_table_trans$gene_start<-as.numeric(as.character(rna_table_trans$gene_start))


##with pvalue
rna_table_trans<-rna_table_trans[rna_table_trans$pval<=0.01,]
rna_table_trans<-rna_table_trans[rna_table_trans$pval<=0.0001,]
rna_table_trans<-rna_table_trans[rna_table_trans$pval<=0.0000000001,]

##with ligh conditions


TestDuplicate4<- function(rownm){
  var<-"no"
  var1<-"no"
  var2<-"no"
  for (i in 1:length(duplication_table$end_x)){
    if((rna_table_trans[rownm,"snpstart"]<=duplication_table[i,"end_x"] && rna_table_trans[rownm,"snpstart"]>=duplication_table[i,"begin_x"] && rna_table_trans[rownm,"snpschrom"]==duplication_table[i,"chromosome_x"])){
      if((rna_table_trans[rownm,"gene_start"]<=duplication_table[i,"end_y"] && rna_table_trans[rownm,"gene_start"]>=duplication_table[i,"begin_y"] && rna_table_trans[rownm,"gene_chr"]==duplication_table[i,"chromosome_y"])){
        var <- "dup_pair"
        break
      }else if((rna_table_trans[rownm,"gene_start"]<=duplication_table[i,"end_x"] && rna_table_trans[rownm,"gene_start"]>=duplication_table[i,"begin_x"] && rna_table_trans[rownm,"gene_chr"]==duplication_table[i,"chromosome_x"])){
        var <- "dup_same_chr"
        break
      }else{
        var1 <- "SNP_dup"
      }
    }else if((rna_table_trans[rownm,"snpstart"]<=duplication_table[i,"end_y"] && rna_table_trans[rownm,"snpstart"]>=duplication_table[i,"begin_y"] && rna_table_trans[rownm,"snpschrom"]==duplication_table[i,"chromosome_y"])){
      if((rna_table_trans[rownm,"gene_start"]<=duplication_table[i,"end_x"] && rna_table_trans[rownm,"gene_start"]>=duplication_table[i,"begin_x"] && rna_table_trans[rownm,"gene_chr"]==duplication_table[i,"chromosome_x"])){
        var <- "dup_pair"
        break
      }else if((rna_table_trans[rownm,"gene_start"]<=duplication_table[i,"end_y"] && rna_table_trans[rownm,"gene_start"]>=duplication_table[i,"begin_y"] && rna_table_trans[rownm,"gene_chr"]==duplication_table[i,"chromosome_y"])){
        var <- "dup_same_chr"
        break
      }else{
        var1 <- "SNP_dup"
      }
    }else if((rna_table_trans[rownm,"gene_start"]<=duplication_table[i,"end_y"] && rna_table_trans[rownm,"gene_start"]>=duplication_table[i,"begin_y"] && rna_table_trans[rownm,"gene_chr"]==duplication_table[i,"chromosome_y"]) | (rna_table_trans[rownm,"gene_start"]<=duplication_table[i,"end_x"] && rna_table_trans[rownm,"gene_start"]>=duplication_table[i,"begin_x"] && rna_table_trans[rownm,"gene_chr"]==duplication_table[i,"chromosome_x"])){
      var2 <- "gene_dup"
    }
  }
  if(var1 == "SNP_dup" && var2 == "gene_dup"){
    var <- "dup_no_pair"
  }else if(var1 == "SNP_dup" && var2 != "gene_dup"){
    var <- "SNP_dup"
  }else if(var1 != "SNP_dup" && var2 == "gene_dup"){
    var <- "gene_dup"
  }else if(var=="no"){
    var <- "nodup"
  }
  return(var)
}

for(y in 1:length(rna_table_trans$snpschrom)){
  rna_table_trans[y,"dupstatus"]<-TestDuplicate4(y)
}

##debug
#i=1
#rownm=2639
#(rna_table_trans[rownm,"snpstart"]<=duplication_table[i,"end_y"] && rna_table_trans[rownm,"snpstart"]>=duplication_table[i,"begin_y"] && rna_table_trans[rownm,"snpschrom"]==duplication_table[i,"chromosome_y"])
#(rna_table_trans[rownm,"gene_start"]<=duplication_table[i,"end_x"] && rna_table_trans[rownm,"gene_start"]>=duplication_table[i,"begin_x"] && rna_table_trans[rownm,"gene_chr"]==duplication_table[i,"chromosome_x"])



barplot(prop.table(table(rna_table_trans$dupstatus)))


write.table(rna_table_trans,file="Synchro_uniq_trans_withdupstatus_detailled.txt",row.names=FALSE)
write.table(rna_table_trans,file="Synchro_uniq_trans_withdupstatus_detailled_p001.txt",row.names=FALSE)
write.table(rna_table_trans,file="Synchro_uniq_trans_withdupstatus_detailled_p00001.txt",row.names=FALSE)
write.table(rna_table_trans,file="Synchro_uniq_trans_withdupstatus_detailled_p00000000001.txt",row.names=FALSE)






ggplot(rna_table_trans,aes(x=log10(pval)))+geom_histogram(binwidth = 2)

#rna_table_trans[rna_table_trans$dupstatus=="dup_pair",]

## shuffle combination of trans snp and gene

## substract table to have only columns wanted
rna_table_trans_sub<-rna_table_trans[,c(3,4,9,10)]


## function to shuffle the rows in a column
shuffle_data<- function(table){
  cols_to_shuffle <- c("gene_chr","gene_start")
  cols_to_shuffle2<-c("snpschrom","snpstart")
  
  # combine the columns you want to shuffle into one matrix
  shuffle_matrix <- table[, cols_to_shuffle]
  shuffle_matrix2 <- table[, cols_to_shuffle2]
  
  
  # shuffle the rows of the matrix
  shuffle_matrix <- shuffle_matrix[sample(nrow(shuffle_matrix)),]
  shuffle_matrix2 <- shuffle_matrix2[sample(nrow(shuffle_matrix2)),]
  
  # split the shuffled matrix back into the original columns
  newtab<- cbind.data.frame(shuffle_matrix, shuffle_matrix2)
  return(newtab)
  
}
rna_table_trans_suffled<-shuffle_data(rna_table_trans_sub)



TestDuplicate5<- function(rownm){
  var<-"no"
  var1<-"no"
  var2<-"no"
  for (i in 1:length(duplication_table$end_x)){
    if((rna_table_trans_suffled[rownm,"snpstart"]<=duplication_table[i,"end_x"] && rna_table_trans_suffled[rownm,"snpstart"]>=duplication_table[i,"begin_x"] && rna_table_trans_suffled[rownm,"snpschrom"]==duplication_table[i,"chromosome_x"])){
      if((rna_table_trans_suffled[rownm,"gene_start"]<=duplication_table[i,"end_y"] && rna_table_trans_suffled[rownm,"gene_start"]>=duplication_table[i,"begin_y"] && rna_table_trans_suffled[rownm,"gene_chr"]==duplication_table[i,"chromosome_y"])){
        var <- "dup_pair"
        break
      }else if((rna_table_trans_suffled[rownm,"gene_start"]<=duplication_table[i,"end_x"] && rna_table_trans_suffled[rownm,"gene_start"]>=duplication_table[i,"begin_x"] && rna_table_trans_suffled[rownm,"gene_chr"]==duplication_table[i,"chromosome_x"])){
        var <- "dup_same_chr"
        break
      }else{
        var1 <- "SNP_dup"
      }
    }else if((rna_table_trans_suffled[rownm,"snpstart"]<=duplication_table[i,"end_y"] && rna_table_trans_suffled[rownm,"snpstart"]>=duplication_table[i,"begin_y"] && rna_table_trans_suffled[rownm,"snpschrom"]==duplication_table[i,"chromosome_y"])){
      if((rna_table_trans_suffled[rownm,"gene_start"]<=duplication_table[i,"end_x"] && rna_table_trans_suffled[rownm,"gene_start"]>=duplication_table[i,"begin_x"] && rna_table_trans_suffled[rownm,"gene_chr"]==duplication_table[i,"chromosome_x"])){
        var <- "dup_pair"
        break
      }else if((rna_table_trans_suffled[rownm,"gene_start"]<=duplication_table[i,"end_y"] && rna_table_trans_suffled[rownm,"gene_start"]>=duplication_table[i,"begin_y"] && rna_table_trans_suffled[rownm,"gene_chr"]==duplication_table[i,"chromosome_y"])){
        var <- "dup_same_chr"
        break
      }else{
        var1 <- "SNP_dup"
      }
    }else if((rna_table_trans_suffled[rownm,"gene_start"]<=duplication_table[i,"end_y"] && rna_table_trans_suffled[rownm,"gene_start"]>=duplication_table[i,"begin_y"] && rna_table_trans_suffled[rownm,"gene_chr"]==duplication_table[i,"chromosome_y"]) | (rna_table_trans_suffled[rownm,"gene_start"]<=duplication_table[i,"end_x"] && rna_table_trans_suffled[rownm,"gene_start"]>=duplication_table[i,"begin_x"] && rna_table_trans_suffled[rownm,"gene_chr"]==duplication_table[i,"chromosome_x"])){
      var2 <- "gene_dup"
    }
  }
  if(var1 == "SNP_dup" && var2 == "gene_dup"){
    var <- "dup_no_pair"
  }else if(var1 == "SNP_dup" && var2 != "gene_dup"){
    var <- "SNP_dup"
  }else if(var1 != "SNP_dup" && var2 == "gene_dup"){
    var <- "gene_dup"
  }else if(var=="no"){
    var <- "nodup"
  }
  return(var)
}

for(y in 1:length(rna_table_trans_suffled$snpschrom)){
  rna_table_trans_suffled[y,"dupstatus"]<-TestDuplicate5(y)
}


write.table(rna_table_trans_suffled,file="Synchro_uniq_trans_withdupstatus_detailled_shuffled.txt",row.names=FALSE)
write.table(rna_table_trans_suffled,file="Synchro_uniq_trans_withdupstatus_detailled_shuffled_p001.txt",row.names=FALSE)
write.table(rna_table_trans_suffled,file="Synchro_uniq_trans_withdupstatus_detailled_shuffled_p00001.txt",row.names=FALSE)
write.table(rna_table_trans_suffled,file="Synchro_uniq_trans_withdupstatus_detailled_shuffled_p00000000001.txt",row.names=FALSE)



barplot(prop.table(table(rna_table_trans_suffled$dupstatus)))


subrnatrans<-rna_table_trans[,c(3,4,9,10,11)]

table_shuffle<-table(rna_table_trans_suffled$dupstatus)
table_suffle<-cbind.data.frame(names(table_shuffle),as.numeric(as.character(table_shuffle)),rep("shuffled",6))
colnames(table_suffle)<-c("dupstatus","values","typedata")
table_real<-table(subrnatrans$dupstatus)
table_real<-cbind.data.frame(names(table_real),as.numeric(as.character(table_real)),rep("realdata",6))
colnames(table_real)<-c("dupstatus","values","typedata")


table_tot<-rbind.data.frame(table_real,table_suffle,stringsAsFactors = FALSE)

library(ggplot2)
ggplot(table_tot,aes(x = typedata,y = values,fill = dupstatus)) +geom_bar(position = "fill", stat = "identity")


## with p_values

rna_table_trans<-read.table(file="Synchro_uniq_trans_withdupstatus_detailled.txt",header=TRUE)
rna_table_transp001<-read.table(file="Synchro_uniq_trans_withdupstatus_detailled_p001.txt",header=TRUE)
rna_table_transp00001<-read.table(file="Synchro_uniq_trans_withdupstatus_detailled_p00001.txt",header=TRUE)
rna_table_transp000000001<-read.table(file="Synchro_uniq_trans_withdupstatus_detailled_p00000000001.txt",header=TRUE)


rna_table_trans$pvalgroup<-"All"
rna_table_transp001$pvalgroup<-"10^-2"
rna_table_transp00001$pvalgroup<-"10^-4"
rna_table_transp000000001$pvalgroup<-"10^-10"

rna_table_trans$typedata<-"realdata"
rna_table_transp001$typedata<-"realdata"
rna_table_transp00001$typedata<-"realdata"
rna_table_transp000000001$typedata<-"realdata"


rna_table_trans_sub<-rna_table_trans[,c(3,4,9,10,11,12,13)]
rna_table_trans_subp001<-rna_table_transp001[,c(3,4,9,10,11,12,13)]
rna_table_trans_subp00001<-rna_table_transp00001[,c(3,4,9,10,11,12,13)]
rna_table_trans_subp000000001<-rna_table_transp000000001[,c(3,4,9,10,11,12,13)]


rna_table_trans_shuffled<-read.table(file="Synchro_uniq_trans_withdupstatus_detailled_shuffled.txt",header=TRUE)
rna_table_trans_shuffled_p001<-read.table(file="Synchro_uniq_trans_withdupstatus_detailled_shuffled_p001.txt",header=TRUE)
rna_table_trans_shuffled_p00001<-read.table(file="Synchro_uniq_trans_withdupstatus_detailled_shuffled_p00001.txt",header=TRUE)
rna_table_trans_shuffled_p00000000001<-read.table(file="Synchro_uniq_trans_withdupstatus_detailled_shuffled_p00000000001.txt",header=TRUE)


rna_table_trans_shuffled$pvalgroup<-"All"
rna_table_trans_shuffled_p001$pvalgroup<-"10^-2"
rna_table_trans_shuffled_p00001$pvalgroup<-"10^-4"
rna_table_trans_shuffled_p00000000001$pvalgroup<-"10^-10"

rna_table_trans_shuffled$typedata<-"shuffled"
rna_table_trans_shuffled_p001$typedata<-"shuffled"
rna_table_trans_shuffled_p00001$typedata<-"shuffled"
rna_table_trans_shuffled_p00000000001$typedata<-"shuffled"



table8complete<-rbind.data.frame(rna_table_trans_sub,rna_table_trans_subp001,rna_table_trans_subp00001,rna_table_trans_subp000000001,rna_table_trans_shuffled,rna_table_trans_shuffled_p001,rna_table_trans_shuffled_p00001,rna_table_trans_shuffled_p00000000001,stringsAsFactors = FALSE)

library(dplyr)

dfcount <- table8complete %>% group_by(pvalgroup,dupstatus,typedata) %>% 
  summarise(total_count=n(),.groups = 'drop') %>%
  as.data.frame()

dfcount$mergeddata<-paste(dfcount$pvalgroup,dfcount$typedata,sep="_")

library(ggplot2)
ggplot(dfcount,aes(x = mergeddata,y = total_count,fill = dupstatus)) +geom_bar(position = "fill", stat = "identity")



## let's check in detail duplicated pair
rna_table_trans_suffled$pair<-paste(rna_table_trans_suffled$gene_chr,rna_table_trans_suffled$snpschrom,sep="_")
table_shuffle_dupair<-table(rna_table_trans_suffled[rna_table_trans_suffled$dupstatus=="dup_pair","pair"])
table_suffle_dupair<-cbind.data.frame(names(table_shuffle_dupair),as.numeric(as.character(table_shuffle_dupair)),rep("shuffled",length(as.numeric(as.character(table_shuffle_dupair)))))
colnames(table_suffle_dupair)<-c("duppair","values","typedata")

subrnatrans$pair<-paste(subrnatrans$gene_chr,subrnatrans$snpschrom,sep="_")
table_real_dupair<-table(subrnatrans[subrnatrans$dupstatus=="dup_pair","pair"])
table_real_dupair<-cbind.data.frame(names(table_real_dupair),as.numeric(as.character(table_real_dupair)),rep("realdata",length(as.numeric(as.character(table_real_dupair)))))
colnames(table_real_dupair)<-c("duppair","values","typedata")

table_tot_pair<-rbind.data.frame(table_real_dupair,table_suffle_dupair,stringsAsFactors = FALSE)

library(ggplot2)
ggplot(table_tot_pair,aes(x = typedata,y = values,fill = duppair)) +geom_bar(position = "fill", stat = "identity")




CreateTables<-function(pvalue,light){
  
  rna_table<-read.table(file="Synchro_eQTL_5Mwindow.csv",sep=",",header=TRUE)
  rna_table$snpschrom<-sapply(strsplit(rna_table$SNP,"_"),FUN = `[[`, 1)
  rna_table$snpschrom<-ifelse(substr(rna_table$snpschrom,4,4)==0,substr(rna_table$snpschrom,5,5),substr(rna_table$snpschrom,4,5))
  rna_table$snpstart<-as.numeric(as.character(sapply(strsplit(rna_table$SNP,"_"),FUN = `[[`, 2)))
  
  
  rna_table_trans<-rna_table[rna_table$cistrans=="trans",]
  
  
  rna_table_trans$snpstart<-as.numeric(as.character(rna_table_trans$snpstart))
  rna_table_trans$gene_start<-as.numeric(as.character(rna_table_trans$gene_start))
  
  rna_table_trans<-rna_table_trans[rna_table_trans$condition==light,]
  rna_table_trans<-rna_table_trans[rna_table_trans$pval<=pvalue,]
  
  for(y in 1:length(rna_table_trans$snpschrom)){
    rna_table_trans[y,"dupstatus"]<-TestDuplicate6(y,rna_table_trans)
  }
  
  ## make shuffled data
  
  rna_table_trans_sub<-rna_table_trans[,c(3,4,9,10)]
  rna_table_trans_suffled<-shuffle_data(rna_table_trans_sub)
  
  for(y in 1:length(rna_table_trans_suffled$snpschrom)){
    rna_table_trans_suffled[y,"dupstatus"]<-TestDuplicate6(y,rna_table_trans_suffled)
  }
  
  ##merge data and register
  rna_table_trans_sub<-rna_table_trans[,c(3,4,9,10,11)]
  
  rna_table_trans_sub$pvalgroup<-pvalue
  rna_table_trans_sub$lightgroup<-light
  rna_table_trans_sub$typedata<-"realdata"
  
  rna_table_trans_suffled$pvalgroup<-pvalue
  rna_table_trans_suffled$lightgroup<-light
  rna_table_trans_suffled$typedata<-"shuffled"
  
  tabdir<-paste("Synchro_uniq_trans_withdupstatus_detailled",pvalue,light,"tot.txt",sep="_")
  tabtot<-rbind.data.frame(rna_table_trans_sub,rna_table_trans_suffled,stringsAsFactors = FALSE)
  write.table(tabtot,file=tabdir,row.names = FALSE)
  #return(tabtot)
}


### create all the tables

CreateTables(pvalue=0.0000000001,light="LL")
CreateTables(pvalue=0.0001,light="LL")
CreateTables(pvalue=0.01,light="LL")
CreateTables(pvalue=1,light="LL")


CreateTables(pvalue=0.0000000001,light="8_16")
CreateTables(pvalue=0.0001,light="8_16")
CreateTables(pvalue=0.01,light="8_16")
CreateTables(pvalue=1,light="8_16")


CreateTables(pvalue=0.0000000001,light="12_12")
CreateTables(pvalue=0.0001,light="12_12")
CreateTables(pvalue=0.01,light="12_12")
CreateTables(pvalue=1,light="12_12")


## debug


## make plot

## LL

llall<-read.table(file = "Synchro_uniq_trans_withdupstatus_detailled_1_LL_tot.txt",header=TRUE)
ll001<-read.table(file = "Synchro_uniq_trans_withdupstatus_detailled_0.01_LL_tot.txt",header=TRUE)
ll00001<-read.table(file = "Synchro_uniq_trans_withdupstatus_detailled_1e-04_LL_tot.txt",header=TRUE)
ll0000000001<-read.table(file = "Synchro_uniq_trans_withdupstatus_detailled_1e-10_LL_tot.txt",header=TRUE)


tablell<-rbind.data.frame(llall,ll001,ll00001,ll0000000001,stringsAsFactors = FALSE)

library(dplyr)

dfcount <- tablell %>% group_by(pvalgroup,dupstatus,typedata,lightgroup) %>% 
  summarise(total_count=n(),.groups = 'drop') %>%
  as.data.frame()

dfcount2<-dfcount[dfcount$dupstatus!="dup_same_chr",]

dfcount3<-dfcount2 %>%
  group_by(pvalgroup, typedata,lightgroup) %>%
  summarise(ratio = total_count[dupstatus == "dup_pair"] / sum(total_count[dupstatus!="dup_pair"]))



library(ggplot2)
ggplot(dfcount3,aes(x = as.factor(pvalgroup),y = ratio,fill = typedata)) +geom_bar(position = "dodge", stat = "identity")+ggtitle("ratio of dup in pair/others dup status, for LL conditions")

## same for 8_16


llall<-read.table(file = "Synchro_uniq_trans_withdupstatus_detailled_1_8_16_tot.txt",header=TRUE)
ll001<-read.table(file = "Synchro_uniq_trans_withdupstatus_detailled_0.01_8_16_tot.txt",header=TRUE)
ll00001<-read.table(file = "Synchro_uniq_trans_withdupstatus_detailled_1e-04_8_16_tot.txt",header=TRUE)
ll0000000001<-read.table(file = "Synchro_uniq_trans_withdupstatus_detailled_1e-10_8_16_tot.txt",header=TRUE)


tablell<-rbind.data.frame(llall,ll001,ll00001,ll0000000001,stringsAsFactors = FALSE)

library(dplyr)

dfcount <- tablell %>% group_by(pvalgroup,dupstatus,typedata,lightgroup) %>% 
  summarise(total_count=n(),.groups = 'drop') %>%
  as.data.frame()

dfcount2<-dfcount[dfcount$dupstatus!="dup_same_chr",]

dfcount3<-dfcount2 %>%
  group_by(pvalgroup, typedata,lightgroup) %>%
  summarise(ratio = total_count[dupstatus == "dup_pair"] / sum(total_count[dupstatus!="dup_pair"]))



library(ggplot2)
ggplot(dfcount3,aes(x = as.factor(pvalgroup),y = ratio,fill = typedata)) +geom_bar(position = "dodge", stat = "identity")+ggtitle("ratio of dup in pair/others dup status, for 8_16 conditions")


## same for 12_12


llall<-read.table(file = "Synchro_uniq_trans_withdupstatus_detailled_1_12_12_tot.txt",header=TRUE)
ll001<-read.table(file = "Synchro_uniq_trans_withdupstatus_detailled_0.01_12_12_tot.txt",header=TRUE)
ll00001<-read.table(file = "Synchro_uniq_trans_withdupstatus_detailled_1e-04_12_12_tot.txt",header=TRUE)
ll0000000001<-read.table(file = "Synchro_uniq_trans_withdupstatus_detailled_1e-10_12_12_tot.txt",header=TRUE)


tablell<-rbind.data.frame(llall,ll001,ll00001,ll0000000001,stringsAsFactors = FALSE)

library(dplyr)

dfcount <- tablell %>% group_by(pvalgroup,dupstatus,typedata,lightgroup) %>% 
  summarise(total_count=n(),.groups = 'drop') %>%
  as.data.frame()

dfcount2<-dfcount[dfcount$dupstatus!="dup_same_chr",]

dfcount3<-dfcount2 %>%
  group_by(pvalgroup, typedata,lightgroup) %>%
  summarise(ratio = total_count[dupstatus == "dup_pair"] / sum(total_count[dupstatus!="dup_pair"]))



library(ggplot2)
ggplot(dfcount3,aes(x = as.factor(pvalgroup),y = ratio,fill = typedata)) +geom_bar(position = "dodge", stat = "identity")+ggtitle("ratio of dup in pair/others dup status, for 12_12 conditions")



##
TestDuplicate6<- function(rownm,table){
  var<-"no"
  var1<-"no"
  var2<-"no"
  for (i in 1:length(duplication_table$end_x)){
    if((table[rownm,"snpstart"]<=duplication_table[i,"end_x"] && table[rownm,"snpstart"]>=duplication_table[i,"begin_x"] && table[rownm,"snpschrom"]==duplication_table[i,"chromosome_x"])){
      if((table[rownm,"gene_start"]<=duplication_table[i,"end_y"] && table[rownm,"gene_start"]>=duplication_table[i,"begin_y"] && table[rownm,"gene_chr"]==duplication_table[i,"chromosome_y"])){
        var <- "dup_pair"
        break
      }else if((table[rownm,"gene_start"]<=duplication_table[i,"end_x"] && table[rownm,"gene_start"]>=duplication_table[i,"begin_x"] && table[rownm,"gene_chr"]==duplication_table[i,"chromosome_x"])){
        var <- "dup_same_chr"
        break
      }else{
        var1 <- "SNP_dup"
      }
    }else if((table[rownm,"snpstart"]<=duplication_table[i,"end_y"] && table[rownm,"snpstart"]>=duplication_table[i,"begin_y"] && table[rownm,"snpschrom"]==duplication_table[i,"chromosome_y"])){
      if((table[rownm,"gene_start"]<=duplication_table[i,"end_x"] && table[rownm,"gene_start"]>=duplication_table[i,"begin_x"] && table[rownm,"gene_chr"]==duplication_table[i,"chromosome_x"])){
        var <- "dup_pair"
        break
      }else if((table[rownm,"gene_start"]<=duplication_table[i,"end_y"] && table[rownm,"gene_start"]>=duplication_table[i,"begin_y"] && table[rownm,"gene_chr"]==duplication_table[i,"chromosome_y"])){
        var <- "dup_same_chr"
        break
      }else{
        var1 <- "SNP_dup"
      }
    }else if((table[rownm,"gene_start"]<=duplication_table[i,"end_y"] && table[rownm,"gene_start"]>=duplication_table[i,"begin_y"] && table[rownm,"gene_chr"]==duplication_table[i,"chromosome_y"]) | (table[rownm,"gene_start"]<=duplication_table[i,"end_x"] && table[rownm,"gene_start"]>=duplication_table[i,"begin_x"] && table[rownm,"gene_chr"]==duplication_table[i,"chromosome_x"])){
      var2 <- "gene_dup"
    }
  }
  if(var1 == "SNP_dup" && var2 == "gene_dup"){
    var <- "dup_no_pair"
  }else if(var1 == "SNP_dup" && var2 != "gene_dup"){
    var <- "SNP_dup"
  }else if(var1 != "SNP_dup" && var2 == "gene_dup"){
    var <- "gene_dup"
  }else if(var=="no"){
    var <- "nodup"
  }
  return(var)
}





### if i look at trans + cis SNPs
duplication_table<-read.table(file="https://salmobase.org/datafiles/TSV/synteny/2021-11/AtlanticSalmon/synteny.tsv",header=TRUE,sep="\t")

rna_table<-read.table(file="Synchro_eQTL_5Mwindow.csv",sep=",",header=TRUE)
rna_table$snpschrom<-sapply(strsplit(rna_table$SNP,"_"),FUN = `[[`, 1)
rna_table$snpschrom<-ifelse(substr(rna_table$snpschrom,4,4)==0,substr(rna_table$snpschrom,5,5),substr(rna_table$snpschrom,4,5))
rna_table$snpstart<-as.numeric(as.character(sapply(strsplit(rna_table$SNP,"_"),FUN = `[[`, 2)))
rna_table$gene_chr<-ifelse(substr(rna_table$gene_chr,1,1)=="s",ifelse(substr(rna_table$gene_chr,4,4)==0,substr(rna_table$gene_chr,5,5),substr(rna_table$gene_chr,4,5)),rna_table$gene_chr)
rna_table$snpstart<-as.numeric(as.character(rna_table$snpstart))
rna_table$gene_start<-as.numeric(as.character(rna_table$gene_start))




TestDuplicate7<- function(rownm,table){
  table_glo<-rna_table
  var<-"no"
  var1<-"no"
  var2<-"no"
  var3<-"nocis"
  for (i in 1:length(duplication_table$end_x)){
    if((table[rownm,"snpstart"]<=duplication_table[i,"end_x"] && table[rownm,"snpstart"]>=duplication_table[i,"begin_x"] && table[rownm,"snpschrom"]==duplication_table[i,"chromosome_x"])){
      if((table[rownm,"gene_start"]<=duplication_table[i,"end_y"] && table[rownm,"gene_start"]>=duplication_table[i,"begin_y"] && table[rownm,"gene_chr"]==duplication_table[i,"chromosome_y"])){
        var <- "dup_pair"
        Pos<-table[rownm,"snpstart"]
        chr<-table[rownm,"snpschrom"]
        subtab<-table_glo[table_glo$snpstart==Pos & table_glo$snpschrom==chr,]
        if("cis" %in% subtab$cistrans){
          var3 <- "alsocis"
        }
        break
      }else if((table[rownm,"gene_start"]<=duplication_table[i,"end_x"] && table[rownm,"gene_start"]>=duplication_table[i,"begin_x"] && table[rownm,"gene_chr"]==duplication_table[i,"chromosome_x"])){
        var <- "dup_same_chr"
        break
      }else{
        var1 <- "SNP_dup"
      }
    }else if((table[rownm,"snpstart"]<=duplication_table[i,"end_y"] && table[rownm,"snpstart"]>=duplication_table[i,"begin_y"] && table[rownm,"snpschrom"]==duplication_table[i,"chromosome_y"])){
      if((table[rownm,"gene_start"]<=duplication_table[i,"end_x"] && table[rownm,"gene_start"]>=duplication_table[i,"begin_x"] && table[rownm,"gene_chr"]==duplication_table[i,"chromosome_x"])){
        var <- "dup_pair"
        Pos<-table[rownm,"snpstart"]
        chr<-table[rownm,"snpschrom"]
        subtab<-table_glo[table_glo$snpstart==Pos & table_glo$snpschrom==chr,]
        if("cis" %in% subtab$cistrans){
          var3 <- "alsocis"
        }
        break
      }else if((table[rownm,"gene_start"]<=duplication_table[i,"end_y"] && table[rownm,"gene_start"]>=duplication_table[i,"begin_y"] && table[rownm,"gene_chr"]==duplication_table[i,"chromosome_y"])){
        var <- "dup_same_chr"
        break
      }else{
        var1 <- "SNP_dup"
      }
    }else if((table[rownm,"gene_start"]<=duplication_table[i,"end_y"] && table[rownm,"gene_start"]>=duplication_table[i,"begin_y"] && table[rownm,"gene_chr"]==duplication_table[i,"chromosome_y"]) | (table[rownm,"gene_start"]<=duplication_table[i,"end_x"] && table[rownm,"gene_start"]>=duplication_table[i,"begin_x"] && table[rownm,"gene_chr"]==duplication_table[i,"chromosome_x"])){
      var2 <- "gene_dup"
    }
  }
  if(var1 == "SNP_dup" && var2 == "gene_dup"){
    var <- "dup_no_pair"
  }else if(var1 == "SNP_dup" && var2 != "gene_dup"){
    var <- "SNP_dup"
  }else if(var1 != "SNP_dup" && var2 == "gene_dup"){
    var <- "gene_dup"
  }else if(var=="no"){
    var <- "nodup"
  }
  return(c(var,var3))
}

rna_table_trans<-rna_table[rna_table$cistrans=="trans",]
for(y in 1:length(rna_table_trans$snpschrom)){
  val<-TestDuplicate7(y,rna_table_trans)
  rna_table_trans[y,"dupstatus"]<-val[1]
  rna_table_trans[y,"cisgene"]<-val[2]
}
dups <- duplicated(subset(rna_table_trans, select = c("snpstart", "snpschrom")))
rna_table_trans<-rna_table_trans[!dups,]

dup_pairtable<-rna_table_trans[rna_table_trans$dupstatus=="dup_pair",]
library(dplyr)


dfcount <-dup_pairtable %>% group_by(cisgene) %>% 
  summarise(total_count=n(),.groups = 'drop') %>%
  as.data.frame()


library(ggplot2)
ggplot(dfcount,aes(x = cisgene,y = total_count,fill = cisgene)) +geom_bar(position = "dodge", stat = "identity")

### check frequency of homeologous SNPs

## for ll 0.05
dups <- duplicated(subset(llall, select = c("snpstart", "snpschrom")))
llalluniq<-llall[!dups & llall$dupstatus=="dup_pair",]
llalluniq_pos<-llalluniq[,3:4]
llalluniq_pos$snpschrom<-ifelse(llalluniq_pos$snpschrom<9,paste("ssa0",llalluniq_pos$snpschrom,sep=""),paste("ssa",llalluniq_pos$snpschrom,sep=""))

write.table(llalluniq_pos,file="/mnt/SCRATCH/cedi/phDSalmon/rna_expression/WGD_geneexp/data/llall_unique_duppair.tab",row.names = FALSE,quote=FALSE,col.names = FALSE,sep="\t")


freq_variants_llalluniq_pos<-read.table(file="/mnt/SCRATCH/cedi/phDSalmon/rna_expression/WGD_geneexp/snps_llall_unique_duppair.vaf.frq",header=FALSE,fill=TRUE)

freq_variants_llalluniq_pos<-freq_variants_llalluniq_pos[-1,]

freq_variants_llalluniq_pos<-freq_variants_llalluniq_pos[freq_variants_llalluniq_pos$V3!="",]
freq_variants_llalluniq_pos$freq<-as.numeric(as.character(sapply(strsplit(freq_variants_llalluniq_pos$V6,split=":"),FUN = `[[`, 2)))

hist(freq_variants_llalluniq_pos$freq)

## for ll 10-4

dups <- duplicated(subset(ll00001, select = c("snpstart", "snpschrom")))
ll00001uniq<-ll00001[!dups & ll00001$dupstatus=="dup_pair",]
ll00001uniq_pos<-ll00001uniq[,3:4]
ll00001uniq_pos$snpschrom<-ifelse(ll00001uniq_pos$snpschrom<9,paste("ssa0",ll00001uniq_pos$snpschrom,sep=""),paste("ssa",ll00001uniq_pos$snpschrom,sep=""))

write.table(ll00001uniq_pos,file="/mnt/SCRATCH/cedi/phDSalmon/rna_expression/WGD_geneexp/data/ll00001_unique_duppair.tab",row.names = FALSE,quote=FALSE,col.names = FALSE,sep="\t")


freq_variants_ll00001uniq_pos<-read.table(file="/mnt/SCRATCH/cedi/phDSalmon/rna_expression/WGD_geneexp/snps_ll00001_unique_duppair.vaf.frq",header=FALSE,fill=TRUE)

freq_variants_ll00001uniq_pos<-freq_variants_ll00001uniq_pos[-1,]

freq_variants_ll00001uniq_pos<-freq_variants_ll00001uniq_pos[freq_variants_ll00001uniq_pos$V3!="",]
freq_variants_ll00001uniq_pos$freq<-as.numeric(as.character(sapply(strsplit(freq_variants_ll00001uniq_pos$V6,split=":"),FUN = `[[`, 2)))

hist(freq_variants_ll00001uniq_pos$freq)



## for ll 10-10
dups <- duplicated(subset(ll0000000001, select = c("snpstart", "snpschrom")))
ll00000001uniq<-ll0000000001[!dups & ll0000000001$dupstatus=="dup_pair",]
ll00000001uniq_pos<-ll00000001uniq[,3:4]
ll00000001uniq_pos$snpschrom<-ifelse(ll00000001uniq_pos$snpschrom<9,paste("ssa0",ll00000001uniq_pos$snpschrom,sep=""),paste("ssa",ll00000001uniq_pos$snpschrom,sep=""))

write.table(ll00000001uniq_pos,file="/mnt/SCRATCH/cedi/phDSalmon/rna_expression/WGD_geneexp/data/ll00000001_unique_duppair.tab",row.names = FALSE,quote=FALSE,col.names = FALSE,sep="\t")


freq_variants_ll00000001uniq_pos<-read.table(file="/mnt/SCRATCH/cedi/phDSalmon/rna_expression/WGD_geneexp/snps_ll00000001_unique_duppair.vaf.frq",header=FALSE,fill=TRUE)

freq_variants_ll00000001uniq_pos<-freq_variants_ll00000001uniq_pos[-1,]

freq_variants_ll00000001uniq_pos<-freq_variants_ll00000001uniq_pos[freq_variants_ll00000001uniq_pos$V3!="",]
freq_variants_ll00000001uniq_pos$freq<-as.numeric(as.character(sapply(strsplit(freq_variants_ll00000001uniq_pos$V6,split=":"),FUN = `[[`, 2)))

hist(freq_variants_ll00000001uniq_pos$freq)

##
for (i in 1:length(tableMarie$CHROM)){
  chrom<-ifelse(substr(tableMarie[i,"CHROM"],4,4)==0,substr(tableMarie[i,"CHROM"],5,5),substr(tableMarie[i,"CHROM"],4,5))
  pos<-tableMarie[i,"POS"]
  tableMarie[i,"minpval"]<-min(rna_table_trans[rna_table_trans$snpschrom==chrom & rna_table_trans$snpstart==pos, "pval"])
}

write.table(tableMarie,file="/mnt/SCRATCH/cedi/phDSalmon/rna_expression/WGD_geneexp/TableMarie_pairsnpunique.txt",row.names = FALSE)

## checking the SNPS in all populations

FE_SV<-read.table(file="/mnt/SCRATCH/cedi/phDSalmon/rna_expression/WGD_geneexp/snps_ll00000001_unique_duppair_FE.vaf.frq",header=FALSE,fill=TRUE)
FE_SV<-FE_SV[-1,]

FE_SV<-FE_SV[FE_SV$V3!="",]
FE_SV$freq<-as.numeric(as.character(sapply(strsplit(FE_SV$V6,split=":"),FUN = `[[`, 2)))

FA_SV<-read.table(file="/mnt/SCRATCH/cedi/phDSalmon/rna_expression/WGD_geneexp/snps_ll00000001_unique_duppair_FA.vaf.frq",header=FALSE,fill=TRUE)
FA_SV<-FA_SV[-1,]

FA_SV<-FA_SV[FA_SV$V3!="",]
FA_SV$freq<-as.numeric(as.character(sapply(strsplit(FA_SV$V6,split=":"),FUN = `[[`, 2)))

WE_SV<-read.table(file="/mnt/SCRATCH/cedi/phDSalmon/rna_expression/WGD_geneexp/snps_ll00000001_unique_duppair_WE.vaf.frq",header=FALSE,fill=TRUE)
WE_SV<-WE_SV[-1,]

WE_SV<-WE_SV[WE_SV$V3!="",]
WE_SV$freq<-as.numeric(as.character(sapply(strsplit(WE_SV$V6,split=":"),FUN = `[[`, 2)))

WA_SV<-read.table(file="/mnt/SCRATCH/cedi/phDSalmon/rna_expression/WGD_geneexp/snps_ll00000001_unique_duppair_WA.vaf.frq",header=FALSE,fill=TRUE)
WA_SV<-WA_SV[-1,]

WA_SV<-WA_SV[WA_SV$V3!="",]
WA_SV$freq<-as.numeric(as.character(sapply(strsplit(WA_SV$V6,split=":"),FUN = `[[`, 2)))

tableMarie<-data.frame(CHROM=FE_SV$V1,POS=FE_SV$V2,freq_FE=FE_SV$freq,freq_FA=FA_SV$freq,freq_WE=WE_SV$freq,freq_WA=WA_SV$freq)

presence_FE<-ifelse(FE_SV$freq>0,1,0)
presence_FA<-ifelse(FA_SV$freq>0,1,0)
presence_WE<-ifelse(WE_SV$freq>0,1,0)
presence_WA<-ifelse(WA_SV$freq>0,1,0)


tableupset<-data.frame(CHROM=FE_SV$V1,POS=FE_SV$V2,Pres_FE=presence_FE,Pres_FA=presence_FA,Pres_WE=presence_WE,Pres_WA=presence_WA)





FE_SV<-read.table(file="/mnt/SCRATCH/cedi/phDSalmon/rna_expression/WGD_geneexp/snps_llall_unique_duppair_FE.vaf.frq",header=FALSE,fill=TRUE)
FE_SV<-FE_SV[-1,]

FE_SV<-FE_SV[FE_SV$V3!="",]
FE_SV$freq<-as.numeric(as.character(sapply(strsplit(FE_SV$V6,split=":"),FUN = `[[`, 2)))

FA_SV<-read.table(file="/mnt/SCRATCH/cedi/phDSalmon/rna_expression/WGD_geneexp/snps_llall_unique_duppair_FA.vaf.frq",header=FALSE,fill=TRUE)
FA_SV<-FA_SV[-1,]

FA_SV<-FA_SV[FA_SV$V3!="",]
FA_SV$freq<-as.numeric(as.character(sapply(strsplit(FA_SV$V6,split=":"),FUN = `[[`, 2)))

WE_SV<-read.table(file="/mnt/SCRATCH/cedi/phDSalmon/rna_expression/WGD_geneexp/snps_llall_unique_duppair_WE.vaf.frq",header=FALSE,fill=TRUE)
WE_SV<-WE_SV[-1,]

WE_SV<-WE_SV[WE_SV$V3!="",]
WE_SV$freq<-as.numeric(as.character(sapply(strsplit(WE_SV$V6,split=":"),FUN = `[[`, 2)))

WA_SV<-read.table(file="/mnt/SCRATCH/cedi/phDSalmon/rna_expression/WGD_geneexp/snps_llall_unique_duppair_WA.vaf.frq",header=FALSE,fill=TRUE)
WA_SV<-WA_SV[-1,]

WA_SV<-WA_SV[WA_SV$V3!="",]
WA_SV$freq<-as.numeric(as.character(sapply(strsplit(WA_SV$V6,split=":"),FUN = `[[`, 2)))

presence_FE<-ifelse(FE_SV$freq>0,1,0)
presence_FA<-ifelse(FA_SV$freq>0,1,0)
presence_WE<-ifelse(WE_SV$freq>0,1,0)
presence_WA<-ifelse(WA_SV$freq>0,1,0)


tableupset<-data.frame(CHROM=FE_SV$V1,POS=FE_SV$V2,Pres_FE=presence_FE,Pres_FA=presence_FA,Pres_WE=presence_WE,Pres_WA=presence_WA)


## SNPs not in dup pair
llall<-read.table(file = "Synchro_uniq_trans_withdupstatus_detailled_1_LL_tot.txt",header=TRUE)
dups <- duplicated(subset(llall, select = c("snpstart", "snpschrom")))
llalluniq<-llall[!dups & llall$dupstatus=="dup_no_pair",]
llalluniq_pos<-llalluniq[,3:4]
llalluniq_pos$snpschrom<-ifelse(llalluniq_pos$snpschrom<9,paste("ssa0",llalluniq_pos$snpschrom,sep=""),paste("ssa",llalluniq_pos$snpschrom,sep=""))

##random sample
llalluniq_pos_sample <- llalluniq_pos[sample(nrow(llalluniq_pos), 1000), ]

write.table(llalluniq_pos_sample,file="/mnt/SCRATCH/cedi/phDSalmon/rna_expression/WGD_geneexp/data/llall_unique_noduppair.tab",row.names = FALSE,quote=FALSE,col.names = FALSE,sep="\t")


##snp cis
cissnp<-rna_table[rna_table$cistrans=="cis",]
dups <- duplicated(subset(cissnp, select = c("snpstart", "snpschrom")))
cissnpuniq<-cissnp[!dups,]
cissnpuniq_pos<-cissnpuniq[,3:4]

cissnpuniq_pos_sample <- cissnpuniq_pos[sample(nrow(cissnpuniq_pos), 1000), ]

write.table(cissnpuniq_pos_sample,file="/mnt/SCRATCH/cedi/phDSalmon/rna_expression/WGD_geneexp/data/cissnp_unique.tab",row.names = FALSE,quote=FALSE,col.names = FALSE,sep="\t")




snpcount<-read.table(file="/mnt/SCRATCH/cedi/phDSalmon/rna_expression/WGD_geneexp/llall_unique_noduppairupsettable.txt",header=TRUE)

snpcount<-read.table(file="/mnt/SCRATCH/cedi/phDSalmon/rna_expression/WGD_geneexp/llall_unique_duppairupsettable.txt",header=TRUE)

snpcount<-read.table(file="/mnt/SCRATCH/cedi/phDSalmon/rna_expression/WGD_geneexp/ll00000001_unique_duppairupsettable.txt",header=TRUE)

snpcount<-read.table(file="/mnt/SCRATCH/cedi/phDSalmon/rna_expression/WGD_geneexp/subSNPupsettable.txt",header=TRUE)
sum(snpcount$total_count)

snpcount<-read.table(file="/mnt/SCRATCH/cedi/phDSalmon/rna_expression/WGD_geneexp/cissnp_uniqueupsettable.txt",header=TRUE)
sum(snpcount$total_count)

##cumulative plot with all data
rna_table_trans$pvalgroup<-ifelse(rna_table_trans$pval<=0.0000000001,10^-10,ifelse(rna_table_trans$pval>0.0000000001 & rna_table_trans$pval<=0.0001,10^-4,ifelse(rna_table_trans$pval>0.0001 & rna_table_trans$pval<=0.01,10^-2,5*10^-2)))

library(dplyr)


dfcount <-rna_table_trans %>% group_by(pvalgroup,dupstatus) %>% 
  summarise(total_count=n(),.groups = 'drop') %>%
  as.data.frame()


library(ggplot2)
ggplot(dfcount,aes(x = as.factor(pvalgroup),y = total_count,fill = dupstatus)) +geom_bar(position = "dodge", stat = "identity")

rna_table_trans$dupstatus[rna_table_trans$dupstatus=="dup_same_chr"]<-"dup_same_syntenic_region"

ggplot(rna_table_trans, aes(log10(pval),color=dupstatus)) + stat_ecdf(geom = "point")


rna_table_trans_sub<-rna_table_trans[rna_table_trans$dupstatus%in%c("dup_same_syntenic_region","dup_pair","dup_no_pair"),]
ggplot(rna_table_trans_sub, aes(log10(pval),color=dupstatus)) + stat_ecdf(geom = "point")+scale_color_manual(values=c("firebrick","gold2","forestgreen"))


ks.test(log10(rna_table_trans_sub[rna_table_trans_sub$dupstatus=="dup_pair","pval"]),log10(rna_table_trans_sub[rna_table_trans_sub$dupstatus=="dup_no_pair","pval"]))

ks.test(log10(rna_table_trans_sub[rna_table_trans_sub$dupstatus=="dup_pair","pval"]),log10(rna_table_trans_sub[rna_table_trans_sub$dupstatus=="dup_same_syntenic_region","pval"]))

ks.test(log10(rna_table_trans_sub[rna_table_trans_sub$dupstatus=="dup_same_syntenic_region","pval"]),log10(rna_table_trans_sub[rna_table_trans_sub$dupstatus=="dup_no_pair","pval"]))



## performing percentage values...
library(dplyr)

llall<-read.table(file = "Synchro_uniq_trans_withdupstatus_detailled_1_LL_tot.txt",header=TRUE)
llall<-read.table(file = "Synchro_uniq_trans_withdupstatus_detailled_1e-10_12_12_tot.txt",header=TRUE)

dfcount <- llall %>% group_by(dupstatus,typedata) %>% 
  summarise(total_count=n(),.groups = 'drop') %>%
  as.data.frame()

#write.table(dfcount,file="dfcount_LL_alldata.txt",row.names=FALSE)

ratio_real_data_dup_pair_vs_all_others<-dfcount[dfcount$typedata=="realdata" & dfcount$dupstatus=="dup_pair","total_count"]/sum(dfcount[dfcount$typedata=="realdata" & dfcount$dupstatus!="dup_pair","total_count"])
ratio_shuffled_dup_pair_vs_all_others<-dfcount[dfcount$typedata=="realdata" & dfcount$dupstatus=="dup_pair","total_count"]/sum(dfcount[dfcount$typedata=="shuffled" & dfcount$dupstatus!="dup_pair","total_count"])

ratio_real_data_dup_pair_vs_all_others/ratio_shuffled_dup_pair_vs_all_others ##1.00994 ; #1.034949

ratio_real_data_dup_pair_vs_yellow_red<-dfcount[dfcount$typedata=="realdata" & dfcount$dupstatus=="dup_pair","total_count"]/sum(dfcount[dfcount$typedata=="realdata" & dfcount$dupstatus%in%c("dup_no_pair","dup_same_chr"),"total_count"])
ratio_shuffled_dup_pair_vs_yellow_red<-dfcount[dfcount$typedata=="realdata" & dfcount$dupstatus=="dup_pair","total_count"]/sum(dfcount[dfcount$typedata=="shuffled" & dfcount$dupstatus%in%c("dup_no_pair","dup_same_chr"),"total_count"])


ratio_real_data_dup_pair_vs_yellow_red/ratio_shuffled_dup_pair_vs_yellow_red ##1.011797 ; #1.040155


ratio_real_data_green_yellow_vs_red<-sum(dfcount[dfcount$typedata=="realdata" & dfcount$dupstatus%in%c("dup_pair","dup_same_chr"),"total_count"])/sum(dfcount[dfcount$typedata=="realdata" & dfcount$dupstatus=="dup_no_pair","total_count"])
ratio_shuffled_green_yellow_vs_red<-sum(dfcount[dfcount$typedata=="realdata" & dfcount$dupstatus%in%c("dup_pair","dup_same_chr"),"total_count"])/sum(dfcount[dfcount$typedata=="shuffled" & dfcount$dupstatus=="dup_no_pair","total_count"])

ratio_real_data_green_yellow_vs_red/ratio_shuffled_green_yellow_vs_red ##1.190715 ; #1.886525


ratio_real_data_dup_pair_vs_all_except_green<-dfcount[dfcount$typedata=="realdata" & dfcount$dupstatus=="dup_pair","total_count"]/sum(dfcount[dfcount$typedata=="realdata" & dfcount$dupstatus%in%c("dup_no_pair","dup_same_chr","gene_dup","SNP_dup"),"total_count"])
ratio_shuffled_dup_pair_vs_all_except_green<-dfcount[dfcount$typedata=="realdata" & dfcount$dupstatus=="dup_pair","total_count"]/sum(dfcount[dfcount$typedata=="shuffled" & dfcount$dupstatus%in%c("dup_no_pair","dup_same_chr","gene_dup","SNP_dup"),"total_count"])

ratio_real_data_dup_pair_vs_all_except_green/ratio_shuffled_dup_pair_vs_all_except_green ##1.01   #1.03



## check some SNPs region


library("gplots")
library(ggplot2); library(reshape2)
library(tidyverse)
#vcf1<-select(vcf,  -c(1,3:9))

ssa08_region<-read.table(file="/mnt/SCRATCH/cedi/phDSalmon/rna_expression/WGD_geneexp/ssa08snp.vcf.gz_258.stat",header=FALSE)
ssa08sample<-read.table(file="/mnt/SCRATCH/cedi/phDSalmon/rna_expression/WGD_geneexp/samplessa08.txt",header=FALSE)
sampname<-ssa08sample$V1
colnames(ssa08_region)<-c("POS",sampname)

ssa16_region<-read.table(file="/mnt/SCRATCH/cedi/phDSalmon/rna_expression/WGD_geneexp/ssa16snp.vcf.gz_258.stat",header=FALSE)
ssa16sample<-read.table(file="/mnt/SCRATCH/cedi/phDSalmon/rna_expression/WGD_geneexp/samplessa16.txt",header=FALSE)
sampname<-ssa16sample$V1
colnames(ssa16_region)<-c("POS",sampname)

ssa05_region<-read.table(file="/mnt/SCRATCH/cedi/phDSalmon/rna_expression/WGD_geneexp/ssa05snp.vcf.gz_258.stat",header=FALSE)
ssa05sample<-read.table(file="/mnt/SCRATCH/cedi/phDSalmon/rna_expression/WGD_geneexp/samplessa05.txt",header=FALSE)
sampname<-ssa05sample$V1
colnames(ssa05_region)<-c("POS",sampname)


ssa05_region_phased<-read.table(file="/mnt/SCRATCH/cedi/phDSalmon/rna_expression/WGD_geneexp/ssa05snp_phased.vcf.gz_258.stat",header=FALSE)
ssa05sample<-read.table(file="/mnt/SCRATCH/cedi/phDSalmon/rna_expression/WGD_geneexp/samplessa05_phased.txt",header=FALSE)
sampname<-ssa05sample$V1
colnames(ssa05_region_phased)<-c("POS",sampname)


ssa16_region_phased<-read.table(file="/mnt/SCRATCH/cedi/phDSalmon/rna_expression/WGD_geneexp/ssa16snp_phased.vcf.gz_258.stat",header=FALSE)
ssa16sample<-read.table(file="/mnt/SCRATCH/cedi/phDSalmon/rna_expression/WGD_geneexp/samplessa16_phased.txt",header=FALSE)
sampname<-ssa16sample$V1
colnames(ssa16_region_phased)<-c("POS",sampname)

#vcf1<-ind13gt[,-c(1,3)]
vcf1<-ssa08_region
vcf1<-ssa16_region
vcf1<-ssa05_region


vcf1[vcf1=="0/0"]<- "0"
vcf1[vcf1=="0/1"]<- "1"
vcf1[vcf1=="1/0"]<- "1"
vcf1[vcf1=="1/1"]<- "2"
vcf1[vcf1=="./."]<- "0"
vcf1[vcf1=="1|0"]<- "1"
vcf1[vcf1=="0|1"]<- "1"
vcf1[vcf1=="1|1"]<- "2"


vcf1<-ssa05_region_phased
vcf1<-ssa16_region_phased

vcf1[vcf1=="0|0"]<- "0"
vcf1[vcf1=="1|0"]<- "1"
vcf1[vcf1=="0|1"]<- "2"
vcf1[vcf1=="1|1"]<- "3"

#library(dplyr)
#colnames(vcf3) <- factor(row.names(vcf))
vcf3 <- data.frame(apply(vcf1, 2, function(x) as.numeric(as.character(x))))
#vcf33<-vcf3[c(1:1000),c(1:271)]

vcf2 <- t(vcf3)
vcf4 <- data.matrix(vcf2)
colnames(vcf4) <- factor(row.names(vcf1))
#vcf4[1:3,1:10]
vcf4 %>% janitor::row_to_names(1) -> vcf4


## subset samples
#select.sample <- subset(color0,  pop!="Canada", select=sample)
## color
dev.off()
color0  <- read.csv("/mnt/SCRATCH/cedi/phDSalmon/evolutionary_analysis/pauline_data/sampleSRinfo.csv", header = T, sep = ",")
pop1 = as.character(color0$color)


my_palette <- colorRampPalette(c("#e7e7e7","#ababab","#000000")) (n=3)
heatmap.2(vcf4, trace="none", na.color = "black",scale="none", margins=c(12,8),
          col = my_palette,density.info="none",dendrogram = "none",RowSideColors=pop1, Rowv=FALSE, Colv=FALSE)

heatmap.2(vcf4, trace="none", na.color = "black",scale="none", margins=c(12,8),
          col = my_palette,density.info="none",dendrogram = "none",RowSideColors=pop1,Colv=FALSE)


my_palette <- colorRampPalette(c("#e7e7e7","#ababab","#616262","#000000")) (n=4)
heatmap.2(vcf4, trace="none", na.color = "black",scale="none", margins=c(12,8),
          col = my_palette,density.info="none",dendrogram = "none",RowSideColors=pop1[1:369], Rowv=FALSE, Colv=FALSE)

heatmap.2(vcf4, trace="none", na.color = "black",scale="none", margins=c(12,8),
          col = my_palette,density.info="none",dendrogram = "none",RowSideColors=pop1[1:369],Colv=FALSE)

abline(v=0.515,col="red")

## gene analysis
rna_table<-read.table(file="Synchro_eQTL_5Mwindow.csv",sep=",",header=TRUE)

rna_table_trans<-read.table(file="Synchro_uniq_trans_withdupstatus_detailled.txt",header=TRUE)

rna_table_trans<-rna_table_trans[rna_table_trans$dupstatus=="dup_pair" & rna_table_trans$condition=="LL",]
genetab<-read.table(file="ss4r_dups_and_singletons_ENSrapid_convPipeline.tsv",header=TRUE)



for (i in 1:length(rna_table_trans$condition)){
  genetocheck<-rna_table_trans[i,"gene"]
  if(genetocheck%in%genetab$gene1){
    tab<-genetab[genetab$gene1==genetocheck,]
    gene<-tab[!is.na(tab$type),"gene2"]
  }else if(genetocheck%in%genetab$gene2){
    tab<-genetab[genetab$gene2==genetocheck,]
    gene<-tab[!is.na(tab$type),"gene1"]
  }else{
    gene<-"NA"
  }
  rna_table_trans[i,"dupgene"]<-gene
}

## check if duplicate is also cis regulated

cistab<-rna_table[rna_table$cistrans=="cis",]
rna_tab_trans_nona<-rna_table_trans
for (i in 1:length(rna_tab_trans_nona$condition)){
  genetotest<-rna_tab_trans_nona[i,"dupgene"]
  if(is.na(genetotest)){
    cisregu<-"none"
  }else{
    
    if(genetotest%in%cistab$gene){
      tabchr<-cistab[cistab$gene==genetotest,]
      chr<-unique(tabchr[!is.na(tabchr$gene),"gene_chr"])
      treatedchr<-ifelse(substr(chr,4,4)==0,substr(chr,5,5),substr(chr,4,5))
      snpvect<-unique(tabchr[!is.na(tabchr$gene),"SNP"])
      snptreatedvect<-as.numeric(as.character(sapply(strsplit(snpvect,"_"),FUN = `[[`, 2)))
      if(treatedchr==rna_tab_trans_nona[i,"snpschrom"]){
        if(rna_tab_trans_nona[i,"snpstart"]%in%snptreatedvect){
          cisregu<-"duplicateallcond"
        }else{
          cisregu<-"notsameSNP"
        }
      }else{
        cisregu<-"other"
      }
    }else{
      cisregu<-"nogene"
    }
  }
  rna_tab_trans_nona[i,"cisregu"]<-cisregu
}

cistab<-rna_table[rna_table$cistrans=="cis" & rna_table$condition=="LL",]
for (i in 1:length(rna_tab_trans_nona$condition)){
  genetotest<-rna_tab_trans_nona[i,"dupgene"]
  if(is.na(genetotest)){
    cisregu<-"none"
  }else{
    if(genetotest%in%cistab$gene){
      tabchr<-cistab[cistab$gene==genetotest,]
      chr<-unique(tabchr[!is.na(tabchr$gene),"gene_chr"])
      treatedchr<-ifelse(substr(chr,4,4)==0,substr(chr,5,5),substr(chr,4,5))
      snpvect<-unique(tabchr[!is.na(tabchr$gene),"SNP"])
      snptreatedvect<-as.numeric(as.character(sapply(strsplit(snpvect,"_"),FUN = `[[`, 2)))
      if(treatedchr==rna_tab_trans_nona[i,"snpschrom"]){
        if(rna_tab_trans_nona[i,"snpstart"]%in%snptreatedvect){
          cisregu<-"duplicatesamecond"
        }else{
          cisregu<-"notsameSNP"
        }
      }else{
        cisregu<-"other"
      }
    }else{
      cisregu<-"nogene"
    }
  }
  rna_tab_trans_nona[i,"cisregu2"]<-cisregu
}


library(dplyr)


dfcount <- rna_tab_trans_nona %>% group_by(cisregu,cisregu2) %>% 
  summarise(total_count=n(),.groups = 'drop') %>%
  as.data.frame()




## continuation of investigation with complete table

rna_full_table<-read.table(file="Synchro_uniq_fulltable_withinfo_03.txt",header=TRUE)

rna_full_table<-rna_full_table[rna_full_table$condition=="LL",]

#rna_full_table[rna_full_table$dupstatus=="dup_diff_chr","dupstatus"]<-"inter_chr_non-ohnolog"
#rna_full_table[rna_full_table$dupstatus=="dup_same_syntenic_region","dupstatus"]<-"intra_chr_syntenic"
#rna_full_table[rna_full_table$dupstatus=="dup_same_chr","dupstatus"]<-"intra_chr_non-syntenic"
#rna_full_table[rna_full_table$dupstatus=="dup_pair","dupstatus"]<-"ohnolog"

#write.table(rna_full_table,file="Synchro_uniq_fulltable_withinfo_03.txt",row.names=FALSE)

library(dplyr)
library(ggplot2)

dfcount <- rna_full_table %>% group_by(dupstatus,cispres) %>% 
  summarise(total_count=n(),.groups = 'drop') %>%
  as.data.frame()

#cispres: 19000 no, 14000 yes
#library(forcats)


#making the right input format
#reorder variable
dfcount$dupstatus<-as.factor(dfcount$dupstatus)
rna_full_table$dupstatus<-as.factor(rna_full_table$dupstatus)

dfcount$dupstatus<-factor(dfcount$dupstatus,levels=c("inter_chr_non-ohnolog","ohnolog","intra_chr_non-syntenic","intra_chr_syntenic","gene_dup_only","SNP_dup_only","none"))
rna_full_table$dupstatus<-factor(rna_full_table$dupstatus,levels=c("inter_chr_non-ohnolog","ohnolog","intra_chr_non-syntenic","intra_chr_syntenic","gene_dup_only","SNP_dup_only","none"))

ggplot(dfcount[dfcount$dupstatus%in%c("inter_chr_non-ohnolog","intra_chr_syntenic","intra_chr_non-syntenic","ohnolog") & !is.na(dfcount$cispres),],aes(x=dupstatus,y=total_count,fill=cispres))+geom_bar(position = "fill", stat = "identity")+theme_bw(base_size = 20)+ggtitle("cis regulating gene presence for the trans-snp PAIRS")


ggplot(rna_full_table[rna_full_table$dupstatus%in%c("inter_chr_non-ohnolog","intra_chr_syntenic","intra_chr_non-syntenic","ohnolog"),], aes(log10(pval),color=dupstatus)) + stat_ecdf(geom = "point")+scale_color_manual(values=c("firebrick","gold2","darkorchid1","forestgreen"))+theme_bw(base_size = 20)+
  ggtitle("cumulative distribution of p-value for different duplication status")

ggplot(rna_full_table[rna_full_table$dupstatus%in%c("inter_chr_non-ohnolog","intra_chr_syntenic","intra_chr_non-syntenic","ohnolog"),], aes(log10(pval),color=dupstatus)) + stat_ecdf(geom = "point")+scale_color_manual(values=c("firebrick","gold2","darkorchid1","forestgreen"))+theme_bw(base_size = 20)+
  ggtitle("cumulative distribution of absolute effect on gene expression for different duplication status")

rna_full_table$counterpart<-ifelse(is.na(rna_full_table$pairID),"no","yes")

dfcount2 <- rna_full_table[rna_full_table$dupstatus=="dup_pair",] %>% group_by(counterpart,cisduplicate) %>% 
  summarise(total_count=n(),.groups = 'drop') %>%
  as.data.frame()

dupstatus_list <- c("inter_chr_non-ohnolog","intra_chr_syntenic","intra_chr_non-syntenic","ohnolog")

rna_full_table$dupstatus<-as.character(rna_full_table$dupstatus)
# Create a new dataframe with counts of each unique combination of snpschrom and snpstart
count_df <- rna_full_table %>% 
  group_by(snpschrom, snpstart) %>% 
  summarise(count = n(), dupstatus = paste0(unique(unlist(strsplit(dupstatus, ",")))), .groups = "drop") %>%
  filter(dupstatus %in% dupstatus_list)

##count only cis

count_df <- rna_full_table %>%
  group_by(snpschrom, snpstart) %>%
  summarise(cistrans_check = any(cistrans == "trans"),
            count_cis = sum(cistrans == "cis"),
            dupstatus = paste0(unique(unlist(strsplit(dupstatus, ",")))), 
            count = n(), 
            .groups = "drop") %>%
  filter(cistrans_check, dupstatus %in% dupstatus_list) %>%
  select(-cistrans_check)

# Create a ggplot histogram of the counts, colored by dupstatus


ggplot(count_df, aes(x = count, fill = dupstatus)) + geom_histogram(position=position_dodge(),binwidth=1) + xlim(-1,30)

#test<-count_df[count_df$count>300,]

percent_df <- count_df %>% 
  filter(dupstatus %in% dupstatus_list) %>% 
  group_by(dupstatus, count) %>% 
  summarise(n = n()) %>% 
  mutate(percent = n / sum(n) * 100)

percent_df$dupstatus<-as.factor(percent_df$dupstatus)
percent_df$dupstatus<-factor(percent_df$dupstatus,levels=c("inter_chr_non-ohnolog","ohnolog","intra_chr_non-syntenic","intra_chr_syntenic"))


# Plot percentages
ggplot(percent_df, aes(x = count, y = percent, fill = dupstatus)) + 
  geom_bar(stat = "identity", position = position_dodge(), alpha = 0.5, color = "black") +
  labs(x = "Count", y = "Percentage", title = "Distribution of gene regulated by SNPs count by duplication status") +
  theme_classic()+xlim(-1,20)+scale_fill_manual(values=c("firebrick","gold2","darkorchid1","forestgreen"))+theme_bw(base_size = 20)

### get gene with count > x to do GO analysis

snps<-count_df[count_df$count>5 & count_df$count<50 & count_df$dupstatus=="ohnolog",]
genes<-rna_full_table[rna_full_table$snpschrom%in%snps$snpschrom & rna_full_table$snpstart%in%snps$snpstart,"gene"]
write.table(data.frame(genes),file="genes_LL_trans.txt",row.names = FALSE,col.names =FALSE,quote=FALSE)
##re do the cis pres plot in another way

rna_full_table <- rna_full_table %>%
  mutate(dist = abs(snpstart - gene_start))

new_table <- rna_full_table %>%
  filter(!is.na(cispres)) %>%
  group_by(snpschrom, snpstart, dupstatus) %>%
  summarize(cispres = ifelse(any(cispres == "yes"), "yes", "no")) %>%
  ungroup() %>% filter(dupstatus %in% dupstatus_list)

#dupstatus = toString(unique(dupstatus))

dfcount <- new_table %>% group_by(dupstatus,cispres) %>% 
  summarise(total_count=n(),.groups = 'drop') %>%
  as.data.frame()

dfcount$dupstatus<-as.factor(dfcount$dupstatus)
dfcount$dupstatus<-factor(dfcount$dupstatus,levels=c("inter_chr_non-ohnolog","ohnolog","intra_chr_non-syntenic","intra_chr_syntenic"))

ggplot(dfcount[dfcount$dupstatus%in%dupstatus_list & !is.na(dfcount$cispres),],aes(x=dupstatus,y=total_count,fill=cispres))+geom_bar(position = "fill", stat = "identity")+theme_bw(base_size = 20)+ggtitle("cis regulating gene presence for unique SNPs")



## makes the differences between close cis or not close cis

rna_full_table<-read.table(file="Synchro_uniq_fulltable_withinfo_04.txt",header=TRUE)

dupstatus_list <- c("inter_chr_non-ohnolog","intra_chr_syntenic","intra_chr_non-syntenic","ohnolog")

library(dplyr)
library(tidyr)
library(ggplot2)

rna_full_table <- rna_full_table %>%
  mutate(dist = abs(snpstart - gene_start))

##add row with min and max dist for SNPs found

new_table <- rna_full_table %>%
  group_by(snpschrom, snpstart) %>%
  mutate(has_cis = any(cistrans == "cis")) %>%
  mutate(max_dist_cis = ifelse(has_cis, max(ifelse(cistrans == "cis", dist, NA), na.rm = TRUE), NA),
         min_dist_cis = ifelse(has_cis, min(ifelse(cistrans == "cis", dist, NA), na.rm = TRUE), NA)) %>%
  ungroup()

#test<-new_table[new_table$snpstart==34078385 & new_table$snpschrom==14,]

new_table2 <- new_table %>%
  filter(!is.na(cispres)) %>%
  group_by(snpschrom, snpstart) %>%
  mutate(cis_status = case_when(
    cispres == "no" ~ "no_cis",
    max_dist_cis > 50000 & min_dist_cis < 50000 ~ "far_close_cis",
    max_dist_cis > 50000 & min_dist_cis >= 50000 ~ "far_cis",
    max_dist_cis <= 50000 & min_dist_cis < 50000 ~ "close_cis",
    TRUE ~ NA_character_
  ))



#test2<-new_table2[new_table2$snpstart==34078385 & new_table2$snpschrom==14,]

count_df <- new_table2 %>% 
  group_by(snpschrom, snpstart, cis_status) %>% 
  summarise(dupstatus = paste0(unique(unlist(strsplit(dupstatus, ",")))), .groups = "drop") %>%
  filter(dupstatus %in% dupstatus_list) %>% group_by(dupstatus, cis_status) %>% summarise(count = n())


count_df$cis_status<-as.factor(count_df$cis_status)

count_df$cis_status<-factor(count_df$cis_status,levels=c("no_cis","close_cis","far_close_cis","far_cis"))



ggplot(count_df,aes(x=dupstatus,y=count,fill=cis_status))+geom_bar(position = "fill", stat = "identity")+theme_bw(base_size = 20)+ggtitle("cis regulating gene presence for unique SNPs")
ggplot(count_df[count_df$cis_status!="no_cis",],aes(x=dupstatus,y=count,fill=cis_status))+geom_bar(position = "fill", stat = "identity")+theme_bw(base_size = 20)+ggtitle("cis regulating gene presence for unique SNPs")


##code that gives different results
# compute distance between snpstart and gene_start
df <- rna_full_table %>%
  mutate(dist = abs(snpstart - gene_start))
# group by unique combination of snpstart and snpschrom
df_grouped <- df %>%
  group_by(snpschrom, snpstart) %>%
  summarise(has_trans = "trans" %in% cistrans,
            has_cis = "cis" %in% cistrans,
            has_close_cis = any(dist < 50000 & cistrans == "cis"),
            has_far_cis = any(dist >= 50000 & cistrans == "cis"),
            dupstatus = toString(unique(dupstatus))) %>%
  separate_rows(dupstatus, sep = ", ") %>% filter(dupstatus %in% dupstatus_list) %>% filter(has_trans != FALSE)


dfcount <- df_grouped %>% group_by(dupstatus,has_cis,has_far_cis,has_close_cis) %>% 
  summarise(total_count=n(),.groups = 'drop') %>%
  as.data.frame()


df <- dfcount %>%
  mutate(categorie = case_when(
    !has_cis ~ "no_cis",
    has_cis & has_close_cis & !has_far_cis ~ "close_cis",
    has_cis & !has_close_cis & has_far_cis ~ "far_cis",
    has_cis & has_close_cis & has_far_cis ~ "close_and_far_cis",
    TRUE ~ NA_character_
  ))

ggplot(df,aes(x=dupstatus,y=total_count,fill=categorie))+geom_bar(position = "fill", stat = "identity")+theme_bw(base_size = 20)+ggtitle("cis regulating gene presence for the unique trans SNPs")




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

#test<-which(IDfile$CHROM==rna_full_table[1250,"snpschrom"] & IDfile$POS<rna_full_table[1250,"snpstart"] & IDfile$END>rna_full_table[1250,"snpstart"])

for (i in 1:length(rna_full_table$condition)){
  matchingID<-which(IDfile$CHROM==rna_full_table[i,"snpschrom"] & IDfile$POS<rna_full_table[i,"snpstart"] & IDfile$END>rna_full_table[i,"snpstart"])
  if(length(matchingID)==0){
    rna_full_table[i,"snp_identity"]<-NA
  }else{
    rna_full_table[i,"snp_identity"]<-IDfile[matchingID[1],"identity"]
  }
  
  matchingID2<-which(IDfile$CHROM==rna_full_table[i,"gene_chr"] & IDfile$POS<rna_full_table[i,"gene_start"] & IDfile$END>rna_full_table[i,"gene_start"])
  if(length(matchingID2)==0){
    rna_full_table[i,"gene_identity"]<-NA
  }else{
    rna_full_table[i,"gene_identity"]<-IDfile[matchingID2[1],"identity"]
  }
  
}

rna_full_table$meanidentity<-(rna_full_table$snp_identity+rna_full_table$gene_identity)/2
rna_full_table$distance<-abs(rna_full_table$snp_identity-rna_full_table$gene_identity)
dupstatus_list <- c("inter_chr_non-ohnolog","intra_chr_syntenic","intra_chr_non-syntenic","ohnolog")

ggplot(rna_full_table[rna_full_table$dupstatus%in%dupstatus_list,],aes(x=dupstatus,y=snp_identity,fill=dupstatus))+geom_boxplot()+theme_bw(base_size = 20)

ggplot(rna_full_table[rna_full_table$dupstatus%in%dupstatus_list,],aes(x=dupstatus,y=gene_identity,fill=dupstatus))+geom_boxplot()+theme_bw(base_size = 20)

ggplot(rna_full_table[rna_full_table$dupstatus%in%dupstatus_list,],aes(x=dupstatus,y=meanidentity,fill=dupstatus))+geom_boxplot()+theme_bw(base_size = 20)

ggplot(rna_full_table[rna_full_table$dupstatus%in%dupstatus_list,],aes(x=dupstatus,y=distance,fill=dupstatus))+geom_boxplot()+theme_bw(base_size = 20)


## shuffled vs non shuffled
rna_full_table<-read.table(file="Synchro_uniq_fulltable_withinfo_04.txt",header=TRUE)
rna_full_table_shuffled<-read.table(file="Synchro_uniq_fulltable_withinfo_shuffled_03.txt",header=TRUE)

rna_full_table<-rna_full_table[rna_full_table$condition=="LL",]

#rna_ohnologs<-rna_full_table[rna_full_table$dupstatus=="ohnolog","gene"]
#write.table(rna_ohnologs,file="ohnologs_genes.txt",row.names = FALSE,col.names = FALSE,quote=FALSE)
#rna_all<-rna_full_table[rna_full_table$dupstatus%in%dupstatus_list,"gene"]
#write.table(rna_all,file="all_genes.txt",row.names = FALSE,col.names = FALSE,quote=FALSE)

rna_full_table$data<-"real"
rna_full_table_shuffled$data<-"shuffled"

fulltab<-rbind.data.frame(rna_full_table,rna_full_table_shuffled,stringsAsFactors = FALSE)

library(dplyr)
library(ggplot2)
library(tidyr)

dupstatus_list <- c("inter_chr_non-ohnolog","intra_chr_syntenic","intra_chr_non-syntenic","ohnolog")

dfcount <- fulltab%>% filter(!is.na(cispres)) %>% group_by(dupstatus,data) %>% 
  summarise(total_count=n(),.groups = 'drop') %>%
  as.data.frame() %>% filter(dupstatus %in% dupstatus_list) %>% group_by(data) %>% 
  mutate(prop = total_count / sum(total_count))

dfcount$dupstatus<-as.factor(dfcount$dupstatus)
dfcount$dupstatus<-factor(dfcount$dupstatus,levels=c("inter_chr_non-ohnolog","ohnolog","intra_chr_non-syntenic","intra_chr_syntenic"))


ggplot(dfcount,aes(x=data,y=prop,fill=dupstatus))+geom_bar(position = "fill", stat = "identity")+scale_fill_manual(values=c("firebrick","gold2","darkorchid1","forestgreen"))+theme_bw(base_size = 20)

ggplot(dfcount,aes(x="",y=prop,fill=dupstatus))+geom_bar(stat="identity",width=1)+
  coord_polar("y",start=0)+geom_text(aes(label=paste0(prop,"%")),position=position_stack(vjust=0.5))+labs(x=NULL,y=NULL,fill=NULL)




dfcount2<-dfcount
dfcount2$dupstatus<-as.character(dfcount2$dupstatus)
dfcount2[dfcount2$dupstatus=="inter_chr_non-ohnolog","dupstatus"]<-"A"
dfcount2[dfcount2$dupstatus=="ohnolog","dupstatus"]<-"B"
dfcount2[dfcount2$dupstatus=="intra_chr_non-syntenic","dupstatus"]<-"C"
dfcount2[dfcount2$dupstatus=="intra_chr_syntenic","dupstatus"]<-"D"

dfcount2[dfcount2$dupstatus=="inter_chr_non-ohnolog","dupstatus"]<-"Inter-chromosomal non ohnolog"
dfcount2[dfcount2$dupstatus=="ohnolog","dupstatus"]<-"Inter-chromosomal ohnolog"
dfcount2[dfcount2$dupstatus=="intra_chr_non-syntenic","dupstatus"]<-"Intra-chromosomal non syntenic"
dfcount2[dfcount2$dupstatus=="intra_chr_syntenic","dupstatus"]<-"Intra-chromosomal syntenic"
ggplot(dfcount2,aes(x=data,y=prop,fill=dupstatus))+geom_bar(position = "fill", stat = "identity")+scale_fill_manual(values=c("#F3B5B3","#F5E69A","#B5AFDE","#B6D8B6"))+theme_bw(base_size = 30)+ylab("Proportion") +xlab("Data type")+guides(fill=guide_legend("Interaction type"))

ggplot(dfcount2,aes(x="",y=prop,fill=dupstatus))+geom_bar(stat="identity",width=1)+
  coord_polar("y",start=0)+theme_void()+scale_fill_manual(values=c("#F3B5B3","#F5E69A","#B5AFDE","#B6D8B6"))+
  theme(legend.text = element_text(size=15))



df_wide <- pivot_wider(dfcount, id_cols = dupstatus, names_from = data, values_from = total_count)
df_wide$enrichment<-df_wide$real/df_wide$shuffled

##use bootstrap values
intra_syntenic<-read.table(file="bootstrap_intra_chr_syntenic_01.txt")
intra_non_syntenic<-read.table(file="bootstrap_intra_chr_non-syntenic_01.txt")
ohnolgs<-read.table(file="bootstrap_ohnolog_01.txt")
non_ohnlogs<-read.table(file="bootstrap_inter_chr_non-ohnolog_01.txt")

df_wide$dupstatus<-as.character(df_wide$dupstatus)
df_wide[df_wide$dupstatus=="inter_chr_non-ohnolog","dupstatus"]<-"Inter-chromosomal non ohnolog"
df_wide[df_wide$dupstatus=="ohnolog","dupstatus"]<-"Inter-chromosomal ohnolog"
df_wide[df_wide$dupstatus=="intra_chr_non-syntenic","dupstatus"]<-"Intra-chromosomal non syntenic"
df_wide[df_wide$dupstatus=="intra_chr_syntenic","dupstatus"]<-"Intra-chromosomal syntenic"

df_wide$enrichment<-c(mean(non_ohnlogs$V1),mean(intra_non_syntenic$V1),mean(intra_syntenic$V1),mean(ohnolgs$V1))

df_wide$dupstatus<-factor(df_wide$dupstatus,levels=c("Inter-chromosomal non ohnolog","Inter-chromosomal ohnolog","Intra-chromosomal non syntenic","Intra-chromosomal syntenic"))

ggplot(df_wide,aes(x=dupstatus,y=enrichment,fill=dupstatus))+geom_bar(stat="identity")+scale_fill_manual(values=c("#F3B5B3","#F5E69A","#B5AFDE","#B6D8B6"))+theme_bw(30)+geom_hline(yintercept=1,color="darkgrey",linetype = "dashed",size=4)+ scale_y_continuous(trans = "log2")+xlab("Interaction type") +ylab("Log2 of enrichment")#+theme(axis.text.x = element_text(angle = 45, hjust=1)) 

##

df_wide$dupstatus<-as.factor(df_wide$dupstatus)
df_wide$dupstatus<-factor(df_wide$dupstatus,levels=c("inter_chr_non-ohnolog","ohnolog","intra_chr_non-syntenic","intra_chr_syntenic"))


ggplot(df_wide,aes(x=dupstatus,y=enrichment,fill=dupstatus))+geom_bar(stat="identity")+scale_fill_manual(values=c("firebrick","gold2","darkorchid1","forestgreen"))+theme_bw(base_size = 20)+geom_hline(yintercept=1,color="darkgrey",linetype = "dashed",size=2)+ scale_y_continuous(trans = "log2")

ggplot(df_wide,aes(x=dupstatus,y=enrichment,fill=dupstatus))+geom_bar(stat="identity")+scale_fill_manual(values=c("#F3CDCC","#F4E6AA","#BAB4D8","#C8DBC8"))+theme_bw(base_size = 20)+geom_hline(yintercept=1,color="darkgrey",linetype = "dashed",size=2)+ scale_y_continuous(trans = "log2")+theme(axis.text.x = element_text(angle = 45, hjust=1))

dfwide2<-df_wide
dfwide2$dupstatus<-as.character(dfwide2$dupstatus)
dfwide2[dfwide2$dupstatus=="inter_chr_non-ohnolog","dupstatus"]<-"A"
dfwide2[dfwide2$dupstatus=="ohnolog","dupstatus"]<-"B"
dfwide2[dfwide2$dupstatus=="intra_chr_non-syntenic","dupstatus"]<-"C"
dfwide2[dfwide2$dupstatus=="intra_chr_syntenic","dupstatus"]<-"D"


ggplot(dfwide2,aes(x=dupstatus,y=enrichment,fill=dupstatus))+geom_bar(stat="identity")+scale_fill_manual(values=c("#F3B5B3","#F5E69A","#B5AFDE","#B6D8B6"))+theme_bw(40)+geom_hline(yintercept=1,color="darkgrey",linetype = "dashed",size=4)+ scale_y_continuous(trans = "log2")+xlab("Interaction type") +ylab("Log2 of enrichment") 


##bootstrap of enrichment
shuffle_data<- function(table){
  cols_to_shuffle <- c("gene_chr","gene_start")
  cols_to_shuffle2<-c("snpschrom","snpstart")
  
  # combine the columns you want to shuffle into one matrix
  shuffle_matrix <- table[, cols_to_shuffle]
  shuffle_matrix2 <- table[, cols_to_shuffle2]
  
  
  # shuffle the rows of the matrix
  shuffle_matrix <- shuffle_matrix[sample(nrow(shuffle_matrix)),]
  shuffle_matrix2 <- shuffle_matrix2[sample(nrow(shuffle_matrix2)),]
  
  # split the shuffled matrix back into the original columns
  newtab<- cbind.data.frame(shuffle_matrix, shuffle_matrix2,table[,c("condition","gene","pval","cistrans")])
  return(newtab)
  
}

library(boot)
duplication_table<-read.table(file="https://salmobase.org/datafiles/TSV/synteny/2021-11/AtlanticSalmon/synteny.tsv",header=TRUE,sep="\t")

rna_table<-read.table(file="Synchro_eQTL_5Mwindow.csv",sep=",",header=TRUE)
rna_table$snpschrom<-sapply(strsplit(rna_table$SNP,"_"),FUN = `[[`, 1)
rna_table$snpschrom<-ifelse(substr(rna_table$snpschrom,4,4)==0,substr(rna_table$snpschrom,5,5),substr(rna_table$snpschrom,4,5))
rna_table$snpstart<-as.numeric(as.character(sapply(strsplit(rna_table$SNP,"_"),FUN = `[[`, 2)))
rna_table$gene_chr<-ifelse(substr(rna_table$gene_chr,1,1)=="s",ifelse(substr(rna_table$gene_chr,4,4)==0,substr(rna_table$gene_chr,5,5),substr(rna_table$gene_chr,4,5)),rna_table$gene_chr)
rna_table$snpstart<-as.numeric(as.character(rna_table$snpstart))
rna_table$gene_start<-as.numeric(as.character(rna_table$gene_start))
rna_table<-rna_table[,-c(5,7)]

rna_table_LL<-rna_table[rna_table$condition=="LL",]


genetab<-read.table(file="ss4r_dups_and_singletons_ENSrapid_convPipeline.tsv",header=TRUE)


Boot_function<-function(data,connectiontype){
  
  
  rna_shuffled_LL<-shuffle_data(data)
  
  rna_shuffled_LL_dupstatus<-AddInfoTable(rna_shuffled_LL)
  
  copyofrnashuffle<-rna_shuffled_LL_dupstatus
  
  rna_shuffled_LL_dupstatus[rna_shuffled_LL_dupstatus$dupstatus=="dup_diff_chr","dupstatus"]<-"inter_chr_non-ohnolog"
  rna_shuffled_LL_dupstatus[rna_shuffled_LL_dupstatus$dupstatus=="dup_same_syntenic_region","dupstatus"]<-"intra_chr_syntenic"
  rna_shuffled_LL_dupstatus[rna_shuffled_LL_dupstatus$dupstatus=="dup_same_chr","dupstatus"]<-"intra_chr_non-syntenic"
  rna_shuffled_LL_dupstatus[rna_shuffled_LL_dupstatus$dupstatus=="dup_pair","dupstatus"]<-"ohnolog"
  
  
  
  rna_full_table<-read.table(file="Synchro_uniq_fulltable_withinfo_04.txt",header=TRUE)
  rna_full_table_shuffled<-rna_shuffled_LL_dupstatus
  
  rna_full_table<-rna_full_table[rna_full_table$condition=="LL",]
  rna_full_table<-rna_full_table[,-c(12:14)]
  
  #rna_ohnologs<-rna_full_table[rna_full_table$dupstatus=="ohnolog","gene"]
  #write.table(rna_ohnologs,file="ohnologs_genes.txt",row.names = FALSE,col.names = FALSE,quote=FALSE)
  #rna_all<-rna_full_table[rna_full_table$dupstatus%in%dupstatus_list,"gene"]
  #write.table(rna_all,file="all_genes.txt",row.names = FALSE,col.names = FALSE,quote=FALSE)
  
  rna_full_table$data<-"real"
  rna_full_table_shuffled$data<-"shuffled"
  
  fulltab<-rbind.data.frame(rna_full_table,rna_full_table_shuffled,stringsAsFactors = FALSE)
  
  library(dplyr)
  library(tidyr)
  
  dupstatus_list <- c("inter_chr_non-ohnolog","intra_chr_syntenic","intra_chr_non-syntenic","ohnolog")
  
  dfcount <- fulltab%>% filter(!is.na(cispres)) %>% group_by(dupstatus,data) %>% 
    summarise(total_count=n(),.groups = 'drop') %>%
    as.data.frame() %>% filter(dupstatus %in% dupstatus_list) %>% group_by(data) %>% 
    mutate(prop = total_count / sum(total_count))
  
  df_wide <- pivot_wider(dfcount, id_cols = dupstatus, names_from = data, values_from = total_count)
  df_wide$enrichment<-df_wide$real/df_wide$shuffled
  
  enrichment<-df_wide[df_wide$dupstatus==connectiontype,"enrichment"]
  return(enrichment)
}

#test<-boot(data=rna_table_LL,statistic = Boot_function,R=5,connectiontype="ohnolog")


# Number of bootstrap iterations
n_iterations <- 2


# Run bootstrap using lapply
bootstrap_results <- lapply(1:n_iterations, function(i) {
  Boot_function(rna_table_LL,connectiontype="ohnolog")
})


##bootstrap simplified function

AddInfoTable<-function(table_rna){
  for(y in 1:length(table_rna$snpschrom)){
    val<-TestDuplicate9(y,table_rna)
    table_rna[y,"dupstatus"]<-val[1]
    table_rna[y,"cispres"]<-val[2]
    table_rna[y,"cisduplicate"]<-val[3]
    #table_rna[y,"pairID"]<-val[4]
    #table_rna[y,"samecondition"]<-val[5]
    #table_rna[y,"duplicate_regu_condition"]<-val[6]
  }
  return(table_rna)
}

TestDuplicate9<- function(rownm,table){
  snpchr<-table[rownm,"snpschrom"]
  snpstart<-table[rownm,"snpstart"]
  ##check if the pair considered is trans and regulaed also a cis gene
  if(table[rownm,"cistrans"]=="cis"){
    genetested<-table[rownm,"gene"]
    cispres<-NA
    cis_duplicate<-"no"
  }else{
    ##extract the cis pair for the snp considered
    table_cis_ofpair<-table[table$snpschrom==snpchr & table$snpstart==snpstart & table$cistrans=="cis",]
    if(length(table_cis_ofpair$condition)==0){ ##if none
      cispres<-"no"
      cis_duplicate<-"no"
    }else{ ##if some are found
      cispres<-"yes"
      ## then check if one of the gene cis regulated is the duplicated copy of the trans gene regulated
      genetested<-table[rownm,"gene"]
      for (z in 1:length(table_cis_ofpair$condition)){ ##loop on each cis gene found
        cisgene<-table_cis_ofpair[z,"gene"]
        ##extract the ID of gene pair if both gene are the trans gene tested and the cis gene tested
        matchID1<-genetab[genetab$gene1==genetested & genetab$gene2==cisgene,"ID"]
        matchID2<-genetab[genetab$gene2==genetested & genetab$gene1==cisgene,"ID"]
        if((length(matchID1)!=0 && !is.na(matchID1)) | (length(matchID2)!=0 && !is.na(matchID2))){ ##if a match is found
          cis_duplicate<-"yes"
          break
        }else{ ##and if it's not found
          cis_duplicate<-"no"
        }
      }
    }
  }
  dupstatus<-"none"
  var1<-"no"
  var2<-"no"
  for (i in 1:length(duplication_table$end_x)){ ## test each region of the duplication table
    if((snpstart<=duplication_table[i,"end_x"] && snpstart>=duplication_table[i,"begin_x"] && snpchr==duplication_table[i,"chromosome_x"])){
      if((table[rownm,"gene_start"]<=duplication_table[i,"end_y"] && table[rownm,"gene_start"]>=duplication_table[i,"begin_y"] && table[rownm,"gene_chr"]==duplication_table[i,"chromosome_y"])){
        ## if snp is in part x of the syntenic region and trans gene is in part y
        dupstatus <- "dup_pair"
        break
        ##else if both are in x
      }else if((table[rownm,"gene_start"]<=duplication_table[i,"end_x"] && table[rownm,"gene_start"]>=duplication_table[i,"begin_x"] && table[rownm,"gene_chr"]==duplication_table[i,"chromosome_x"])){
        dupstatus <- "dup_same_syntenic_region"
        break
      }else{
        var1 <- "SNP_dup"
      }
    }else if((snpstart<=duplication_table[i,"end_y"] && snpstart>=duplication_table[i,"begin_y"] && snpchr==duplication_table[i,"chromosome_y"])){
      if((table[rownm,"gene_start"]<=duplication_table[i,"end_x"] && table[rownm,"gene_start"]>=duplication_table[i,"begin_x"] && table[rownm,"gene_chr"]==duplication_table[i,"chromosome_x"])){
        ## if snp is in part y of the syntenic region and trans gene is in part x
        dupstatus <- "dup_pair"
        break
        ##else if both are in y
      }else if((table[rownm,"gene_start"]<=duplication_table[i,"end_y"] && table[rownm,"gene_start"]>=duplication_table[i,"begin_y"] && table[rownm,"gene_chr"]==duplication_table[i,"chromosome_y"])){
        dupstatus <- "dup_same_syntenic_region"
        break
      }else{
        var1 <- "SNP_dup"
      }
    }else if((table[rownm,"gene_start"]<=duplication_table[i,"end_y"] && table[rownm,"gene_start"]>=duplication_table[i,"begin_y"] && table[rownm,"gene_chr"]==duplication_table[i,"chromosome_y"]) | (table[rownm,"gene_start"]<=duplication_table[i,"end_x"] && table[rownm,"gene_start"]>=duplication_table[i,"begin_x"] && table[rownm,"gene_chr"]==duplication_table[i,"chromosome_x"])){
      var2 <- "gene_dup"
    }
  } ##end of the loop
  ##test if both snp and gene are in dup pair or not
  if(var1 == "SNP_dup" && var2 == "gene_dup"){
    ##if in this case, check if the trans gene is on the same chr but not the same syntenic region
    if(snpchr==table[rownm,"gene_chr"]){
      dupstatus <- "dup_same_chr"
    }else{
      dupstatus <- "dup_diff_chr"
    }
  }else if(var1 == "SNP_dup" && var2 != "gene_dup"){
    dupstatus <- "SNP_dup_only"
  }else if(var1 != "SNP_dup" && var2 == "gene_dup"){
    dupstatus <- "gene_dup_only"
  }else if(dupstatus=="no"){
    dupstatus <- "nodup"
  }
  finalvect<-c(dupstatus,cispres,cis_duplicate)
  return(finalvect)
}


##complete table with transcription factors status

rna_full_table<-read.table(file="Synchro_uniq_fulltable_withinfo_04.txt",header=TRUE)

dupstatus_list <- c("inter_chr_non-ohnolog","intra_chr_syntenic","intra_chr_non-syntenic","ohnolog")

library(tidyverse)
library(ggplot2)

##for genes regulated
##data from "general_rna_analysis" script
rna_full_table$geneTF<-ifelse(rna_full_table$gene%in%tfgenes,"TF",ifelse(rna_full_table$gene%in%genewithdesc,"noTF","Unknown"))

dfcount <- rna_full_table %>% group_by(dupstatus,geneTF) %>% 
  summarise(total_count=n(),.groups = 'drop') %>%
  as.data.frame()  %>% filter(dupstatus %in% dupstatus_list)

ggplot(dfcount,aes(x=dupstatus,y=total_count,fill=geneTF))+geom_bar(position = "fill", stat = "identity")+theme_bw(base_size = 20)+ggtitle("Proportion of genes regulated that are TF")

##for SNPs
rna_table_snp<-read.table(file="Synchro_uniq_fulltable_withregulation_snpTF_03.txt",header=TRUE)

dfcount <- rna_table_snp %>% group_by(dupstatus,snpTF) %>% 
  summarise(total_count=n(),.groups = 'drop') %>%
  as.data.frame()  %>% filter(dupstatus %in% dupstatus_list)

ggplot(dfcount,aes(x=dupstatus,y=total_count,fill=snpTF))+geom_bar(position = "fill", stat = "identity")+theme_bw(base_size = 20)+ggtitle("Proportion of SNPs that are in a TF")

sum(dfcount$total_count)

## make the most complete table possible

duplication_table<-read.table(file="https://salmobase.org/datafiles/TSV/synteny/2021-11/AtlanticSalmon/synteny.tsv",header=TRUE,sep="\t")

rna_table<-read.table(file="Synchro_eQTL_5Mwindow.csv",sep=",",header=TRUE)
rna_table$snpschrom<-sapply(strsplit(rna_table$SNP,"_"),FUN = `[[`, 1)
rna_table$snpschrom<-ifelse(substr(rna_table$snpschrom,4,4)==0,substr(rna_table$snpschrom,5,5),substr(rna_table$snpschrom,4,5))
rna_table$snpstart<-as.numeric(as.character(sapply(strsplit(rna_table$SNP,"_"),FUN = `[[`, 2)))
rna_table$gene_chr<-ifelse(substr(rna_table$gene_chr,1,1)=="s",ifelse(substr(rna_table$gene_chr,4,4)==0,substr(rna_table$gene_chr,5,5),substr(rna_table$gene_chr,4,5)),rna_table$gene_chr)
rna_table$snpstart<-as.numeric(as.character(rna_table$snpstart))
rna_table$gene_start<-as.numeric(as.character(rna_table$gene_start))
rna_table<-rna_table[,-c(5,7)]

genetab<-read.table(file="ss4r_dups_and_singletons_ENSrapid_convPipeline.tsv",header=TRUE)

## function to add info

##need the rna table, the duplication table and the gene duplicate table
AddInfoTable<-function(table_rna){
  for(y in 1:length(table_rna$snpschrom)){
    val<-TestDuplicate8(y,table_rna)
    table_rna[y,"dupstatus"]<-val[1]
    table_rna[y,"cispres"]<-val[2]
    table_rna[y,"cisduplicate"]<-val[3]
    table_rna[y,"pairID"]<-val[4]
    table_rna[y,"samecondition"]<-val[5]
    table_rna[y,"duplicate_regu_condition"]<-val[6]
  }
  return(table_rna)
}

#test1<-rna_table[rna_table$snpschrom==29 & rna_table$snpstart==9517949 & rna_table$cistrans=="cis",]
#test2<-rna_table[rna_table$snpschrom==29 & rna_table$snpstart==27489688 & rna_table$cistrans=="cis",]

#length(test2$condition)
#condmatching<-test2[test2$gene=="ENSSSAG00000116217","condition"]

##do a shuffling if needed
## CAREFUL : DO NOT take in account informations others than dupstatus
#subtable<-rna_table[rna_table$cistrans=="trans",]
#rna_table_trans_suffled<-shuffle_data(subtable)
#rna_table_shuffled<-AddInfoTable(rna_table_trans_suffled)
#rna_table_shuffled[rna_table_shuffled$dupstatus=="dup_diff_chr","dupstatus"]<-"inter_chr_non-ohnolog"
#rna_table_shuffled[rna_table_shuffled$dupstatus=="dup_same_syntenic_region","dupstatus"]<-"intra_chr_syntenic"
#rna_table_shuffled[rna_table_shuffled$dupstatus=="dup_same_chr","dupstatus"]<-"intra_chr_non-syntenic"
#rna_table_shuffled[rna_table_shuffled$dupstatus=="dup_pair","dupstatus"]<-"ohnolog"
#write.table(rna_table_shuffled,file="Synchro_uniq_fulltable_withinfo_shuffled_02.txt",row.names=FALSE)


#subtable<-rna_table[11700:12000,]

#subtable<-AddInfoTable(subtable)

rna_table<-AddInfoTable(rna_table)


#write.table(rna_table,file="Synchro_uniq_fulltable_withinfo_01.txt",row.names=FALSE)
write.table(rna_table,file="Synchro_uniq_fulltable_withinfo_02.txt",row.names=FALSE)



TestDuplicate8<- function(rownm,table){
  snpchr<-table[rownm,"snpschrom"]
  snpstart<-table[rownm,"snpstart"]
  ##check if the pair considered is trans and regulaed also a cis gene
  if(table[rownm,"cistrans"]=="cis"){
    ##check if they have a duplicate counterpart
    genetested<-table[rownm,"gene"]
    matchID1<-genetab[genetab$gene1==genetested,"ID"]
    matchtype1<-genetab[genetab$gene1==genetested,"type"]
    matchID2<-genetab[genetab$gene2==genetested,"ID"]
    matchtype2<-genetab[genetab$gene2==genetested,"type"]
    if((length(matchID1)!=0 && !is.na(matchID1))){
      if(matchtype1=="ss4r"){
        pairID<-matchID1
      }else{
        pairID<-NA
      }
    }else if((length(matchID2)!=0 && !is.na(matchID2))){
      if(matchtype2=="ss4r"){
        pairID<-matchID2
      }else{
        pairID<-NA
      }
    }else{
      pairID<-NA
    }
    cispres<-NA
    cis_duplicate<-"no"
    samecondition<-NA
    duplicate_regu_condition<-NA
  }else{
    ##extract the cis pair for the snp considered
    table_cis_ofpair<-table[table$snpschrom==snpchr & table$snpstart==snpstart & table$cistrans=="cis",]
    if(length(table_cis_ofpair$condition)==0){ ##if none
      cispres<-"no"
      cis_duplicate<-"no"
      pairID<-NA
      samecondition<-NA
      duplicate_regu_condition<-NA
    }else{ ##if some are found
      cispres<-"yes"
      ## then check if one of the gene cis regulated is the duplicated copy of the trans gene regulated
      genetested<-table[rownm,"gene"]
      for (z in 1:length(table_cis_ofpair$condition)){ ##loop on each cis gene found
        cisgene<-table_cis_ofpair[z,"gene"]
        ##extract the ID of gene pair if both gene are the trans gene tested and the cis gene tested
        matchID1<-genetab[genetab$gene1==genetested & genetab$gene2==cisgene,"ID"]
        matchID2<-genetab[genetab$gene2==genetested & genetab$gene1==cisgene,"ID"]
        if((length(matchID1)!=0 && !is.na(matchID1)) | (length(matchID2)!=0 && !is.na(matchID2))){ ##if a match is found
          cis_duplicate<-"yes"
          if((length(matchID1)!=0 && !is.na(matchID1))){
            pairID<-matchID1
          }else if((length(matchID2)!=0 && !is.na(matchID2))){
            pairID<-matchID2
          }
          ##get conditions
          condmatching<-table_cis_ofpair[table_cis_ofpair$gene==cisgene,"condition"]
          condition_gene_tested<-table[rownm,"condition"]
          if(condition_gene_tested%in%condmatching){ ##if the condition of the trans gene is also found in the cis dupicate gene
            duplicate_regu_condition<-condition_gene_tested
            samecondition<-"yes"
          }else{ ##or else
            samecondition<-"no"
            if(length(condmatching)>1){
              duplicate_regu_condition<-paste(condmatching[1],condmatching[2],sep=";")
            }else{
              duplicate_regu_condition<-condmatching
            }
          }
          ##get out of the loop
          break
        }else{ ##and if it's not found
          cis_duplicate<-"no"
          ##check if there is a counterpart gene or not
          genetested<-table[rownm,"gene"]
          matchID1<-genetab[genetab$gene1==genetested,"ID"]
          matchtype1<-genetab[genetab$gene1==genetested,"type"]
          matchID2<-genetab[genetab$gene2==genetested,"ID"]
          matchtype2<-genetab[genetab$gene2==genetested,"type"]
          if((length(matchID1)!=0 && !is.na(matchID1))){
            if(matchtype1=="ss4r"){
              pairID<-matchID2
            }else{
              pairID<-NA
            }
          }else if((length(matchID2)!=0 && !is.na(matchID2))){
            if(matchtype2=="ss4r"){
              pairID<-matchID2
            }else{
              pairID<-NA
            }
          }else{
            pairID<-NA
          }
          samecondition<-NA
          duplicate_regu_condition<-NA
        }
      }
    }
  }
  dupstatus<-"none"
  var1<-"no"
  var2<-"no"
  for (i in 1:length(duplication_table$end_x)){ ## test each region of the duplication table
    if((snpstart<=duplication_table[i,"end_x"] && snpstart>=duplication_table[i,"begin_x"] && snpchr==duplication_table[i,"chromosome_x"])){
      if((table[rownm,"gene_start"]<=duplication_table[i,"end_y"] && table[rownm,"gene_start"]>=duplication_table[i,"begin_y"] && table[rownm,"gene_chr"]==duplication_table[i,"chromosome_y"])){
        ## if snp is in part x of the syntenic region and trans gene is in part y
        dupstatus <- "dup_pair"
        break
        ##else if both are in x
      }else if((table[rownm,"gene_start"]<=duplication_table[i,"end_x"] && table[rownm,"gene_start"]>=duplication_table[i,"begin_x"] && table[rownm,"gene_chr"]==duplication_table[i,"chromosome_x"])){
        dupstatus <- "dup_same_syntenic_region"
        break
      }else{
        var1 <- "SNP_dup"
      }
    }else if((snpstart<=duplication_table[i,"end_y"] && snpstart>=duplication_table[i,"begin_y"] && snpchr==duplication_table[i,"chromosome_y"])){
      if((table[rownm,"gene_start"]<=duplication_table[i,"end_x"] && table[rownm,"gene_start"]>=duplication_table[i,"begin_x"] && table[rownm,"gene_chr"]==duplication_table[i,"chromosome_x"])){
        ## if snp is in part y of the syntenic region and trans gene is in part x
        dupstatus <- "dup_pair"
        break
        ##else if both are in y
      }else if((table[rownm,"gene_start"]<=duplication_table[i,"end_y"] && table[rownm,"gene_start"]>=duplication_table[i,"begin_y"] && table[rownm,"gene_chr"]==duplication_table[i,"chromosome_y"])){
        dupstatus <- "dup_same_syntenic_region"
        break
      }else{
        var1 <- "SNP_dup"
      }
    }else if((table[rownm,"gene_start"]<=duplication_table[i,"end_y"] && table[rownm,"gene_start"]>=duplication_table[i,"begin_y"] && table[rownm,"gene_chr"]==duplication_table[i,"chromosome_y"]) | (table[rownm,"gene_start"]<=duplication_table[i,"end_x"] && table[rownm,"gene_start"]>=duplication_table[i,"begin_x"] && table[rownm,"gene_chr"]==duplication_table[i,"chromosome_x"])){
      var2 <- "gene_dup"
    }
  } ##end of the loop
  ##test if both snp and gene are in dup pair or not
  if(var1 == "SNP_dup" && var2 == "gene_dup"){
    ##if in this case, check if the trans gene is on the same chr but not the same syntenic region
    if(snpchr==table[rownm,"gene_chr"]){
      dupstatus <- "dup_same_chr"
    }else{
      dupstatus <- "dup_diff_chr"
    }
  }else if(var1 == "SNP_dup" && var2 != "gene_dup"){
    dupstatus <- "SNP_dup_only"
  }else if(var1 != "SNP_dup" && var2 == "gene_dup"){
    dupstatus <- "gene_dup_only"
  }else if(dupstatus=="no"){
    dupstatus <- "nodup"
  }
  finalvect<-c(dupstatus,cispres,cis_duplicate,pairID,samecondition,duplicate_regu_condition)
  return(finalvect)
}





##debug
rna_table[11747,]

rownm<-140
table<-subtable

snpchr<-table[rownm,"snpschrom"]
snpstart<-table[rownm,"snpstart"]
##check if the pair considered is trans and regulaed also a cis gene
if(table[rownm,"cistrans"]=="cis"){
  cispres<-NA
  cis_duplicate<-"no"
  pairID<-NA
  samecondition<-NA
  duplicate_regu_condition<-NA
}else{
  ##extract the cis pair for the snp considered
  table_cis_ofpair<-table[table$snpschrom==snpchr & table$snpstart==snpstart & table$cistrans=="cis",]
  if(length(table_cis_ofpair$condition)==0){ ##if none
    cispres<-"no"
    cis_duplicate<-"no"
    pairID<-NA
    samecondition<-NA
    duplicate_regu_condition<-NA
  }else{ ##if some are found
    cispres<-"yes"
    ## then check if one of the gene cis regulated is the duplicated copy of the trans gene regulated
    genetested<-table[rownm,"gene"]
    for (z in 1:length(table_cis_ofpair$condition)){ ##loop on each cis gene found
      cisgene<-table_cis_ofpair[z,"gene"]
      ##extract the ID of gene pair if both gene are the trans gene tested and the cis gene tested
      matchID<-genetab[(genetab$gene1==genetested & genetab$gene2==cisgene) | (genetab$gene2==genetested & genetab$gene1==cisgene),"ID"]
      if(length(matchID)==0){ ##if a match is not found
        cis_duplicate<-"no"
        pairID<-NA
        samecondition<-NA
        duplicate_regu_condition<-NA
      }else{ ##and if it's found
        cis_duplicate<-"yes"
        pairID<-matchID
        ##get conditions
        condmatching<-table_cis_ofpair[table_cis_ofpair$gene==cisgene,"condition"]
        condition_gene_tested<-table[rownm,"condition"]
        if(condition_gene_tested%in%condmatching){ ##if the condition of the trans gene is also found in the cis dupicate gene
          duplicate_regu_condition<-condition_gene_tested
          samecondition<-"yes"
        }else{ ##or else
          samecondition<-"no"
          if(length(condmatching)>1){
            duplicate_regu_condition<-paste(condmatching[1],condmatching[2],sep=";")
          }else{
            duplicate_regu_condition<-condmatching
          }
        }
        ##get out of the loop
        break
      }
    }
  }
}


#test<-rna_table[1:200,]

shuffle_data<- function(table){
  cols_to_shuffle <- c("gene_chr","gene_start")
  cols_to_shuffle2<-c("snpschrom","snpstart")
  
  # combine the columns you want to shuffle into one matrix
  shuffle_matrix <- table[, cols_to_shuffle]
  shuffle_matrix2 <- table[, cols_to_shuffle2]
  
  
  # shuffle the rows of the matrix
  shuffle_matrix <- shuffle_matrix[sample(nrow(shuffle_matrix)),]
  shuffle_matrix2 <- shuffle_matrix2[sample(nrow(shuffle_matrix2)),]
  
  # split the shuffled matrix back into the original columns
  newtab<- cbind.data.frame(shuffle_matrix, shuffle_matrix2,table[,c("condition","gene","pval","cistrans")])
  return(newtab)
  
}
#rna_table_trans_suffled<-shuffle_data(test)
rna_table_LL<-rna_table[rna_table$condition=="LL",]
rna_shuffled_LL<-shuffle_data(rna_table_LL)

rna_shuffled_LL_dupstatus<-AddInfoTable(rna_shuffled_LL)

copyofrnashuffle<-rna_shuffled_LL_dupstatus

rna_shuffled_LL_dupstatus[rna_shuffled_LL_dupstatus$dupstatus=="dup_diff_chr","dupstatus"]<-"inter_chr_non-ohnolog"
rna_shuffled_LL_dupstatus[rna_shuffled_LL_dupstatus$dupstatus=="dup_same_syntenic_region","dupstatus"]<-"intra_chr_syntenic"
rna_shuffled_LL_dupstatus[rna_shuffled_LL_dupstatus$dupstatus=="dup_same_chr","dupstatus"]<-"intra_chr_non-syntenic"
rna_shuffled_LL_dupstatus[rna_shuffled_LL_dupstatus$dupstatus=="dup_pair","dupstatus"]<-"ohnolog"

write.table(rna_shuffled_LL_dupstatus,file="Synchro_uniq_fulltable_withinfo_shuffled_03.txt",row.names = FALSE)

### Get the duplicate of each gene and see if it is regulated

genetabss<-genetab[genetab$type=="ss4r",]
for ( i in 1:length(rna_full_table$gene)){
  gene_tested<-rna_full_table[i,"gene"]
  ##get duplicated gene
  subtab1<-genetabss[(genetabss$gene1==gene_tested),"gene2"]
  subtab2<-genetabss[(genetabss$gene2==gene_tested),"gene1"]
  if((length(subtab1)==0 | is.null(subtab1)) && (length(subtab2)==0 | is.null(subtab2)) ){
    dupID<-NA
  }else if((length(subtab1)==0 | is.null(subtab1)) && (length(subtab2)!=0 & !is.null(subtab2)) ){
    dupID<-subtab2
  }else if((length(subtab1)!=0 & !is.null(subtab1)) && (length(subtab2)==0 | is.null(subtab2)) ){
    dupID<-subtab1
  }
  ##check if the gene is present and regulated
  if(is.na(dupID)){
    regulation<-NA
  }else{
    if(dupID%in%rna_full_table$gene){
      subtabID<-rna_full_table[rna_full_table$gene==dupID,]
      for (y in 1:length(subtabID$gene)){
        chromy<-subtabID[y,"snpschrom"]
        starty<-subtabID[y,"snpstart"]
        chromi<-rna_full_table[i,"snpschrom"]
        starti<-rna_full_table[i,"snpstart"]
        if(chromy==chromi && starty==starti){
          regulation<-"samesnp"
          break
        }else{
          regulation<-"notsamesnp"
        }
      }
    }else{
      regulation<-"notregulated"
    }
  }
  rna_full_table[i,"dupID"]<-dupID
  rna_full_table[i,"regulation"]<-regulation
}

write.table(rna_full_table,file="Synchro_uniq_fulltable_withregulation_03.txt",row.names=FALSE)
rna_full_table_trimmed<-rna_full_table[!is.na(rna_full_table$regulation),]


library(dplyr)
library(ggplot2)

rna_full_table<-read.table(file="Synchro_uniq_fulltable_withregulation_03.txt",header=TRUE)
rna_full_table_trimmed<-rna_full_table[!is.na(rna_full_table$regulation),]

dupstatus_list <- c("inter_chr_non-ohnolog","intra_chr_syntenic","intra_chr_non-syntenic","ohnolog")


dfcount <- rna_full_table_trimmed %>% group_by(dupstatus,regulation,cistrans) %>% 
  summarise(total_count=n(),.groups = 'drop') %>%
  as.data.frame() 

ciscountnotregu<-sum(dfcount[dfcount$cistrans=="cis" & dfcount$regulation=="notregulated","total_count"])
ciscountnotsamesnp<-sum(dfcount[dfcount$cistrans=="cis" & dfcount$regulation=="notsamesnp","total_count"])
ciscountsamesnp<-sum(dfcount[dfcount$cistrans=="cis" & dfcount$regulation=="samesnp","total_count"])
cisdata<-data.frame(dupstatus="cis",regulation=c("notregulated","notsamesnp","samesnp"),total_count=c(ciscountnotregu,ciscountnotsamesnp,ciscountsamesnp))

dfcount<- dfcount %>% filter(dupstatus %in% dupstatus_list)
dfcount2<-dfcount[dfcount$cistrans=="trans",c(1,2,4)]
dfcount3<-rbind.data.frame(dfcount2,cisdata,stringsAsFactors = FALSE)

dfcount3$dupstatus<-as.factor(dfcount3$dupstatus)
dfcount3$dupstatus<-factor(dfcount3$dupstatus,levels=c("ohnolog","inter_chr_non-ohnolog","intra_chr_non-syntenic","intra_chr_syntenic","cis"))


ggplot(dfcount3,aes(x=dupstatus,y=total_count,fill=regulation))+geom_bar(position="fill",stat="identity")+ggtitle("Regulation of the counterpart duplicate gene")+theme_bw(base_size = 20)


##see if regulated gene by same snp are TF
#rna_full_table_trimmed$mergedvar<-paste(rna_full_table_trimmed$dupstatus,rna_full_table_trimmed$regulation,sep="_")

dfcount <- rna_full_table_trimmed %>% group_by(dupstatus,regulation,geneTF,cistrans) %>% 
  summarise(total_count=n(),.groups = 'drop') %>%
  as.data.frame() 

cistab<-dfcount[dfcount$cistrans=="cis",]
dfcountcis <- cistab %>% group_by(regulation,geneTF) %>% 
  summarise(total_count=n(),.groups = 'drop') %>%
  as.data.frame() 
dfcountcis$dupstatus<-"cis"

dfcount2<-dfcount[dfcount$cistrans=="trans",c(-4)]
dfcount3<- dfcount2 %>% filter(dupstatus %in% dupstatus_list)
dfcount4<-rbind.data.frame(dfcount3,dfcountcis,stringsAsFactors = FALSE)

dfcount4$dupstatus<-as.factor(dfcount4$dupstatus)
dfcount4$dupstatus<-factor(dfcount4$dupstatus,levels=c("ohnolog","inter_chr_non-ohnolog","intra_chr_non-syntenic","intra_chr_syntenic","cis"))

dfcount4$mergedvar<-paste(dfcount4$dupstatus,dfcount4$regulation,sep="_")

library(gridExtra)
p<-ggplot(dfcount4,aes(x=mergedvar,y=total_count,fill=geneTF))+geom_bar(position="fill",stat="identity")

ggplot(dfcount4,aes(x=regulation,y=total_count,fill=geneTF))+geom_bar(position="fill",stat="identity")+facet_wrap(~dupstatus)+theme_bw(base_size = 20)+ggtitle("proportion of TF gene regulated by snp for each type of regulation of the duplicated gene")+theme(axis.text.x = element_text(angle = 45, hjust=1,size=16))


p2<-ggplot(dfcount4,aes(x=mergedvar,y=total_count))+geom_bar(stat="identity")
grid.arrange(p2, p,ncol=1, nrow=2)







## make a table of duplicate gene pairs and others informations

genetab<-read.table(file="ss4r_dups_and_singletons_ENSrapid_convPipeline.tsv",header=TRUE)
genetab<-genetab[genetab$type=="ss4r",]
genetab<-genetab[,-c(1,2,6)]

rna_full_table<-read.table(file="Synchro_uniq_fulltable_withinfo_04.txt",header=TRUE)
rna_full_tab_LL<-rna_full_table[rna_full_table$condition=="LL",]

## check different scenarios:
## no regulator: genes are not found in the rna_full_table
## one gene of the pair is regulated not the other
## both are regulated but by different snps
## both are regulated by the same snp
## note : a pair can be the 2 lasts

#gene1row1<-genetab[180,"gene1"]
#tabpair1<-rna_full_tab_LL[rna_full_tab_LL$gene==gene1row1,]
#gene2row1<-genetab[180,"gene2"]
#tabpair2<-rna_full_tab_LL[rna_full_tab_LL$gene==gene2row1,]

for(y in 1:length(genetab$ID)){
  
  gene1row1<-genetab[y,"gene1"]
  tabpair1<-rna_full_tab_LL[rna_full_tab_LL$gene==gene1row1,]
  gene2row1<-genetab[y,"gene2"]
  tabpair2<-rna_full_tab_LL[rna_full_tab_LL$gene==gene2row1,]  
  
  ltab1<-length(tabpair1$condition)
  ltab2<-length(tabpair2$condition)
  
  noreg<-"regulated"
  singreg<-"notsingreg"
  bothdiffreg<-"nobothdiffreg"
  bothsamereg<-"nobothsamereg"
  if(ltab1==0 && ltab2==0){
    noreg<-"noreg"
  }
  if(ltab1==0 && ltab2!=0){
    singreg<-"singreg"
  }
  if(ltab1!=0 && ltab2==0){
    singreg<-"singreg"
  }
  if(ltab1!=0 && ltab2!=0){
    tabpair1$IDreg<-paste(tabpair1$snpschrom,tabpair1$snpstart,sep="_")
    tabpair2$IDreg<-paste(tabpair2$snpschrom,tabpair2$snpstart,sep="_")
    issamesnp<-ifelse(tabpair1$IDreg%in%tabpair2$IDreg,"yes","no")
    if("no"%in%issamesnp){
      bothdiffreg<-"bothdiffreg"
    }
    if("yes"%in%issamesnp){
      bothsamereg<-"bothsamereg"
    }
  }
  
  genetab[y,"noreg"]<-noreg
  genetab[y,"singreg"]<-singreg
  ##say that if regulated by same snp then we keep only this result
  if(bothdiffreg=="bothdiffreg" && bothsamereg=="bothsamereg"){
    bothdiffreg<-"nobothdiffreg"
  }
  genetab[y,"bothdiffreg"]<-bothdiffreg
  genetab[y,"bothsamereg"]<-bothsamereg
}

write.table(genetab,file="genetabgenepair01.txt",row.names=FALSE)



library(ggplot2)
library(dplyr)

##make proportion table

genetab<-read.table(file="genetabgenepair01.txt",header=TRUE)

regprop<-length(genetab[genetab$noreg=="regulated","noreg"])/length(genetab$gene1)
singregprop<-length(genetab[genetab$singreg=="singreg","noreg"])/length(genetab$gene1)
bothdiffregprop<-length(genetab[genetab$bothdiffreg=="bothdiffreg","noreg"])/length(genetab$gene1)
bothsameregprop<-length(genetab[genetab$bothsamereg=="bothsamereg","noreg"])/length(genetab$gene1)

proptable<-data.frame(categorie=c("regulated","one copy regulated","both copies regulated by different SNPs","both copies by the same SNPs"),prop=c(regprop,singregprop,bothdiffregprop,bothsameregprop))
ggplot(proptable,aes(x=categorie,y=prop))+geom_bar(stat = "identity")


regprop<-length(genetab[genetab$noreg=="noreg","noreg"])/length(genetab$gene1)  ##52.1%
singregprop<-length(genetab[genetab$singreg=="singreg","noreg"])/length(genetab$gene1)   ##30.7%
bothdiffregprop<-length(genetab[genetab$bothdiffreg=="bothdiffreg","noreg"])/length(genetab$gene1) ##16.3%
bothsameregprop<-length(genetab[genetab$bothsamereg=="bothsamereg","noreg"])/length(genetab$gene1)  ##0.8%

proptable<-data.frame(data="all",group=c("No eQTL","Single copy eQTL","Distinct eQTL","Shared eQTL"),prop=c(regprop,singregprop,bothdiffregprop,bothsameregprop))
ggplot(proptable,aes(x=data,y=prop,fill=group))+geom_bar(stat = "identity",position="fill")+theme_bw(base_size = 40)+scale_fill_manual(values=c("#F9B27C","#99B8EA","#8CE86D","#D186D8"))



##another version
for(y in 1:length(genetab$ID)){
  
  gene1row1<-genetab[y,"gene1"]
  tabpair1<-rna_full_tab_LL[rna_full_tab_LL$gene==gene1row1,]
  gene2row1<-genetab[y,"gene2"]
  tabpair2<-rna_full_tab_LL[rna_full_tab_LL$gene==gene2row1,]  
  
  ltab1<-length(tabpair1$condition)
  ltab2<-length(tabpair2$condition)
  
  reg<-"noreg"
  bothdiffreg<-"nobothdiffreg"
  bothsamereg<-"nobothsamereg"
  if(ltab1==0 && ltab2==0){
  }
  if(ltab1==0 && ltab2!=0){
    reg<-"singreg"
  }
  if(ltab1!=0 && ltab2==0){
    reg<-"singreg"
  }
  if(ltab1!=0 && ltab2!=0){
    tabpair1$IDreg<-paste(tabpair1$snpschrom,tabpair1$snpstart,sep="_")
    tabpair2$IDreg<-paste(tabpair2$snpschrom,tabpair2$snpstart,sep="_")
    issamesnp<-ifelse(tabpair1$IDreg%in%tabpair2$IDreg,"yes","no")
    if("no"%in%issamesnp){
      bothdiffreg<-"bothdiffreg"
    }
    if("yes"%in%issamesnp){
      bothsamereg<-"bothsamereg"
    }
  }
  
  
  
  ##say that if regulated by same snp then we keep only this result
  if(bothdiffreg=="bothdiffreg" && bothsamereg=="bothsamereg"){
    reg<-"bothsamereg"
  }else if(bothdiffreg=="nobothdiffreg" && bothsamereg=="bothsamereg"){
    reg<-"bothsamereg"
  }else if(bothdiffreg=="bothdiffreg" && bothsamereg=="nobothsamereg"){
    reg<-"bothdiffreg"
  }
  genetab[y,"regulation"]<-reg
}

write.table(genetab,file="genetabgenepair02.txt",row.names=FALSE)




##now include differentiate cis and trans
## check different scenarios:
## nb of gene regulated : none / one / both not shared / both shared
## if one, regulation is cis or trans or both (and if trans, is it in which of the 4 scenario are we)
## if both are regulated : is it cis + trans or trans + trans (and which type of trans)
## in the case of 2 trans, paste both trans type
## see drawing for expectations  

genetab<-read.table(file="ss4r_dups_and_singletons_ENSrapid_convPipeline.tsv",header=TRUE)
genetab<-genetab[genetab$type=="ss4r",]
genetab<-genetab[,-c(1,2,6)]

rna_full_table<-read.table(file="Synchro_uniq_fulltable_withinfo_04.txt",header=TRUE)
rna_full_tab_LL<-rna_full_table[rna_full_table$condition=="LL",]

for(y in 1:length(genetab$ID)){
  
  gene1row1<-genetab[y,"gene1"]
  tabpair1<-rna_full_tab_LL[rna_full_tab_LL$gene==gene1row1,]
  gene2row1<-genetab[y,"gene2"]
  tabpair2<-rna_full_tab_LL[rna_full_tab_LL$gene==gene2row1,] 
  
  ltab1<-length(tabpair1$condition)
  ltab2<-length(tabpair2$condition)
  
  reg<-"noreg"
  regtype<-"none"
  regscenario<-"none"
  bothdiffreg<-"nobothdiffreg"
  bothsamereg<-"nobothsamereg"
  regscenario2<-""
  regscenario1<-""
  if(ltab1==0 && ltab2==0){
  }
  if(ltab1==0 && ltab2!=0){
    reg<-"singreg"
    if("cis"%in%tabpair2$cistrans){
      if("trans"%in%tabpair2$cistrans){
        regtype<-"both"
        tabpair2trans<-tabpair2[tabpair2$cistrans=="trans","dupstatus"]
        regscenario <- paste(tabpair2trans, collapse = "|")
      }else{
        regtype<-"cis"
      }
    }else if("trans"%in%tabpair2$cistrans){
      regtype<-"trans"
      tabpair2trans<-tabpair2[tabpair2$cistrans=="trans","dupstatus"]
      regscenario <- paste(tabpair2trans, collapse = "|")
    }
  }
  if(ltab1!=0 && ltab2==0){
    reg<-"singreg"
    if("cis"%in%tabpair1$cistrans){
      if("trans"%in%tabpair1$cistrans){
        regtype<-"both"
        tabpair1trans<-tabpair1[tabpair1$cistrans=="trans","dupstatus"]
        regscenario <- paste(tabpair1trans, collapse = "|")
      }else{
        regtype<-"cis"
      }
    }else if("trans"%in%tabpair1$cistrans){
      regtype<-"trans"
      tabpair1trans<-tabpair1[tabpair1$cistrans=="trans","dupstatus"]
      regscenario <- paste(tabpair1trans, collapse = "|")
    }
  }
  if(ltab1!=0 && ltab2!=0){
    tabpair1$IDreg<-paste(tabpair1$snpschrom,tabpair1$snpstart,sep="_")
    tabpair2$IDreg<-paste(tabpair2$snpschrom,tabpair2$snpstart,sep="_")
    issamesnp<-ifelse(tabpair1$IDreg%in%tabpair2$IDreg,"yes","no")
    if("no"%in%issamesnp){
      bothdiffreg<-"bothdiffreg"
    }
    if("yes"%in%issamesnp){
      bothsamereg<-"bothsamereg"
    }
    
    ##additional information
    if("cis"%in%tabpair1$cistrans){
      if("trans"%in%tabpair1$cistrans){
        regtype1<-"both"
        tabpair1trans<-tabpair1[tabpair1$cistrans=="trans","dupstatus"]
        regscenario1 <- paste(tabpair1trans, collapse = "|")
      }else{
        regtype1<-"cis"
      }
    }else if("trans"%in%tabpair1$cistrans){
      regtype1<-"trans"
      tabpair1trans<-tabpair1[tabpair1$cistrans=="trans","dupstatus"]
      regscenario1 <- paste(tabpair1trans, collapse = "|")
    }
    
    if("cis"%in%tabpair2$cistrans){
      if("trans"%in%tabpair2$cistrans){
        regtype2<-"both"
        tabpair2trans<-tabpair2[tabpair2$cistrans=="trans","dupstatus"]
        regscenario2 <- paste(tabpair2trans, collapse = "|")
      }else{
        regtype2<-"cis"
      }
    }else if("trans"%in%tabpair2$cistrans){
      regtype2<-"trans"
      tabpair2trans<-tabpair2[tabpair2$cistrans=="trans","dupstatus"]
      regscenario2 <- paste(tabpair2trans, collapse = "|")
    }
    
    regtype<-paste(regtype1,regtype2,sep="_")
    regscenario<-paste(regscenario1,regscenario2,sep="#")
    
  }
  
  
  
  ##say that if regulated by same snp then we keep only this result
  if(bothdiffreg=="bothdiffreg" && bothsamereg=="bothsamereg"){
    reg<-"bothsamereg"
  }else if(bothdiffreg=="nobothdiffreg" && bothsamereg=="bothsamereg"){
    reg<-"bothsamereg"
  }else if(bothdiffreg=="bothdiffreg" && bothsamereg=="nobothsamereg"){
    reg<-"bothdiffreg"
  }
  genetab[y,"regulation"]<-reg
  genetab[y,"regtype"]<-regtype
  genetab[y,"reegscenario"]<-regscenario
}

#write.table(genetab,file="genetabgenepairwithtrans01.txt",row.names=FALSE)

genetab<-read.table(file="genetabgenepairwithtrans01.txt",header=TRUE)

genetabtest<-genetab[genetab$regulation=="bothsamereg",]


genetabtest[genetabtest$reegscenario=="intra_chr_syntenic#ohnolog" | genetabtest$reegscenario=="ohnolog#intra_chr_syntenic",]
genetabtest[genetabtest$reegscenario=="inter_chr_non-ohnolog#intra_chr_non-syntenic" | genetabtest$reegscenario=="intra_chr_non-syntenic#inter_chr_non-ohnolog",]
genetabtest[genetabtest$reegscenario=="gene_dup_only#gene_dup_only","regtype"]

test<-genetabtest[genetabtest$reegscenario=="inter_chr_non-ohnolog#inter_chr_non-ohnolog","regtype"]
table(test)
#genetabtest2<-genetab[genetab$regulation=="bothdiffreg",]
#genetabtest3<-genetab[genetab$regulation=="bothdiffreg" & genetab$regtype=="cis_cis",]


effect<-c(6,2,5,9,67)
scenario<-c("A","E","D","C","B")
scenario_detailled<-c("A1","A2","A3","B1","B2","B3","C1","C2","C3","D1","D2","E1","E2")
effect_detailled<-c(1,4,1,41,19,7,2,3,4,3,2,1,1)

tableglobal<-data.frame(effect=effect,scenario=scenario,data="global")

ggplot(tableglobal,aes(x=data,y=effect,fill=scenario))+geom_bar(position="fill",stat="identity")+geom_text(aes(label = effect), position = position_fill(vjust = 0.5),size=10,color="black",fontface = "bold")+theme_bw(20)+scale_fill_manual(values=c("firebrick1","dodgerblue3","olivedrab2","darkorchid1","gold1"))

tabledetailled<-data.frame(effect=effect_detailled,scenario=scenario_detailled,data="detailled")

ggplot(tabledetailled,aes(x=data,y=effect,fill=scenario))+geom_bar(position="fill",stat="identity")+geom_text(aes(label = effect), position = position_fill(vjust = 0.5),size=6,color="lemonchiffon",fontface = "bold")+theme_bw(20)+scale_fill_manual(values=c("brown1","red","firebrick4","deepskyblue1","dodgerblue4","blue","olivedrab3","chartreuse2","darkolivegreen4","darkorchid1","darkorchid4","yellow3","goldenrod3"))





### new data // same analysis
newdata<-read.csv(file="cistrans_gene_August.csv",header=TRUE)

newdata$snpchrom<-as.numeric(as.character(ifelse(substr(sapply(strsplit(newdata$lead_SNP,"_"),FUN = `[[`, 1),4,4)==0,substr(sapply(strsplit(newdata$lead_SNP,"_"),FUN = `[[`, 1),5,5),substr(sapply(strsplit(newdata$lead_SNP,"_"),FUN = `[[`, 1),4,5))))
newdata$snppos<-as.numeric(as.character(sapply(strsplit(newdata$lead_SNP,"_"),FUN = `[[`, 2)))

duplication_table<-read.table(file="https://salmobase.org/datafiles/TSV/synteny/2021-11/AtlanticSalmon/synteny.tsv",header=TRUE,sep="\t")




## convert position to number

lengthvector<-c(174498729,95481959,105780080,90536438,92788608,96060288,68862998,28860523,161282225,125877811,111868677,101677876,
                114417674,101980477,110670232,96486271,87489397,84084598,88107222,96847506,59819933,63823863,52460201,49354470,54385492,
                55994222,45305548,41468476,43051128)
chromadd<-c(0,cumsum(lengthvector)[1:28])

## convert for duplication table
chromname<-c(1:29)
chrommid<-chromadd+lengthvector/2
listval<-cbind.data.frame(chromname,lengthvector,chromadd,chrommid)

duplication_table$startxnb<-duplication_table$begin_x+listval[match(duplication_table$chromosome_x,listval$chromname),3]

duplication_table$endxnb<-duplication_table$end_x+listval[match(duplication_table$chromosome_x,listval$chromname),3]

duplication_table$startynb<-duplication_table$begin_y+listval[match(duplication_table$chromosome_y,listval$chromname),3]

duplication_table$endynb<-duplication_table$end_y+listval[match(duplication_table$chromosome_y,listval$chromname),3]

##and for the data

newdata$snpposnb<-newdata$snppos+listval[match(newdata$snpchrom,listval$chromname),3]

newdata$genestartnb<-newdata$start+listval[match(newdata$chr,listval$chromname),3]

newdata$geneendnb<-newdata$end+listval[match(newdata$chr,listval$chromname),3]


####end of transformation
#duplication_tablerow1<-duplication_table[1,]

#test2<-ifelse(newdata$snpposnb>=duplication_tablerow1$startxnb & newdata$snpposnb<=duplication_tablerow1$endxnb,ifelse(newdata$genestartnb>=duplication_tablerow1$startynb & newdata$genestartnb<=duplication_tablerow1$endynb,"Inter-chromosomal ohnolog",ifelse(newdata$geneendnb>=duplication_tablerow1$startynb & newdata$geneendnb<=duplication_tablerow1$endynb,"Inter-chromosomal ohnolog",ifelse(newdata$genestartnb>=chromadd[newdata$snpchrom] & newdata$genestartnb<=chromadd[newdata$snpchrom+1],ifelse(newdata$genestartnb>=duplication_tablerow1$startxnb & newdata$genestartnb<=duplication_tablerow1$endxnb,"Intra-chromosomal syntenic",ifelse(newdata$geneendnb>=duplication_tablerow1$startxnb & newdata$geneendnb<=duplication_tablerow1$endxnb,"Intra-chromosomal syntenic","Intra-chromosomal non syntenic")),"eQTL present"))),ifelse(newdata$snpposnb>=duplication_tablerow1$startynb & newdata$snpposnb<=duplication_tablerow1$endynb,ifelse(newdata$genestartnb>=duplication_tablerow1$startxnb & newdata$genestartnb<=duplication_tablerow1$endxnb,"Inter-chromosomal ohnolog",ifelse(newdata$geneendnb>=duplication_tablerow1$startxnb & newdata$geneendnb<=duplication_tablerow1$endxnb,"Inter-chromosomal ohnolog",ifelse(newdata$genestartnb>=chromadd[newdata$snpchrom] & newdata$genestartnb<=chromadd[newdata$snpchrom+1],ifelse(newdata$genestartnb>=duplication_tablerow1$startynb & newdata$genestartnb<=duplication_tablerow1$endynb,"Intra-chromosomal syntenic",ifelse(newdata$geneendnb>=duplication_tablerow1$startynb & newdata$geneendnb<=duplication_tablerow1$endynb,"Intra-chromosomal syntenic","Intra-chromosomal non syntenic")),"eQTL present"))),ifelse(newdata$genestartnb>=duplication_tablerow1$startynb & newdata$genestartnb<=duplication_tablerow1$endynb,"Gene present",ifelse(newdata$geneendnb>=duplication_tablerow1$startynb & newdata$geneendnb<=duplication_tablerow1$endynb,"Gene present",ifelse(newdata$genestartnb>=duplication_tablerow1$startxnb & newdata$genestartnb<=duplication_tablerow1$endxnb,"Gene present",ifelse(newdata$geneendnb>=duplication_tablerow1$startxnb & newdata$geneendnb<=duplication_tablerow1$endxnb,"Gene present","NA"))))))
#test2
#test<-ifelse(newdata$snpposnb>=duplication_tablerow1$startxnb & newdata$snpposnb<=duplication_tablerow1$endxnb,ifelse(newdata$genestartnb>=duplication_tablerow1$startynb & newdata$genestartnb<=duplication_tablerow1$endynb,"Inter-chromosomal ohnolog",ifelse(newdata$geneendnb>=duplication_tablerow1$startynb & newdata$geneendnb<=duplication_tablerow1$endynb,"Inter-chromosomal ohnolog",ifelse(newdata$genestartnb>=chromadd[newdata$snpchrom] & newdata$genestartnb<=chromadd[newdata$snpchrom+1],ifelse(newdata$genestartnb>=duplication_tablerow1$startxnb & newdata$genestartnb<=duplication_tablerow1$endxnb,"Intra-chromosomal syntenic","Intra-chromosomal non syntenic"),"NA"))),ifelse(newdata$snpposnb>=duplication_tablerow1$startynb & newdata$snpposnb<=duplication_tablerow1$endynb,ifelse(newdata$genestartnb>=duplication_tablerow1$startxnb & newdata$genestartnb<=duplication_tablerow1$endxnb,"Inter-chromosomal ohnolog",ifelse(newdata$geneendnb>=duplication_tablerow1$startxnb & newdata$geneendnb<=duplication_tablerow1$endxnb,"Inter-chromosomal ohnolog",ifelse(newdata$genestartnb>=chromadd[newdata$snpchrom] & newdata$genestartnb<=chromadd[newdata$snpchrom+1],ifelse(newdata$genestartnb>=duplication_tablerow1$startynb & newdata$genestartnb<=duplication_tablerow1$endynb,"Intra-chromosomal syntenic","Intra-chromosomal non syntenic"),"NA"))),"NA"))

#test3<-ifelse(newdata$genestartnb>=duplication_tablerow1$startynb & newdata$genestartnb<=duplication_tablerow1$endynb,"Inter-chromosomal ohnolog",ifelse(newdata$geneendnb>=duplication_tablerow1$startynb & newdata$geneendnb<=duplication_tablerow1$endynb,"Inter-chromosomal ohnolog",ifelse(newdata$genestartnb>=chromadd[newdata$snpchrom] & newdata$genestartnb<=chromadd[newdata$snpchrom+1],ifelse(newdata$genestartnb>=duplication_tablerow1$startxnb & newdata$genestartnb<=duplication_tablerow1$endxnb,"Intra-chromosomal syntenic","Intra-chromosomal non syntenic"),"NA")))

#ifelse(newdata$snpposnb>=duplication_tablerow1$startynb & newdata$snpposnb<=duplication_tablerow1$endynb,ifelse(newdata$genestartnb>=duplication_tablerow1$startxnb & newdata$genestartnb<=duplication_tablerow1$endxnb,"Inter-chromosomal ohnolog",ifelse(newdata$geneendnb>=duplication_tablerow1$startxnb & newdata$geneendnb<=duplication_tablerow1$endxnb,"Inter-chromosomal ohnolog",ifelse(newdata$genestartnb>=chromadd[newdata$snpchrom] & newdata$genestartnb<=chromadd[newdata$snpchrom+1],ifelse(newdata$genestartnb>=duplication_tablerow1$startynb & newdata$genestartnb<=duplication_tablerow1$endynb,"Intra-chromosomal syntenic","Intra-chromosomal non syntenic"),"NA"))),"NA")



#ifelse(newdata$genestartnb>=duplication_tablerow1$startynb & newdata$genestartnb<=duplication_tablerow1$endynb)

#ifelse(newdata$geneendnb>=duplication_tablerow1$startynb & newdata$geneendnb<=duplication_tablerow1$endynb,"Inter-chromosomal ohnolog",ifelse(newdata$genestartnb>=chromadd[newdata$snpchrom] & newdata$genestartnb<=chromadd[newdata$snpchrom+1],ifelse(newdata$genestartnb>=duplication_tablerow1$startxnb & newdata$genestartnb<=duplication_tablerow1$endxnb,"Intra-chromosomal syntenic","Intra-chromosomal non syntenic"),"NA"))

#ifelse(newdata$genestartnb>=chromadd[newdata$snpchrom] & newdata$genestartnb<=chromadd[newdata$snpchrom+1],ifelse(newdata$genestartnb>=duplication_tablerow1$startxnb & newdata$genestartnb<=duplication_tablerow1$endxnb,"Intra-chromosomal syntenic","Intra-chromosomal non syntenic"),"NA")

#ifelse(newdata$genestartnb>=duplication_tablerow1$startxnb & newdata$genestartnb<=duplication_tablerow1$endxnb,"Intra-chromosomal syntenic","Intra-chromosomal non syntenic")


newdata2<-newdata
lengthdata<-length(newdata2)
for(i in 1:length(duplication_table$id)){
  duplication_tablerow1<-duplication_table[i,]
  
  test<-ifelse(newdata$snpposnb>=duplication_tablerow1$startxnb & newdata$snpposnb<=duplication_tablerow1$endxnb,ifelse(newdata$genestartnb>=duplication_tablerow1$startynb & newdata$genestartnb<=duplication_tablerow1$endynb,"Inter-chromosomal ohnolog",ifelse(newdata$geneendnb>=duplication_tablerow1$startynb & newdata$geneendnb<=duplication_tablerow1$endynb,"Inter-chromosomal ohnolog",ifelse(newdata$genestartnb>=chromadd[newdata$snpchrom] & newdata$genestartnb<=chromadd[newdata$snpchrom+1],ifelse(newdata$genestartnb>=duplication_tablerow1$startxnb & newdata$genestartnb<=duplication_tablerow1$endxnb,"Intra-chromosomal syntenic",ifelse(newdata$geneendnb>=duplication_tablerow1$startxnb & newdata$geneendnb<=duplication_tablerow1$endxnb,"Intra-chromosomal syntenic","Intra-chromosomal non syntenic")),"eQTL present"))),ifelse(newdata$snpposnb>=duplication_tablerow1$startynb & newdata$snpposnb<=duplication_tablerow1$endynb,ifelse(newdata$genestartnb>=duplication_tablerow1$startxnb & newdata$genestartnb<=duplication_tablerow1$endxnb,"Inter-chromosomal ohnolog",ifelse(newdata$geneendnb>=duplication_tablerow1$startxnb & newdata$geneendnb<=duplication_tablerow1$endxnb,"Inter-chromosomal ohnolog",ifelse(newdata$genestartnb>=chromadd[newdata$snpchrom] & newdata$genestartnb<=chromadd[newdata$snpchrom+1],ifelse(newdata$genestartnb>=duplication_tablerow1$startynb & newdata$genestartnb<=duplication_tablerow1$endynb,"Intra-chromosomal syntenic",ifelse(newdata$geneendnb>=duplication_tablerow1$startynb & newdata$geneendnb<=duplication_tablerow1$endynb,"Intra-chromosomal syntenic","Intra-chromosomal non syntenic")),"eQTL present"))),ifelse(newdata$genestartnb>=duplication_tablerow1$startynb & newdata$genestartnb<=duplication_tablerow1$endynb,"Gene present",ifelse(newdata$geneendnb>=duplication_tablerow1$startynb & newdata$geneendnb<=duplication_tablerow1$endynb,"Gene present",ifelse(newdata$genestartnb>=duplication_tablerow1$startxnb & newdata$genestartnb<=duplication_tablerow1$endxnb,"Gene present",ifelse(newdata$geneendnb>=duplication_tablerow1$startxnb & newdata$geneendnb<=duplication_tablerow1$endxnb,"Gene present","NA"))))))
  newdata2[,lengthdata+i]<-test
}


## then process to get one value output

list_of_values <- c("Inter-chromosomal ohnolog", "Intra-chromosomal syntenic","Intra-chromosomal non syntenic")

# Define a function to apply row-wise
process_row <- function(row) {
  if (any(row %in% list_of_values)) {
    return(row[row %in% list_of_values][1])  # Output the first found value
  } else if ("Gene present" %in% row && "eQTL present" %in% row) {
    return("Inter-chromosomal non ohnolog")
  } else {
    return(NA)
  }
}

# Apply the function to each row using apply()
newdata2$result <- apply(newdata2, 1, process_row)

newdata$duplicatestatus<-newdata2[,103]

write.table(newdata,file="rna_LL_all_august_version_01.txt",row.names = FALSE)

table(newdata[newdata$kind=="trans","duplicatestatus"])
newdatatest<-newdata[!is.na(newdata$duplicatestatus),]
data<-newdatatest[newdatatest$kind=="cis" & newdatatest$duplicatestatus=="Inter-chromosomal ohnolog",]



## shuffle code

shuffle_data<- function(table){
  cols_to_shuffle <- c("chr","start","end")
  cols_to_shuffle2<-c("snpchrom","snppos")
  
  # combine the columns you want to shuffle into one matrix
  shuffle_matrix <- table[, cols_to_shuffle]
  shuffle_matrix2 <- table[, cols_to_shuffle2]
  
  
  # shuffle the rows of the matrix
  shuffle_matrix <- shuffle_matrix[sample(nrow(shuffle_matrix)),]
  shuffle_matrix2 <- shuffle_matrix2[sample(nrow(shuffle_matrix2)),]
  
  # split the shuffled matrix back into the original columns
  newtab<- cbind.data.frame(shuffle_matrix, shuffle_matrix2)
  return(newtab)
  
}

shuffleddata<-shuffle_data(newdata)

## convert position to number

lengthvector<-c(174498729,95481959,105780080,90536438,92788608,96060288,68862998,28860523,161282225,125877811,111868677,101677876,
                114417674,101980477,110670232,96486271,87489397,84084598,88107222,96847506,59819933,63823863,52460201,49354470,54385492,
                55994222,45305548,41468476,43051128)
chromadd<-c(0,cumsum(lengthvector)[1:28])

chromname<-c(1:29)
chrommid<-chromadd+lengthvector/2
listval<-cbind.data.frame(chromname,lengthvector,chromadd,chrommid)

shuffleddata$snpposnb<-shuffleddata$snppos+listval[match(shuffleddata$snpchrom,listval$chromname),3]

shuffleddata$genestartnb<-shuffleddata$start+listval[match(shuffleddata$chr,listval$chromname),3]

shuffleddata$geneendnb<-shuffleddata$end+listval[match(shuffleddata$chr,listval$chromname),3]



shuffleddata2<-shuffleddata
lengthdata<-length(shuffleddata2)
for(i in 1:length(duplication_table$id)){
  duplication_tablerow1<-duplication_table[i,]
  
  test<-ifelse(shuffleddata$snpposnb>=duplication_tablerow1$startxnb & shuffleddata$snpposnb<=duplication_tablerow1$endxnb,ifelse(shuffleddata$genestartnb>=duplication_tablerow1$startynb & shuffleddata$genestartnb<=duplication_tablerow1$endynb,"Inter-chromosomal ohnolog",ifelse(shuffleddata$geneendnb>=duplication_tablerow1$startynb & shuffleddata$geneendnb<=duplication_tablerow1$endynb,"Inter-chromosomal ohnolog",ifelse(shuffleddata$genestartnb>=chromadd[shuffleddata$snpchrom] & shuffleddata$genestartnb<=chromadd[shuffleddata$snpchrom+1],ifelse(shuffleddata$genestartnb>=duplication_tablerow1$startxnb & shuffleddata$genestartnb<=duplication_tablerow1$endxnb,"Intra-chromosomal syntenic",ifelse(shuffleddata$geneendnb>=duplication_tablerow1$startxnb & shuffleddata$geneendnb<=duplication_tablerow1$endxnb,"Intra-chromosomal syntenic","Intra-chromosomal non syntenic")),"eQTL present"))),ifelse(shuffleddata$snpposnb>=duplication_tablerow1$startynb & shuffleddata$snpposnb<=duplication_tablerow1$endynb,ifelse(shuffleddata$genestartnb>=duplication_tablerow1$startxnb & shuffleddata$genestartnb<=duplication_tablerow1$endxnb,"Inter-chromosomal ohnolog",ifelse(shuffleddata$geneendnb>=duplication_tablerow1$startxnb & shuffleddata$geneendnb<=duplication_tablerow1$endxnb,"Inter-chromosomal ohnolog",ifelse(shuffleddata$genestartnb>=chromadd[shuffleddata$snpchrom] & shuffleddata$genestartnb<=chromadd[shuffleddata$snpchrom+1],ifelse(shuffleddata$genestartnb>=duplication_tablerow1$startynb & shuffleddata$genestartnb<=duplication_tablerow1$endynb,"Intra-chromosomal syntenic",ifelse(shuffleddata$geneendnb>=duplication_tablerow1$startynb & shuffleddata$geneendnb<=duplication_tablerow1$endynb,"Intra-chromosomal syntenic","Intra-chromosomal non syntenic")),"eQTL present"))),ifelse(shuffleddata$genestartnb>=duplication_tablerow1$startynb & shuffleddata$genestartnb<=duplication_tablerow1$endynb,"Gene present",ifelse(shuffleddata$geneendnb>=duplication_tablerow1$startynb & shuffleddata$geneendnb<=duplication_tablerow1$endynb,"Gene present",ifelse(shuffleddata$genestartnb>=duplication_tablerow1$startxnb & shuffleddata$genestartnb<=duplication_tablerow1$endxnb,"Gene present",ifelse(shuffleddata$geneendnb>=duplication_tablerow1$startxnb & shuffleddata$geneendnb<=duplication_tablerow1$endxnb,"Gene present","NA"))))))
  shuffleddata2[,lengthdata+i]<-test
}


## then process to get one value output

list_of_values <- c("Inter-chromosomal ohnolog", "Intra-chromosomal syntenic","Intra-chromosomal non syntenic")

# Define a function to apply row-wise
process_row <- function(row) {
  if (any(row %in% list_of_values)) {
    return(row[row %in% list_of_values][1])  # Output the first found value
  } else if ("Gene present" %in% row && "eQTL present" %in% row) {
    return("Inter-chromosomal non ohnolog")
  } else {
    return(NA)
  }
}

# Apply the function to each row using apply()
shuffleddata2$result <- apply(shuffleddata2, 1, process_row)

shuffleddata$duplicatestatus<-shuffleddata2[,96]

###plots and prop


newdata$data<-"real"
shuffleddata$data<-"shuffled"
shuffleddata<-cbind.data.frame(shuffleddata,newdata[,c(1:6,10)])

fulltab<-rbind.data.frame(newdata,shuffleddata,stringsAsFactors = FALSE)

library(dplyr)
library(ggplot2)
library(tidyr)

dupstatus_list <- c("Inter-chromosomal ohnolog", "Intra-chromosomal syntenic","Intra-chromosomal non syntenic","Inter-chromosomal non ohnolog")

dfcount <- fulltab%>% filter(!is.na(duplicatestatus))  %>% filter(kind!="cis")  %>% group_by(duplicatestatus,data) %>% 
  summarise(total_count=n(),.groups = 'drop') %>%
  as.data.frame() %>% filter(duplicatestatus %in% dupstatus_list) %>% group_by(data) %>% 
  mutate(prop = total_count / sum(total_count))

#dfcount$dupstatus<-as.factor(dfcount$dupstatus)
#dfcount$dupstatus<-factor(dfcount$dupstatus,levels=c("inter_chr_non-ohnolog","ohnolog","intra_chr_non-syntenic","intra_chr_syntenic"))




ggplot(dfcount,aes(x=data,y=prop,fill=duplicatestatus))+geom_bar(position = "fill", stat = "identity")+scale_fill_manual(values=c("#F3B5B3","#F5E69A","#B5AFDE","#B6D8B6"))+theme_bw(base_size = 30)+ylab("Proportion") +xlab("Data type")+guides(fill=guide_legend("Interaction type"))

ggplot(dfcount,aes(x="",y=prop,fill=duplicatestatus))+geom_bar(stat="identity",width=1)+
  coord_polar("y",start=0)+theme_void()+scale_fill_manual(values=c("#F3B5B3","#F5E69A","#B5AFDE","#B6D8B6"))+
  theme(legend.text = element_text(size=15))


df_wide <- pivot_wider(dfcount, id_cols = duplicatestatus, names_from = data, values_from = total_count)
df_wide$enrichment<-df_wide$real/df_wide$shuffled




### bootstrap
real_data<-read.table(file="rna_LL_all_august_version_01.txt",header=TRUE)
real_data$data<-"real"

shuffle_data<- function(table){
  cols_to_shuffle <- c("chr","start","end")
  cols_to_shuffle2<-c("snpchrom","snppos")
  
  # combine the columns you want to shuffle into one matrix
  shuffle_matrix <- table[, cols_to_shuffle]
  shuffle_matrix2 <- table[, cols_to_shuffle2]
  
  
  # shuffle the rows of the matrix
  shuffle_matrix <- shuffle_matrix[sample(nrow(shuffle_matrix)),]
  shuffle_matrix2 <- shuffle_matrix2[sample(nrow(shuffle_matrix2)),]
  
  # split the shuffled matrix back into the original columns
  newtab<- cbind.data.frame(shuffle_matrix, shuffle_matrix2)
  return(newtab)
  
}

Bootstrap_trans<-function(nb,table_init){
  
  library(dplyr)
  library(tidyr)
  
  duplication_table<-read.table(file="https://salmobase.org/datafiles/TSV/synteny/2021-11/AtlanticSalmon/synteny.tsv",header=TRUE,sep="\t")
  
  
  
  lengthvector<-c(174498729,95481959,105780080,90536438,92788608,96060288,68862998,28860523,161282225,125877811,111868677,101677876,
                  114417674,101980477,110670232,96486271,87489397,84084598,88107222,96847506,59819933,63823863,52460201,49354470,54385492,
                  55994222,45305548,41468476,43051128)
  chromadd<-c(0,cumsum(lengthvector)[1:28])
  
  chromname<-c(1:29)
  chrommid<-chromadd+lengthvector/2
  listval<-cbind.data.frame(chromname,lengthvector,chromadd,chrommid)
  
  # Define a function to apply row-wise
  process_row <- function(row) {
    if (any(row %in% list_of_values)) {
      return(row[row %in% list_of_values][1])  # Output the first found value
    } else if ("Gene present" %in% row && "eQTL present" %in% row) {
      return("Inter-chromosomal non ohnolog")
    } else {
      return(NA)
    }
  }
  
  dupstatus_list <- c("Inter-chromosomal ohnolog", "Intra-chromosomal syntenic","Intra-chromosomal non syntenic","Inter-chromosomal non ohnolog")
  df_bootstrap<-c()
  for (y in 1:nb){
    
    shuffleddata<-shuffle_data(table_init)
    
    ## convert position to number
    
    duplication_table$startxnb<-duplication_table$begin_x+listval[match(duplication_table$chromosome_x,listval$chromname),3]
    
    duplication_table$endxnb<-duplication_table$end_x+listval[match(duplication_table$chromosome_x,listval$chromname),3]
    
    duplication_table$startynb<-duplication_table$begin_y+listval[match(duplication_table$chromosome_y,listval$chromname),3]
    
    duplication_table$endynb<-duplication_table$end_y+listval[match(duplication_table$chromosome_y,listval$chromname),3]
    
    shuffleddata$snpposnb<-shuffleddata$snppos+listval[match(shuffleddata$snpchrom,listval$chromname),3]
    
    shuffleddata$genestartnb<-shuffleddata$start+listval[match(shuffleddata$chr,listval$chromname),3]
    
    shuffleddata$geneendnb<-shuffleddata$end+listval[match(shuffleddata$chr,listval$chromname),3]
    
    
    
    shuffleddata2<-shuffleddata
    lengthdata<-length(shuffleddata2)
    for(i in 1:length(duplication_table$id)){
      duplication_tablerow1<-duplication_table[i,]
      
      test<-ifelse(shuffleddata$snpposnb>=duplication_tablerow1$startxnb & shuffleddata$snpposnb<=duplication_tablerow1$endxnb,ifelse(shuffleddata$genestartnb>=duplication_tablerow1$startynb & shuffleddata$genestartnb<=duplication_tablerow1$endynb,"Inter-chromosomal ohnolog",ifelse(shuffleddata$geneendnb>=duplication_tablerow1$startynb & shuffleddata$geneendnb<=duplication_tablerow1$endynb,"Inter-chromosomal ohnolog",ifelse(shuffleddata$genestartnb>=chromadd[shuffleddata$snpchrom] & shuffleddata$genestartnb<=chromadd[shuffleddata$snpchrom+1],ifelse(shuffleddata$genestartnb>=duplication_tablerow1$startxnb & shuffleddata$genestartnb<=duplication_tablerow1$endxnb,"Intra-chromosomal syntenic",ifelse(shuffleddata$geneendnb>=duplication_tablerow1$startxnb & shuffleddata$geneendnb<=duplication_tablerow1$endxnb,"Intra-chromosomal syntenic","Intra-chromosomal non syntenic")),"eQTL present"))),ifelse(shuffleddata$snpposnb>=duplication_tablerow1$startynb & shuffleddata$snpposnb<=duplication_tablerow1$endynb,ifelse(shuffleddata$genestartnb>=duplication_tablerow1$startxnb & shuffleddata$genestartnb<=duplication_tablerow1$endxnb,"Inter-chromosomal ohnolog",ifelse(shuffleddata$geneendnb>=duplication_tablerow1$startxnb & shuffleddata$geneendnb<=duplication_tablerow1$endxnb,"Inter-chromosomal ohnolog",ifelse(shuffleddata$genestartnb>=chromadd[shuffleddata$snpchrom] & shuffleddata$genestartnb<=chromadd[shuffleddata$snpchrom+1],ifelse(shuffleddata$genestartnb>=duplication_tablerow1$startynb & shuffleddata$genestartnb<=duplication_tablerow1$endynb,"Intra-chromosomal syntenic",ifelse(shuffleddata$geneendnb>=duplication_tablerow1$startynb & shuffleddata$geneendnb<=duplication_tablerow1$endynb,"Intra-chromosomal syntenic","Intra-chromosomal non syntenic")),"eQTL present"))),ifelse(shuffleddata$genestartnb>=duplication_tablerow1$startynb & shuffleddata$genestartnb<=duplication_tablerow1$endynb,"Gene present",ifelse(shuffleddata$geneendnb>=duplication_tablerow1$startynb & shuffleddata$geneendnb<=duplication_tablerow1$endynb,"Gene present",ifelse(shuffleddata$genestartnb>=duplication_tablerow1$startxnb & shuffleddata$genestartnb<=duplication_tablerow1$endxnb,"Gene present",ifelse(shuffleddata$geneendnb>=duplication_tablerow1$startxnb & shuffleddata$geneendnb<=duplication_tablerow1$endxnb,"Gene present","NA"))))))
      shuffleddata2[,lengthdata+i]<-test
    }
    
    
    ## then process to get one value output
    
    list_of_values <- c("Inter-chromosomal ohnolog", "Intra-chromosomal syntenic","Intra-chromosomal non syntenic")
    
    # Apply the function to each row using apply()
    shuffleddata2$result <- apply(shuffleddata2, 1, process_row)
    
    shuffleddata$duplicatestatus<-shuffleddata2[,96]
    
    shuffleddata$data<-"shuffled"
    shuffleddata<-cbind.data.frame(shuffleddata,table_init[,c(1:6,10)])
    
    fulltab<-rbind.data.frame(table_init,shuffleddata,stringsAsFactors = FALSE)
    
    
    dfcount <- fulltab%>% filter(!is.na(duplicatestatus))  %>% filter(kind!="cis")  %>% group_by(duplicatestatus,data) %>% 
      summarise(total_count=n(),.groups = 'drop') %>%
      as.data.frame() %>% filter(duplicatestatus %in% dupstatus_list) %>% group_by(data) %>% 
      mutate(prop = total_count / sum(total_count))
    
    df_wide <- pivot_wider(dfcount, id_cols = duplicatestatus, names_from = data, values_from = total_count)
    df_wide$enrichment<-df_wide$real/df_wide$shuffled
    df_wide$bootstrapnb<-y
    
    df_bootstrap<-rbind.data.frame(df_bootstrap,df_wide,stringsAsFactors = FALSE)
    
  }
  return(df_bootstrap)
}

test<-Bootstrap_trans(nb=1000,table_init = real_data)


write.table(test,file="bootstrap_rnaLL_1000_01.txt",row.names = FALSE)


library(dplyr)

# Compute mean and standard error
result <- test %>%
  group_by(duplicatestatus) %>%
  summarize(
    mean_enrichment = mean(enrichment),
    se_enrichment = sd(enrichment) / sqrt(n())
  )

library(ggplot2)

ggplot(result,aes(x=duplicatestatus,y=mean_enrichment,fill=duplicatestatus))+geom_bar(stat="identity")+scale_fill_manual(values=c("#F3B5B3","#F5E69A","#B5AFDE","#B6D8B6"))+theme_bw(30)+geom_hline(yintercept=1,color="darkgrey",linetype = "dashed",size=4)+ scale_y_continuous(trans = "log2")+
  xlab("Interaction type") +ylab("Log2 of enrichment")+geom_errorbar(aes(x=duplicatestatus, ymin = mean_enrichment - 1.96*se_enrichment, ymax = mean_enrichment + 1.96*se_enrichment), color = "black")#+theme(axis.text.x = element_text(angle = 45, hjust=1)) 
