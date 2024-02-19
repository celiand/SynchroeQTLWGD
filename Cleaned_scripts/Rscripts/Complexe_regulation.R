### Celian Diblasi
### 19/02/2024
### Investigation of regulation of complexes

##### First make data for co regulation (this part is quite heavy so run it on cluster if possible) ####


real_data<-read.table(file="/mnt/users/cedi/rna_LL_all_august_version_01.txt",header=TRUE) ## data from XXX

## first import data
file_path<-"/mnt/users/cedi/8030.protein.links.v12.0.txt.gz" ##data from string database
data <- read.table(file = gzfile(file_path), header = TRUE)

# Close the connection to the gzfile
close(gzfile(file_path))

### Format data

data$protein1<-sapply(strsplit(data$protein1,".",fixed=TRUE),FUN = `[[`, 2)
data$protein2<-sapply(strsplit(data$protein2,".",fixed=TRUE),FUN = `[[`, 2)

###  Get genes that code for protein info - can be found in personalized info on ensembl ###
protgene<-read.table(file="/mnt/users/cedi/mart_export.txt",header=TRUE,sep="\t")

protgene<-protgene[,c(1,5)]

eqtl_prot<-merge(real_data,protgene,by.x="gene",by.y="Gene.stable.ID")


##in case they are several regulator, just split cis and trans so we will have maximum 1 on each
library(tidyverse)

data_pluseqtlp1<-data%>%
  mutate(cis_eqtl = ifelse((protein1 %in% eqtl_prot$Protein.stable.ID) & eqtl_prot$kind=="cis" ,eqtl_prot$lead_SNP , "NA"),
         trans_eqtl = ifelse((protein1 %in% eqtl_prot$Protein.stable.ID) & eqtl_prot$kind=="trans" ,eqtl_prot$lead_SNP , "NA"))

##first create data frames with match cis and trans eQTl for protein 1 and 2
eqtlprot_plusp1_cis<-eqtl_prot %>% filter(kind=="cis") %>% mutate(cis_eQTL_P1=ifelse(Protein.stable.ID %in% data$protein1,lead_SNP,NA))
eqtlprot_plusp1_trans<-eqtl_prot %>% filter(kind=="trans") %>% mutate(trans_eQTL_P1=ifelse(Protein.stable.ID %in% data$protein1,lead_SNP,NA))

eqtlprot_plusp2_cis<-eqtl_prot %>% filter(kind=="cis") %>% mutate(cis_eQTL=ifelse(Protein.stable.ID %in% data$protein2,lead_SNP,NA))
eqtlprot_plusp2_trans<-eqtl_prot %>% filter(kind=="trans") %>% mutate(trans_eQTL=ifelse(Protein.stable.ID %in% data$protein2,lead_SNP,NA))

##then make the matches between the column

##first with P1
data_P1_ciseqtl<-rep(NA, length(data$protein1))
matches <- match(data$protein1, eqtlprot_plusp1_cis$Protein.stable.ID)

data_P1_ciseqtl[!is.na(matches)] <- eqtlprot_plusp1_cis$cis_eQTL[matches[!is.na(matches)]]
data$cis_eQTL_P1<-data_P1_ciseqtl


data_P1_transeqtl<-rep(NA, length(data$protein1))
matches <- match(data$protein1, eqtlprot_plusp1_trans$Protein.stable.ID)

data_P1_transeqtl[!is.na(matches)] <- eqtlprot_plusp1_trans$trans_eQTL[matches[!is.na(matches)]]
data$trans_eQTL_P1<-data_P1_transeqtl


##then with P2
data_P2_ciseqtl<-rep(NA, length(data$protein2))
matches <- match(data$protein2, eqtlprot_plusp2_cis$Protein.stable.ID)

data_P2_ciseqtl[!is.na(matches)] <- eqtlprot_plusp2_cis$cis_eQTL[matches[!is.na(matches)]]
data$cis_eQTL_P2<-data_P2_ciseqtl


data_P2_transeqtl<-rep(NA, length(data$protein2))
matches <- match(data$protein2, eqtlprot_plusp2_trans$Protein.stable.ID)

data_P2_transeqtl[!is.na(matches)] <- eqtlprot_plusp2_trans$trans_eQTL[matches[!is.na(matches)]]
data$trans_eQTL_P2<-data_P2_transeqtl


data$co_reg<-ifelse((data$cis_eQTL_P1==data$cis_eQTL_P2) | (data$cis_eQTL_P1==data$trans_eQTL_P2), data$cis_eQTL_P1,ifelse((data$trans_eQTL_P1==data$cis_eQTL_P2) | (data$trans_eQTL_P1==data$trans_eQTL_P2),data$trans_eQTL_P1,NA))

##write out results
data2<-data[!is.na(data$co_reg),]
write.table(data2,"/mnt/users/cedi/coreg_data_filtered.txt",row.names=FALSE)

write.table(data,"/mnt/users/cedi/coreg_data_full.txt",row.names=FALSE)


##### Then make random expectations (again this part is quite heavy) #####


#### make a data with no double (keep only one of P1 - P2 and P2 - P1)

file_path<-"/mnt/users/cedi/8030.protein.links.v12.0.txt.gz"
data <- read.table(file = gzfile(file_path), header = TRUE)


close(gzfile(file_path))

data$protein1<-sapply(strsplit(data$protein1,".",fixed=TRUE),FUN = `[[`, 2)
data$protein2<-sapply(strsplit(data$protein2,".",fixed=TRUE),FUN = `[[`, 2)


### clean double

datafilter <- data %>%
mutate(
protein1_new = pmin(protein1, protein2),
protein2_new = pmax(protein1, protein2)
) %>%
distinct(protein1_new, protein2_new, .keep_all = TRUE)

datafilter<-datafilter[,1:2]

write.table(datafilter,file="/mnt/users/cedi/coregdata_nodouble.txt",row.names=FALSE)

### then make the random expectation
library(tidyverse)

real_data<-read.table(file="/mnt/users/cedi/rna_LL_all_august_version_01.txt",header=TRUE)


datafilter<-read.table(file="/mnt/users/cedi/coregdata_nodouble.txt",header=TRUE) ### get data without double

## make a vector with all protein present (keeping number of occurences they appear)
allprotein<-c(datafilter[,1],datafilter[,2])

nbprot<-length(allprotein)

result_shuffle<-c()
for(i in 1:100){
  shuffled_vector <- sample(allprotein)
  
  
  
  association_table <- data.frame(protein1=shuffled_vector[1:(nbprot/2)],protein2=shuffled_vector[((nbprot/2)+1):nbprot])
  
  association_table_bis<-association_table[!(association_table$protein1 == association_table$protein2),]
  
  
  data<-association_table_bis
  
  
  protgene<-read.table(file="/mnt/users/cedi/mart_export.txt",header=TRUE,sep="\t")
  
  protgene<-protgene[,c(1,5)]
  
  eqtl_prot<-merge(real_data,protgene,by.x="gene",by.y="Gene.stable.ID")
  
  
  ##in case they are several regulator, just split cis and trans so we will have maximum 1 on each
  library(tidyverse)
  
  data_pluseqtlp1<-data%>%
    mutate(cis_eqtl = ifelse((protein1 %in% eqtl_prot$Protein.stable.ID) & eqtl_prot$kind=="cis" ,eqtl_prot$lead_SNP , "NA"),
           trans_eqtl = ifelse((protein1 %in% eqtl_prot$Protein.stable.ID) & eqtl_prot$kind=="trans" ,eqtl_prot$lead_SNP , "NA"))
  
  ##first create data frames with match cis and trans eQTl for protein 1 and 2
  eqtlprot_plusp1_cis<-eqtl_prot %>% filter(kind=="cis") %>% mutate(cis_eQTL_P1=ifelse(Protein.stable.ID %in% data$protein1,lead_SNP,NA))
  eqtlprot_plusp1_trans<-eqtl_prot %>% filter(kind=="trans") %>% mutate(trans_eQTL_P1=ifelse(Protein.stable.ID %in% data$protein1,lead_SNP,NA))
  
  eqtlprot_plusp2_cis<-eqtl_prot %>% filter(kind=="cis") %>% mutate(cis_eQTL=ifelse(Protein.stable.ID %in% data$protein2,lead_SNP,NA))
  eqtlprot_plusp2_trans<-eqtl_prot %>% filter(kind=="trans") %>% mutate(trans_eQTL=ifelse(Protein.stable.ID %in% data$protein2,lead_SNP,NA))
  
  ##then make the matches between the column
  
  ##first with P1
  data_P1_ciseqtl<-rep(NA, length(data$protein1))
  matches <- match(data$protein1, eqtlprot_plusp1_cis$Protein.stable.ID)
  
  data_P1_ciseqtl[!is.na(matches)] <- eqtlprot_plusp1_cis$cis_eQTL[matches[!is.na(matches)]]
  data$cis_eQTL_P1<-data_P1_ciseqtl
  
  
  data_P1_transeqtl<-rep(NA, length(data$protein1))
  matches <- match(data$protein1, eqtlprot_plusp1_trans$Protein.stable.ID)
  
  data_P1_transeqtl[!is.na(matches)] <- eqtlprot_plusp1_trans$trans_eQTL[matches[!is.na(matches)]]
  data$trans_eQTL_P1<-data_P1_transeqtl
  
  
  ##then with P2
  data_P2_ciseqtl<-rep(NA, length(data$protein2))
  matches <- match(data$protein2, eqtlprot_plusp2_cis$Protein.stable.ID)
  
  data_P2_ciseqtl[!is.na(matches)] <- eqtlprot_plusp2_cis$cis_eQTL[matches[!is.na(matches)]]
  data$cis_eQTL_P2<-data_P2_ciseqtl
  
  
  data_P2_transeqtl<-rep(NA, length(data$protein2))
  matches <- match(data$protein2, eqtlprot_plusp2_trans$Protein.stable.ID)
  
  data_P2_transeqtl[!is.na(matches)] <- eqtlprot_plusp2_trans$trans_eQTL[matches[!is.na(matches)]]
  data$trans_eQTL_P2<-data_P2_transeqtl
  
  
  data$co_reg<-ifelse((data$cis_eQTL_P1==data$cis_eQTL_P2) | (data$cis_eQTL_P1==data$trans_eQTL_P2), data$cis_eQTL_P1,ifelse((data$trans_eQTL_P1==data$cis_eQTL_P2) | (data$trans_eQTL_P1==data$trans_eQTL_P2),data$trans_eQTL_P1,NA))
  
  ## analyze and save result
  
  ##keep only coreg cases
  data2<-data[!is.na(data$co_reg),]
  
  # Create a table using the table() function
  my_table <- table(data2$co_reg)
  
  # Convert the table to a data frame
  df <- as.data.frame(my_table)
  
  # Sort the data frame in descending order by the frequency column
  df_sorted <- df[order(-df$Freq), ]
  
  df_sorted$run<-i
  
  result_shuffle<-rbind.data.frame(result_shuffle,df_sorted,stringsAsFactors = FALSE)
  
}

write.table(result_shuffle,file="/mnt/users/cedi/result_shuffle_prot_V3_100it.txt",row.names=FALSE)

#### Then make analysis and figures

randomcoreg<-read.table(file="result_shuffle_prot_V3_100it.txt",header=TRUE)


datareg<-read.table(file="coreg_data_filtered.txt",header=TRUE)


##remove duplicate in coreg data (rows where prot 1 and prot 2 are swaped)
library(tidyverse)

dataregfilter <- datareg %>%
  mutate(
    protein1_new = pmin(protein1, protein2),
    protein2_new = pmax(protein1, protein2)
  ) %>%
  distinct(protein1_new, protein2_new, .keep_all = TRUE)

#table(dataregfilter$chr)
#dataregfilter[dataregfilter$chr=="ssa29",]

# Create a table using the table() function
my_table <- table(dataregfilter$co_reg)

# Convert the table to a data frame
df <- as.data.frame(my_table)

# Sort the data frame in descending order by the frequency column
df_sorted <- df[order(-df$Freq), ]

df_sorted$Var1<-as.character(df_sorted$Var1)

## compute mean, sd and run were is appears for all lead_SNP present and all other SNP
for ( i in unique(df_sorted$Var1)){
  subtabeqtl<-randomcoreg[randomcoreg$Var1==i,]
  nbrun<-length(subtabeqtl$run)
  meaneqtl<-mean(subtabeqtl$Freq)
  sdeqtl<-sd(subtabeqtl$Freq)
  df_sorted[df_sorted$Var1==i,"mean"]<-meaneqtl
  df_sorted[df_sorted$Var1==i,"sd"]<-sdeqtl
  df_sorted[df_sorted$Var1==i,"nbrun"]<-nbrun
}

nnpresenteqtl<-unique(randomcoreg[!(randomcoreg$Var1%in%df_sorted$Var1),"Var1"])

for ( i in nnpresenteqtl){
  subtabeqtl<-randomcoreg[randomcoreg$Var1==i,]
  nbrun<-length(subtabeqtl$run)
  meaneqtl<-mean(subtabeqtl$Freq)
  sdeqtl<-sd(subtabeqtl$Freq)
  row<-c(as.character(i),0,meaneqtl,sdeqtl,nbrun)
  df_sorted<-rbind.data.frame(df_sorted,row,stringsAsFactors = FALSE)
}

## format data
df_sorted2<-df_sorted[df_sorted$mean!="NaN",]

df_sorted2$mean<-as.numeric(as.character(df_sorted2$mean))

df_sorted2$Freq<-as.numeric(as.character(df_sorted2$Freq))

df_sorted2$sd<-as.numeric(as.character(df_sorted2$sd))

df_sorted2$nbrun<-as.numeric(as.character(df_sorted2$nbrun))

df_sorted2$ratio<-df_sorted2$Freq/df_sorted2$mean

df_sorted2$presence<-ifelse(df_sorted2$Freq==0,"no","yes")

df_sorted_3<-df_sorted2[df_sorted2$Freq>=1,]
## 154 eQTLs regulating 1 protein complex at least
res<-df_sorted_3[df_sorted_3$Freq>df_sorted_3$mean,]
## 44 eQTLs where observed > expected

library(ggplot2)


#ggplot(df_sorted2,aes(x=nbrun,y=ratio,color=presence))+geom_point()

#ggplot(df_sorted2,aes(x=Var1,y=ratio,color=presence))+geom_point()

df_sorted2$label<-""
df_sorted2[c(2,3,4,7),"label"]<-df_sorted2[c(2,3,4,7),"Var1"]

#ggplot(df_sorted2,aes(x=Freq,y=mean,color=presence))+geom_point(size=3)+geom_errorbar(aes(x=Freq, ymin = mean - sd, ymax = mean+ sd), color = "black")+geom_abline(intercept = 0, slope = 1, color = "red")+geom_text(aes(label = label), hjust=-0.3,color="black")



## perform statistical test
#pvalvec_wilc<-c()
pval_ttest<-c()

for ( i in unique(df_sorted2$Var1)){
  subtabeqtl<-randomcoreg[randomcoreg$Var1==i,"Freq"]
  freqobs<-df_sorted2[df_sorted2$Var1==i,"Freq"]
  #wilctest<-wilcox.test(freqobs,subtabeqtl)
  #wilcpval<-wilctest$p.value
  #pvalvec_wilc<-c(pvalvec_wilc,wilcpval)
  if(df_sorted2[df_sorted2$Var1==i,"nbrun"]>1 & df_sorted2[df_sorted2$Var1==i,"sd"]!=0){
    t_test_result <- t.test(subtabeqtl, mu = freqobs)
    t.test_pval<-t_test_result$p.value
    pval_ttest<-c(pval_ttest, t.test_pval)
  }else{
    pval_ttest<-c(pval_ttest, NA)
  }
  
}



adjusted_p_values <- p.adjust(pval_ttest, method = "bonferroni")

df_sorted2$pval_ttest<-pval_ttest
df_sorted2$adjPval_ttest<-adjusted_p_values


df_sorted2$significant_tt<-ifelse(df_sorted2$adjPval_ttest<=0.05,"yes","no")

df_sorted2_nona<-df_sorted2[!is.na(df_sorted2$significant_tt),]


df_sorted2_nona$label<-""
df_sorted2_nona[c(2,3,4,7),"label"]<-df_sorted2_nona[c(2,3,4,7),"Var1"]

ggplot(df_sorted2_nona,aes(x=Freq,y=mean,color=significant_tt))+geom_point(size=3)+geom_errorbar(aes(x=Freq, ymin = mean - sd, ymax = mean+ sd), color = "black")+geom_abline(intercept = 0, slope = 1, color = "#ee6b44")+theme_bw(25)+ylab("Expected number of complexes regulated")+xlab("Observed number of complexes regulated")+geom_text(aes(label = label), hjust=-0.2,color="black",size=5)+
  scale_color_manual(values=c("#bd31de","#22ba6e"))



## investigate which gene are the 4 outliers variants regulating (put obtained genes in GO enrichment)
real_data<-read.table(file="rna_LL_all_august_version_01.txt",header=TRUE)

real_datagene_n1<-real_data[real_data$lead_SNP=="ssa14_50870510","gene"]
write.table(real_datagene_n1,file="genesGO_ssa14_50870510.txt",row.names=FALSE)

real_datagene_n2<-real_data[real_data$lead_SNP=="ssa14_51370830","gene"]
write.table(real_datagene_n2,file="genesGO_ssa14_51370830.txt",row.names=FALSE)

real_datagene_n3<-real_data[real_data$lead_SNP=="ssa15_26623522","gene"]
write.table(real_datagene_n3,file="genesGO_sssa15_26623522.txt",row.names=FALSE)

real_datagene_n4<-real_data[real_data$lead_SNP=="ssa14_40168953","gene"]
write.table(real_datagene_n4,file="genesGO_ssa14_40168953.txt",row.names=FALSE)