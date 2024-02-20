#### Celian Diblasi
#### 20/02/2024
#### Investigate ohnlogs gene expression



## -- function to find pair -- ##

data<-read.table(file="synchroLL_finemap.txt",header=TRUE)

cisdata<-read.table(file="synchroLL_cis_finemap.txt",header=TRUE)

all_data<-read.table(file="rna_LL_all_august_version_01.txt",header=TRUE) ## data from Region_shuffle_bootstrap.R

all_cis<-all_data[all_data$kind=="cis",c(1,2,5)]
all_cis$posterior_prob<-NA

colnames(all_cis)[2]<-"SNP"
data$kind<-"trans"

cistransdata<-rbind.data.frame(all_cis,data,stringsAsFactors = FALSE)

cistransdata<-rbind.data.frame(cisdata,data,stringsAsFactors = FALSE)

genetab<-read.table(file="ss4r_dups_and_singletons_ENSrapid_convPipeline.tsv",header=TRUE) ## data of ohnologs
genetab<-genetab[genetab$type=="ss4r",]
genetab<-genetab[,-c(1,2,6)]


getPairstatus<-function(threshold){
  #data_subset<- cistransdata[is.na(cistransdata$posterior_prob) | cistransdata$posterior_prob >= threshold, ]
  data_subset<- cistransdata[cistransdata$posterior_prob >= threshold, ]
  for(y in 1:length(genetab$ID)){
    
    gene1row1<-genetab[y,"gene1"]
    tabpair1<-data_subset[data_subset$gene==gene1row1,]
    gene2row1<-genetab[y,"gene2"]
    tabpair2<-data_subset[data_subset$gene==gene2row1,]  
    
    ltab1<-length(tabpair1$SNP)
    ltab2<-length(tabpair2$SNP)
    
    reg<-"No eQTL"
    bothdiffreg<-"nobothdiffreg"
    bothsamereg<-"nobothsamereg"
    if(ltab1==0 && ltab2==0){
    }
    if(ltab1==0 && ltab2!=0){
      reg<-"Single copy eQTL"
    }
    if(ltab1!=0 && ltab2==0){
      reg<-"Single copy eQTL"
    }
    if(ltab1!=0 && ltab2!=0){
      issamesnp<-ifelse(tabpair1$SNP%in%tabpair2$SNP,"yes","no")
      if("no"%in%issamesnp){
        bothdiffreg<-"Distinct eQTL"
      }
      if("yes"%in%issamesnp){
        bothsamereg<-"Shared eQTL"
      }
    }
    
    
    
    ##say that if regulated by same snp then we keep only this result
    if(bothdiffreg=="Distinct eQTL" && bothsamereg=="Shared eQTL"){
      reg<-"Shared eQTL"
    }else if(bothdiffreg=="nobothdiffreg" && bothsamereg=="Shared eQTL"){
      reg<-"Shared eQTL"
    }else if(bothdiffreg=="Distinct eQTL" && bothsamereg=="nobothsamereg"){
      reg<-"Distinct eQTL"
    }
    
    genetab[y,"regulation"]<-reg
  }
  return(genetab$regulation)
}

for(y in seq(from=0,to=1,by=0.05)){
  result<-getPairstatus(y)
  namecol<-paste("regulation",y,sep="_")
  genetab[,namecol]<-result
}

write.table(genetab,file="genetabgenepair04.txt",row.names=FALSE)

### Get gene expression for both copies -- ###

library(tidyverse)

Expression <- read_table2("/mnt/project/MSLab/Marie/2022_Synchrosmolt/QTLTools/Synchro_mt_TPM.bed")    %>% select(-c(1,2,3,5,6))
Genotype <- read_table2("/mnt/project/MSLab/Marie/2022_Synchrosmolt/QTLTools/imputed_0.01_2730_smolts.traw")   %>% select(-c(1,3,4,5,6))
Covariate <-  read_table2("covariate_synchro.txt") 
## fix the name in genotype
names(Genotype)[2:ncol(Genotype)] <- str_remove(names(Genotype)[2:ncol(Genotype)], "0_")
# scale the expression
Expression[, 2:ncol(Expression)] <- t(apply(Expression[, 2:ncol(Expression)], 1, scale))



genepair<-read.table(file="genetabgenepair04.txt",header=TRUE)


genepair1<-Expression[Expression$pid%in%genepair$gene1,]
genepair2<-Expression[Expression$pid%in%genepair$gene2,]



# Pivot genotype table to long format
#geno_table_long_p1 <- genepair1 %>%
#pivot_longer(cols = -SNP, names_to = "id", values_to = "genotype")

#Pivot expression table to long format
expr_table_long1 <- genepair1 %>%
  pivot_longer(cols = -pid, names_to = "id", values_to = "expression")

expr_table_long2 <- genepair2 %>%
  pivot_longer(cols = -pid, names_to = "id", values_to = "expression")

covar_kept<-Covariate[c(3,5),]
id<-colnames(covar_kept)
tcovar<-t(covar_kept)
long_covar<-cbind.data.frame(id,tcovar)
long_covar<-long_covar[-1,]
colnames(long_covar)<-c("id","sex","LL")

long_covar_LL<-long_covar[,c(1,3)]

#longcovarll<-long_covar_LL[long_covar_LL$LL=="LL",] ## 907 individuals



mergedpair1 <- merge(expr_table_long1 , long_covar, by = "id")
mergedpair2 <- merge(expr_table_long2 , long_covar, by = "id")

##filter

mergedpair1_filter<- mergedpair1[mergedpair1$LL=="LL",]
mergedpair2_filter<- mergedpair2[mergedpair2$LL=="LL",]

##can keep more column if wanted
mergedpair1_filter<-mergedpair1_filter[,1:3]
mergedpair2_filter<-mergedpair2_filter[,1:3]


library(dplyr)

# Left join table 1 with table 2 based on gene1 and pid
joined_table1_2 <- left_join(genepair, mergedpair1_filter, by = c("gene1" = "pid"))
joined_table1_2<-joined_table1_2[!is.na(joined_table1_2$expression),]

# Left join result with table 3 based on gene2 and pid
joined_table1_2_3 <- left_join(joined_table1_2, mergedpair2_filter, by = c("id"="id","gene2" = "pid"))
joined_table1_2_3 <- joined_table1_2_3 [!is.na(joined_table1_2_3$expression.y),]

colnames(joined_table1_2_3)[9:10]<-c("Expression_gene1","Expression_gene2")
write.table(joined_table1_2_3,file="expression_gene_pair02.txt",row.names=FALSE)


### -- Analysis and plots -- ###

library(ggplot2)
library(dplyr)
library(tidyr)

gene_pair_tab<-read.table(file="expression_gene_pair02.txt",header=TRUE)

##coefficient

# calculate spearman coefficient for each condition
result <- gene_pair_tab %>% 
  group_by(ID) %>% 
  summarize(spearman_coeff = cor(Expression_gene1, Expression_gene2, method = "spearman"))

genepair<-read.table(file="genetabgenepair04.txt",header=TRUE)


results2<-left_join(result, genepair, by = c("ID"="ID"))

# Define the desired order of levels
desired_order <- c("Shared eQTL", "Distinct eQTL", "Single copy eQTL", "No eQTL")

# Convert regulation_0 to a factor with the desired order
results2$regulation_0 <- factor(results2$regulation_0, levels = desired_order)


df2 <- results2 %>%
  group_by(regulation_0) %>%
  summarise(Mean = mean(spearman_coeff))

df2$regulation_0<- factor(df2$regulation_0, levels = desired_order)

ggplot(results2,aes(x=spearman_coeff,fill=regulation_0))+geom_histogram(color="black")+theme_bw(20)+facet_wrap(~regulation_0,scales = 'free_y',ncol=4)+xlab("Spearman coefficient")+
  scale_fill_manual(values=c("#8CE86D","#F9B27C","#D186D8","#99B8EA"))+geom_vline(data = df2, mapping = aes(xintercept = Mean),color="firebrick3",size=1.5)+coord_flip()




results2<-left_join(result, genepair, by = c("ID"="ID"))

dfcount2 <- results2%>% group_by(regulation_0) %>% 
  summarise(total_count=n(),.groups = 'drop') %>%
  as.data.frame()  %>% 
  mutate(prop = total_count / sum(total_count))

ggplot(dfcount2,aes(x="",y=prop,fill=regulation_0))+geom_bar(stat="identity",width=1)+
  coord_polar("y",start=0)+theme_void()+scale_fill_manual(values=c("#F9B27C","#99B8EA","#8CE86D","#D186D8"))+
  theme(legend.text = element_text(size=25))


### prop:

#Distinct: 29%
#No eQTL: 24%
## Single eQTL: 42%
## Shared eQTL: 5%



### make ks test to see differences between distribution

results3<-as.data.frame(results2)

ks.test(results3[results3$regulation_0=="Shared eQTL","spearman_coeff"],results3[results3$regulation_0=="Distinct eQTL","spearman_coeff"])
ks.test(results3[results3$regulation_0=="Shared eQTL","spearman_coeff"],results3[results3$regulation_0=="Single copy eQTL","spearman_coeff"])
ks.test(results3[results3$regulation_0=="Shared eQTL","spearman_coeff"],results3[results3$regulation_0=="No eQTL","spearman_coeff"])
ks.test(results3[results3$regulation_0=="Distinct eQTL","spearman_coeff"],results3[results3$regulation_0=="Single copy eQTL","spearman_coeff"])
ks.test(results3[results3$regulation_0=="Distinct eQTL","spearman_coeff"],results3[results3$regulation_0=="No eQTL","spearman_coeff"])
ks.test(results3[results3$regulation_0=="Single copy eQTL","spearman_coeff"],results3[results3$regulation_0=="No eQTL","spearman_coeff"])


ks.test(abs(results3[results3$regulation_0=="Shared eQTL","spearman_coeff"]),abs(results3[results3$regulation_0=="Distinct eQTL","spearman_coeff"]))
ks.test(abs(results3[results3$regulation_0=="Shared eQTL","spearman_coeff"]),abs(results3[results3$regulation_0=="Single copy eQTL","spearman_coeff"]))
ks.test(abs(results3[results3$regulation_0=="Shared eQTL","spearman_coeff"]),abs(results3[results3$regulation_0=="No eQTL","spearman_coeff"]))
ks.test(abs(results3[results3$regulation_0=="Distinct eQTL","spearman_coeff"]),abs(results3[results3$regulation_0=="Single copy eQTL","spearman_coeff"]))
ks.test(abs(results3[results3$regulation_0=="Distinct eQTL","spearman_coeff"]),abs(results3[results3$regulation_0=="No eQTL","spearman_coeff"]))
ks.test(abs(results3[results3$regulation_0=="Single copy eQTL","spearman_coeff"]),abs(results3[results3$regulation_0=="No eQTL","spearman_coeff"]))


x<-abs(results3[results3$regulation_0=="Shared eQTL","spearman_coeff"])
x<-abs(results3[results3$regulation_0=="Distinct eQTL","spearman_coeff"])
x<-abs(results3[results3$regulation_0=="Single copy eQTL","spearman_coeff"])
x<-abs(results3[results3$regulation_0=="No eQTL","spearman_coeff"]) 

mean(x)
sd(x)/sqrt(length(x))