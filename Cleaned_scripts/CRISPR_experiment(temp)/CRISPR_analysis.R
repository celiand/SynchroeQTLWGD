#### Celian Diblasi
#### 22/10/2024
#### CRISPR experiment plot

##import data from qpcr
data<-read.table(file="qpcr_data_clean.csv",header=TRUE)

library(ggplot2)

data$chr_gene_factor<-paste("chr",as.character(data$chr_gene),sep="")
data$chr_gene_expression_factor<-paste("chr",as.character(data$Chr_gene_expression),sep="")

data$KOstatut<-ifelse(data$Code_gene_edited==data$Code_gene_expression,"KO","notKO")
data$Code_gene_expression_delta<-paste("??",data$Code_gene_edited,sep="") ### encoding issues here, but it's to paste the "delta" sign with the name
#colnames(data)



ggplot(data[data$Gene_expression=="mterf3",],aes(x=as.factor(Code_gene_expression),y=relative_expression,fill=as.factor(KOstatut)))+geom_boxplot()+facet_wrap(~Code_gene_expression_delta)+
  geom_hline(yintercept=1,color="darkred",linetype = "dashed",size=1.5)+scale_fill_manual(values=c("darkgrey","white"))+theme_bw(15)+xlab("Copy measured")+labs(fill="copy KO")+scale_y_continuous(breaks=c(seq(0.6,1.6,by=0.2)))+ylab("Relative expression")

ggplot(data[data$Gene_expression=="E3",],aes(x=as.factor(Code_gene_expression),y=relative_expression,fill=as.factor(KOstatut)))+geom_boxplot()+facet_wrap(~Code_gene_expression_delta)+
  geom_hline(yintercept=1,color="darkred",linetype = "dashed",size=1.5)+scale_fill_manual(values=c("darkgrey","white"))+theme_bw(15)+xlab("Copy measured")+labs(fill="copy KO")+scale_y_continuous(breaks=c(seq(0.4,1.4,by=0.2)))+ylab("Relative expression")+ylim(0.4,1.4)

### difference of expression
t.test(data[data$Code_gene_edited=="mterf3_14" & data$Code_gene_expression=="mterf3_27","relative_expression" ], mu = 1)
t.test(data[data$Code_gene_edited=="mterf3_27" & data$Code_gene_expression=="mterf3_14","relative_expression" ], mu = 1)

t.test(data[data$Code_gene_edited=="E3_14" & data$Code_gene_expression=="E3_27","relative_expression" ], mu = 1)
t.test(data[data$Code_gene_edited=="E3_27" & data$Code_gene_expression=="E3_14","relative_expression" ], mu = 1)







