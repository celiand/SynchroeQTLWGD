### -- check overlap of eQTL with pauline's results -- ##

## first format xpehh data from pauline
xpehhtab<-read.table(file="/mnt/SCRATCH/pabu/bedtools/input2/xp_ehh_euro.bed")
colnames(xpehhtab)<-c("CHROM","POS","END")
xpehhtab$CHROM<-ifelse(xpehhtab$CHROM<10,paste("ssa0",xpehhtab$CHROM,sep=""),paste("ssa",xpehhtab$CHROM,sep=""))

## then format rna data // see script Region_shuffle_bootstrap.R
rna_full_table<-read.table(file="rna_LL_all_august_version_01.txt",header=TRUE)
rna_full_table$CHROM<-ifelse(rna_full_table$snpchrom<10,paste("ssa0",rna_full_table$snpchrom,sep=""),paste("ssa",rna_full_table$snpchrom,sep=""))

## then export data
options("scipen"=100, "digits"=4) ##disable scientific notation for later exportation
write.table(xpehhtab,file="xpehh_edit_pauline.bed",row.names = FALSE,col.names = FALSE,quote = FALSE,sep="\t")

rna_full_table_bed<-data.frame(rna_full_table$CHROM,rna_full_table$snppos,rna_full_table$snppos)
write.table(rna_full_table_bed,file="rna_eQTL.bed",row.names = FALSE,col.names = FALSE,quote = FALSE,sep="\t")


## then  run a bedtools intersect in bash using that command:
# bedtools intersect -a xpehh_edit_pauline.bed -b rna_eQTL.bed -wb > intersect_eQTL_xpehh.txt


## then check the results

tab_intersect<-read.table(file="intersect_eQTL_xpehh.txt")

tab_intersect$snpid<-paste(tab_intersect$V1,tab_intersect$V2,sep="_")
rna_full_table$snpid<-paste(rna_full_table$CHROM,rna_full_table$snppos,sep="_")
subsetrna<-rna_full_table[rna_full_table$snpid%in%tab_intersect$snpid,]




##list genes
## cleavage and polyadenilation specific factor : https://www.ensembl.org/Salmo_salar/Gene/Summary?db=core;g=ENSSSAG00000041135;r=27:16940619-16963211
## choline phosphotransferase : https://www.ensembl.org/Salmo_salar/Gene/Summary?db=core;g=ENSSSAG00000089632;r=7:46103701-46155453
## makorin, ring finger protein : https://www.ensembl.org/Salmo_salar/Gene/Summary?db=core;g=ENSSSAG00000009919;r=17:77305030-77335674
## tuftelin interacting protein : https://www.ensembl.org/Salmo_salar/Gene/Summary?db=core;g=ENSSSAG00000074395;r=24:15578622-15587289;t=ENSSSAT00000134571
## myoneurin like gene: https://www.ensembl.org/Salmo_salar/Gene/Summary?db=core;g=ENSSSAG00000077916;r=29:27073639-27096456;t=ENSSSAT00000144957 / https://www.ncbi.nlm.nih.gov/gene/106590491
## unknown : https://www.ensembl.org/Salmo_salar/Gene/Summary?db=core;g=ENSSSAG00000003771;r=6:47321227-47327207 
## small subunit processome component: https://www.ensembl.org/Salmo_salar/Gene/Summary?db=core;g=ENSSSAG00000099529;r=27:6322783-6332626;t=ENSSSAT00000188575
## proteasome subunit : https://www.ensembl.org/Salmo_salar/Gene/Summary?db=core;g=ENSSSAG00000085574;r=27:10518241-10521084;t=ENSSSAT00000143990
## unknown : https://www.ensembl.org/Salmo_salar/Gene/Summary?db=core;g=ENSSSAG00000086391;r=27:12832006-12846200
## VPS52 subunit of GARP complex : https://www.ensembl.org/Salmo_salar/Gene/Summary?db=core;g=ENSSSAG00000078638;r=27:10789984-10840875 / https://www.ncbi.nlm.nih.gov/gene/6293
## unknown: https://www.ensembl.org/Salmo_salar/Gene/Summary?db=core;g=ENSSSAG00000009103;r=27:11417438-11456842
## creb regulated transcription coactivator : https://www.ensembl.org/Salmo_salar/Gene/Summary?db=core;g=ENSSSAG00000041719;r=27:12402729-12420179 / https://www.ncbi.nlm.nih.gov/gene/106588329
## unknown : https://www.ensembl.org/Salmo_salar/Gene/Summary?db=core;g=ENSSSAG00000077337;r=27:10386765-10391729
## lantibiotic synthetase component: https://www.ensembl.org/Salmo_salar/Gene/Summary?db=core;g=ENSSSAG00000078655;r=29:18443597-18510367
## dis3l2 exoribonuclease : https://www.ensembl.org/Salmo_salar/Gene/Summary?db=core;g=ENSSSAG00000004442;r=14:10851627-10867844 / https://www.ncbi.nlm.nih.gov/gene/129563
## sf3a1 : https://www.ensembl.org/Salmo_salar/Gene/Summary?db=core;g=ENSSSAG00000007482;r=20:21080066-21097383
## Unknown : https://www.ensembl.org/Salmo_salar/Gene/Summary?db=core;g=ENSSSAG00000102113;r=1:29134268-29137370;t=ENSSSAT00000250174
## mediator complex subunit 6 : https://www.ensembl.org/Salmo_salar/Gene/Summary?db=core;g=ENSSSAG00000000096;r=1:5108913-5125037
