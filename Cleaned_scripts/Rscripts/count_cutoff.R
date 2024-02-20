### Celian Diblasi
### 20/02/2024
### check count depending on PiP cutoff

##check how count change
library(tidyverse)
##first do the same analysis to keep the same genes
gene_pair_tab<-read.table(file="expression_gene_pair02.txt",header=TRUE) ### data from Ohnologs_analysis.R

##coefficient

# calculate spearman coefficient for each condition
result <- gene_pair_tab %>% 
  group_by(ID) %>% 
  summarize(spearman_coeff = cor(Expression_gene1, Expression_gene2, method = "spearman"))

genepair<-read.table(file="genetabgenepair04.txt",header=TRUE) ### data from Ohnologs_analysis.R


results2<-left_join(result, genepair, by = c("ID"="ID"))

##then check count change
genetab<-results2

library(dplyr)

# Extract regulation column names
regulation_columns <- grep("^regulation_", colnames(genetab), value = TRUE)

# Extract unique eQTL types
eqtl_types <- unique(genetab$regulation_0)

# Initialize an empty data frame to store results
summary_table <- data.frame(regulation_number = numeric(0))

# Add columns for each eQTL type
for (eqtl_type in eqtl_types) {
  summary_table[[eqtl_type]] <- integer(0)
}

# Loop through regulation columns
for (reg_col in regulation_columns) {
  regulation_number <- as.numeric(gsub("^regulation_", "", reg_col))
  row_data <- list(regulation_number)
  
  for (eqtl_type in eqtl_types) {
    nb_eqtl_type <- sum(genetab[[reg_col]] == eqtl_type, na.rm = TRUE)
    row_data[[eqtl_type]] <- nb_eqtl_type
  }
  
  # Add a row to the summary table
  summary_table <- rbind(summary_table, do.call(cbind, row_data))
}


library(dplyr)
library(tidyr)


# Reshape the summary_table to the desired format
summary_table_long <- summary_table %>%
  pivot_longer(cols = -V1, names_to = "eQTL_type", values_to = "count")

# Rename the eQTL types to be more human-readable
eqtl_type_names <- c("No eQTL", "Distinct eQTL", "Single copy eQTL", "Shared eQTL")
summary_table_long$eQTL_type <- factor(summary_table_long$eQTL_type, levels = eqtl_type_names)
colnames(summary_table_long)[1]<-"threshold"

summary_table_long$eQTL_type<-as.factor(summary_table_long$eQTL_type)
summary_table_long$eQTL_type<-factor(summary_table_long$eQTL_type,levels=c("Distinct eQTL","No eQTL","Shared eQTL", "Single copy eQTL"))


library(ggplot2)

ggplot(summary_table_long,aes(x=threshold,y=count,color=eQTL_type))+geom_line(size=2)+geom_point(size=4)+scale_color_manual(values=c("#F9B27C","#99B8EA","#8CE86D","#D186D8"))+theme_bw(20)

