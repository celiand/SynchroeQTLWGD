#### Marie Saitou - 12/10/2023

library(tidyverse)
#####
df <- read_csv("geneLL_class.csv")
#https://www.dropbox.com/s/776gg2bwmv3evly/geneLL_class.csv?dl=0
df %>%
  mutate(bin = ntile(median_all_TPM, n=10)) %>% 
  dplyr::count(cistrans,bin)  -> df2
ggplot(df2, aes(x = bin, y = n, fill = cistrans)) + 
  geom_bar(position = "fill",stat = "identity") +theme_minimal()+ 
  scale_fill_manual(values = c("orange", "#FCEBB6", "#76c8c8","#e2e2e2"))+
  labs(x = "Gene expression fraction",y="eQTL presence")
#####
