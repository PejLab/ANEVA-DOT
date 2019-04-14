setwd("Dropbox/Lappalainen Lab/Projects/ANEVA-DOT/data/")

library(ANEVADOT)
library(tidyverse) 

Vg_GTEx_v7[1:10,1:10]

# load("sample_ase.rda")
# 
# my_tiss_ASE %>%
#   filter(LOW_MAPABILITY == 0 & MAPPING_BIAS_SIM == 0) %>% # quality filter, don't remove genotype warning
#   group_by(GENE_ID) %>%
#   filter(TOTAL_COUNT == max(TOTAL_COUNT)) -> # Use the highest expressed ASEsnp
#   my_input

# Make the output columns
output_columns <- c("GENE_ID", "TISSUE_ID",  "REF_COUNT", "ALT_COUNT", "TOTAL_COUNT", "NULL_RATIO")

#sample_ASE <- as.data.frame(my_input %>% select(output_columns))
load("sample_ASE_data.rda")

# re organize the tables
tiss <- "MSCLSK" # The data comes from a skeletal muscle sample
covered_genes <- intersect(Vg_GTEx_v7$IDs, sample_ASE$GENE_ID)
covered_gene_Vgs <- Vg_GTEx_v7[match(covered_genes, Vg_GTEx_v7$IDs), tiss] 
covered_gene_ASE_data <- sample_ASE[match(covered_genes, sample_ASE$GENE_ID),]

# Apply the correction for variance but DON'T correct for log2 space in aneva Vg scores
# The documentation needs to be updated. 
covered_gene_SDgs <- sqrt(covered_gene_Vgs) 

# Now run ANEVA-DOT
ANEVADOT_scores <- ANEVA_DOT(covered_gene_ASE_data, output_columns = output_columns, 
                          eh1 = "REF_COUNT", eh2 = "ALT_COUNT", coverage = 10, 
                          r0 = covered_gene_ASE_data$NULL_RATIO,
                          Eg_std = covered_gene_SDgs, plot = T)
