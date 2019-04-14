# ANEVA-DOT
ANEVA Dosage Outlier Test

Author: Pejman Mohammadi (firstname@Scripps.edu), Christina Sousa (sousa.19@osu.edu)

This test is designed to detect if ASE data reveals sufficient imbalance to induce an outlier in total gene expression. The test takes in in two vectors of ASE count data from an R data frame, as well as a vector of population estimates for standard deviations (sample estimates are provided). The output is a data frame containing raw and adjusted p-values using Benjamini-Hoschberg method, and a plot showing outlier data points in red.

To get started, here's an example:

```r
library(ANEVADOT)

# Define the output columns
output_columns <- c("GENE_ID", "TISSUE_ID",  "REF_COUNT", "ALT_COUNT", "TOTAL_COUNT", "NULL_RATIO")

# Load the Sample ASE data
load("data/sample_ASE_data.rda")

# re organize the tables by:
# 1: Selecting only genes that have Vg scores available
# 2: Reordering ASE data and Vg scores so they align
tiss <- "MSCLSK" # The data comes from a skeletal muscle sample
covered_genes <- intersect(Vg_GTEx_v7$IDs, sample_ASE$GENE_ID)
covered_gene_Vgs <- Vg_GTEx_v7[match(covered_genes, Vg_GTEx_v7$IDs), tiss] 
covered_gene_ASE_data <- sample_ASE[match(covered_genes, sample_ASE$GENE_ID),]

# Take the square root of the Vg scores to the get the Standard Deviation (SDg)
covered_gene_SDgs <- sqrt(covered_gene_Vgs) 

# Run ANEVA-DOT
ANEVADOT_scores <- ANEVA_DOT(covered_gene_ASE_data, output_columns = output_columns, 
                          eh1 = "REF_COUNT", eh2 = "ALT_COUNT", coverage = 10, 
                          r0 = covered_gene_ASE_data$NULL_RATIO,
                          Eg_std = covered_gene_SDgs, plot = T)
```

BugReports: `https://github.com/PejLab/ANEVA-DOT/issues`
