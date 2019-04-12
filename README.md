# ANEVA-DOT
ANEVA Dosage Outlier Test

Author: Pejman Mohammadi (firstname@Scripps.edu), Christina Sousa (sousa.19@buckeyemail.osu.edu)

This test is designed to detect if ASE data reveals sufficient imbalance to induce an outlier in total gene expression. The test takes in in two vectors of ASE count data from an R data frame, as well as a vector of population estimates for standard deviations (sample estimates are provided). The output is a data frame containing raw and adjusted p-values using Benjamini-Hoschberg method, and a plot showing outlier data points in red.

BugReports: `https://github.com/PejLab/ANEVA-DOT/issues`
