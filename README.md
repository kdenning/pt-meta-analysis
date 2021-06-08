# pt-meta-analysis-prelims

The format of the meta-analysis files is as follows:
- There are 2 data files: the original coded data (prelims data) and the converted data. You can reproduce the converted data file using the "Effect size conversion.RMD" file.
- All functions that we wrote are saved in a separate file (functions/meta_functions.R) and sourced in at the top of each RMD file. 
- The analyses are split into 4 separate documents that correspond to our pre-registered analyses (will eventually be publicaly available):

1. Overall multivariate model (Multivartive MLM Results.RMD)
2. Overlap/merging model (Overlap Results.RMD)
3. Stereotyping model (Stereotyping Results.RMD)
4. Proscoial outcomes categorized as impacting "Interpersonal feelings" (Interpersonal Results.RMD)

- Analyses looking into the data overall, outliers, p-value/p-curve analysis, etc. can be found in the Overall Multivariate Model document.
- Due to the number of studies that had to be included in each Forest Plot (may still need changed), they do not reproduce well in the knitted document. The plots have been saved as images in this repository. The categories with more studies are still difficult to read. The plots are more useful for general trends.
