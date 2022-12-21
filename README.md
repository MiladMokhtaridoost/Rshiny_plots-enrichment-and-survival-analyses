# Rshiny_plots-enrichment-and-survival-analyses

I applied gene set/pathway analysis and survival analysis and generated all related heatmaps in the following papers, using R shiny plots in this repository. 
Publications list:

[1] Identifying Tissue-and Cohort-Specific RNA Regulatory Modules in Cancer Cells Using Multitask Learning, Cancers, vol. 14, p. 4939, 2022.

[2] An efficient framework to identify key mirna–mrna regulatory modules in cancer, Bioinformatics, vol. 36, no. 26, pp. i592–i600, 2020.

[3] Identifying key miRNA–mRNA regulatory modules in cancer using sparse multivariate factor regression, Proceedings of the 6th International Conference on Machine Learning, Optimization, and Data Science (LOD 2020), pp. 422–433, 2020.

step 1: Copy all four R scripts (prepare_cache_files.R, global.R, server.R and ui.R) in a same folder.

Step 2: Modify cohhorts and the underlying ranks(# of modules) in  prepare_cache_files.R and global.R

step 3: Run prepare_cache_files.R

step 4: Open global.R code and run the shiny application by Run App option in Rstudio

*Open the output in browser for better performance.

The gene set enrichment result, miRNA-mRNA inrteraction heatmap, and Kaplan-Meier survival curve for the selected cohort, module and gene set/pathway will appear in the browser.
