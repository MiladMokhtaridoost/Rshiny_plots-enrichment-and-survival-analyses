library(shiny)
library(shinydashboard)
library(survival)
library(pheatmap)

load_cohort <- function(cohort_name) {
  # load TCGA cohort if it is different than loaded cohort
  if (is.null(selected_cohort) == TRUE || selected_cohort != cohort_name) {
    print(sprintf("loading cohort %s", cohort_name))
    # set selected cohort name
    selected_cohort <<- cohort_name
    selected_module <<- NULL
    
    # load TCGA cohort
    load(sprintf("./data/cache_%d/%s/%s.RData", filter_ratio, selected_cohort, selected_cohort))
    tcga_cohort <<- tcga_cohort
    
    # load survival matrix
    load(sprintf("./data/cache_%d/%s/%s_survival_matrix.RData", filter_ratio, selected_cohort, selected_cohort))
    survival_matrix <<- survival_matrix
    
    # load pathway membership indicator matrix for existing mRNAs
    load(sprintf("./data/cache_%d/%s/%s_pathway_dictionary.RData", filter_ratio, selected_cohort, selected_cohort))
    pathway_dictionary <<- pathway_dictionary
    
    # load result table for module weights
    load(sprintf("./data/cache_%d/%s/%s_module_weights_table.RData", filter_ratio, selected_cohort, selected_cohort))
    module_weights_table <<- module_weights_table
  }
}

load_module <- function(module_index) {
  # load regulatory module if it is different than loaded module
  if (is.null(selected_module) == TRUE || selected_module != module_index) {
    print(sprintf("loading module %d for %s", module_index, selected_cohort))
    # set selected module index
    selected_module <<- module_index
    
    # load miRNA table for selected module
    load(sprintf("./data/cache_%d/%s/%s_mirna_weights_table_module_%d.RData", filter_ratio, selected_cohort, selected_cohort, selected_module))
    mirna_weights_table <<- mirna_weights_table
    
    # load mRNA table for selected module
    load(sprintf("./data/cache_%d/%s/%s_mrna_weights_table_module_%d.RData", filter_ratio, selected_cohort, selected_cohort, selected_module))
    mrna_weights_table <<- mrna_weights_table
    
    # load enrichment result table for selected module
    load(sprintf("./data/cache_%d/%s/%s_enrichment_result_table_module_%d.RData", filter_ratio, selected_cohort, selected_cohort, selected_module))
    enrichment_result_table <<- enrichment_result_table
    
    # load miRNA and mRNA gene expression profiles for selected module
    load(sprintf("./data/cache_%d/%s/%s_heatmap_mirna_weights_module_%d.RData", filter_ratio, selected_cohort, selected_cohort, selected_module))
    heatmap_mirna_weights <<- heatmap_mirna_weights
    load(sprintf("./data/cache_%d/%s/%s_heatmap_mrna_weights_module_%d.RData", filter_ratio, selected_cohort, selected_cohort, selected_module))
    heatmap_mrna_weights <<- heatmap_mrna_weights
    
    # load k-means clustering on mRNA expression profiles for selected module
    load(sprintf("./data/cache_%d/%s/%s_patient_clustering_module_%d.RData", filter_ratio, selected_cohort, selected_cohort, selected_module))
    patient_clustering <<- patient_clustering
    
    # load expression matrix to draw heat maps for selected module
    load(sprintf("./data/cache_%d/%s/%s_heatmap_expression_module_%d.RData", filter_ratio, selected_cohort, selected_cohort, selected_module))
    heatmap_expression <<- heatmap_expression
    
    # load association table for selected module
    load(sprintf("./data/cache_%d/%s/%s_association_result_table_module_%d.RData", filter_ratio, selected_cohort, selected_cohort, selected_module))
    association_result_table <<- association_result_table
  }
}

load_pan_cancer_results <- function(gene_set, cut_off) {
  # load pan-cancer analysis results if it is different than loaded gene set
  if (is.null(selected_gene_set) == TRUE || selected_gene_set != gene_set || is.null(selected_cut_off) == TRUE || selected_cut_off != cut_off) {
    print(sprintf("loading pan-cancer analysis results for %s and %g", gene_set, cut_off))
    # set selected gene set
    selected_gene_set <<- gene_set
    selected_cut_off <<- cut_off
    
    # load pan-cancer analysis results
    load(file = sprintf("./data/cache_%d/pan_cancer_enrichment_result.RData", filter_ratio))
    
    pan_cancer_enrichment_result <- pan_cancer_enrichment_result[which(sapply(rownames(pan_cancer_enrichment_result), FUN = function (x) {strsplit(x, fixed = TRUE, split = "_")[[1]][1]}) %in% gene_set),]
    for (c in 2:ncol(pan_cancer_enrichment_result)) {
      pan_cancer_enrichment_result[,c] <- p.adjust(p = pan_cancer_enrichment_result[,c], method = "fdr") 
    }
    pan_cancer_enrichment_result_thresholded <- 1 * (pan_cancer_enrichment_result < selected_cut_off)
    pan_cancer_enrichment_result_thresholded <- pan_cancer_enrichment_result_thresholded[rowSums(pan_cancer_enrichment_result_thresholded) > 0, colSums(pan_cancer_enrichment_result_thresholded) > 0]
    pan_cancer_enrichment_result <- cbind("Pathway/gene set" = rownames(pan_cancer_enrichment_result), pan_cancer_enrichment_result)
    pan_cancer_enrichment_result <<- pan_cancer_enrichment_result
    pan_cancer_enrichment_result_thresholded <<- pan_cancer_enrichment_result_thresholded
  }
}

cohort_index <- read.table("./data/cohort_codes.csv", 
                           row.names = 1, sep = ",", header = TRUE,
                           stringsAsFactors = FALSE)

cohort_names <- rownames(cohort_index)
names(cohort_names) <- sprintf("%s: %s", rownames(cohort_index), cohort_index[,1])

#"TCGA-Kidney" = 10, "TCGA-KIRC" = 9, "TCGA-KIRP" = 5, "TCGA-LUAD" = 10,"TCGA-Lung" = 10, "TCGA-LUSC" = 9,"TCGA-DLBC" = 4, "TCGA-LAML"=8, "TCGA-DLBCLAML"=9,
tumor_ranks_10 <- c(
                    "TCGA-Kidney2" = 10,"TCGA-KIRC2" = 8, "TCGA-KIRP2" = 6,  
                    "TCGA-LUAD2" = 10,"TCGA-Lung2" = 10, "TCGA-LUSC2" = 9,
                    "TCGA-DLBC2" = 3, "TCGA-LAML2"=9, "TCGA-DLBCLAML2"=10,
                    "TCGA-LUAD" = 10,"TCGA-Lung" = 10, "TCGA-LUSC" = 9,
                    "TCGA-DLBC" = 4, "TCGA-LAML"=8, "TCGA-DLBCLAML"=9)

tumor_ranks_20 <- c( "TCGA-ACC" =  20, "TCGA-BLCA" = 20, "TCGA-BRCA" = 20, "TCGA-CESC" = 20, "TCGA-CHOL" =  20,
                     "TCGA-COAD" = 20, "TCGA-DLBC" =  20, "TCGA-ESCA" =  20, "TCGA-HNSC" = 20, "TCGA-KICH" =  20, 
                     "TCGA-KIRC" = 20, "TCGA-KIRP" = 20, "TCGA-LAML" =  20,  "TCGA-LGG" = 20, "TCGA-LIHC" = 20, 
                     "TCGA-LUAD" = 20, "TCGA-LUSC" = 20, "TCGA-MESO" =  20,   "TCGA-OV" = 20, "TCGA-PAAD" =  20, 
                     "TCGA-PCPG" =  20, "TCGA-PRAD" = 20, "TCGA-READ" =  20, "TCGA-SARC" = 20, "TCGA-SKCM" =  20,
                     "TCGA-STAD" = 20, "TCGA-TGCT" =  20, "TCGA-THCA" = 20, #"TCGA-THYM" =  20, 
                     "TCGA-UCS"  =  20,  "TCGA-UVM" =  20, "TCGA-UCEC" = 20)

tumor_ranks_50 <- c( "TCGA-ACC" =  20, "TCGA-BLCA" = 250, "TCGA-BRCA" = 250, "TCGA-CESC" = 250, "TCGA-CHOL" =  10, 
                    "TCGA-COAD" = 250, "TCGA-DLBC" =  10, "TCGA-ESCA" =  30, "TCGA-HNSC" = 250, "TCGA-KICH" =  10, 
                    "TCGA-KIRC" = 200, "TCGA-KIRP" = 250, "TCGA-LAML" =  30,  "TCGA-LGG" = 250, "TCGA-LIHC" = 250,
                    "TCGA-LUAD" = 250, "TCGA-LUSC" = 250, "TCGA-MESO" =  20,   "TCGA-OV" = 250, "TCGA-PAAD" =  50,
                    "TCGA-PCPG" =  50, "TCGA-PRAD" = 250, "TCGA-READ" =  30, "TCGA-SARC" = 250, "TCGA-SKCM" =  20,
                    "TCGA-STAD" = 250, "TCGA-TGCT" =  20, "TCGA-THCA" = 250, "TCGA-THYM" =  20, "TCGA-UCEC" = 250,
                    "TCGA-UCS"  =  10,  "TCGA-UVM" =  10)

tumor_ranks_100 <- c( "TCGA-ACC" =  20, "TCGA-BLCA" = 250, "TCGA-BRCA" = 150, "TCGA-CESC" = 250, "TCGA-CHOL" =  10, 
                     "TCGA-COAD" =  50, "TCGA-DLBC" =  10, "TCGA-ESCA" =  40, "TCGA-HNSC" = 250, "TCGA-KICH" =  10,
                     "TCGA-KIRC" = 250, "TCGA-KIRP" = 250, "TCGA-LAML" =  20,  "TCGA-LGG" = 250, "TCGA-LIHC" = 250, 
                     "TCGA-LUAD" = 250, "TCGA-LUSC" = 250, "TCGA-MESO" =  10,   "TCGA-OV" = 250, "TCGA-PAAD" =  30, 
                     "TCGA-PCPG" =  40, "TCGA-PRAD" = 200, "TCGA-READ" =  30, "TCGA-SARC" = 250, "TCGA-SKCM" =  20,
                     "TCGA-STAD" = 250, "TCGA-TGCT" =  20, "TCGA-THCA" = 250, "TCGA-THYM" =  10, "TCGA-UCEC" = 250, 
                     "TCGA-UCS"  =  10,  "TCGA-UVM" =  10)

filter_ratio <- 10
if (filter_ratio == 10) {
  tumor_ranks <- tumor_ranks_10
}
if (filter_ratio == 20) {
  tumor_ranks <- tumor_ranks_20
}
if (filter_ratio == 50) {
  tumor_ranks <- tumor_ranks_50
}
if (filter_ratio == 100) {
  tumor_ranks <- tumor_ranks_100
}

survival_colors <- c(rgb(red = 217 / 255, green = 95 / 255, blue = 2 / 255), 
                     rgb(red = 117 / 255, green = 112 / 255, blue = 179 / 255))

pan_cancer_colors <- c("#e0e0e0", "#b2182b")

heatmap_cell_height <- 12
heatmap_mirna_max_size <- 10
heatmap_mrna_max_size <- 50
