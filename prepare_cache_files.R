cohort_index <- read.table("./data/cohort_codes.csv", 
                           row.names = 1, sep = ",", header = TRUE,
                           stringsAsFactors = FALSE)

cohort_names <- rownames(cohort_index)
names(cohort_names) <- sprintf("%s: %s", rownames(cohort_index), cohort_index[,1])

tumor_ranks_10 <- c("TCGA-Kidney2" = 10, 
                    "TCGA-KIRC2" = 8, "TCGA-KIRP2" = 6,  
                    "TCGA-LUAD2" = 10,"TCGA-Lung2" = 10, "TCGA-LUSC2" = 9,
                    "TCGA-DLBC2" = 3, "TCGA-LAML2"=9, "TCGA-DLBCLAML2"=10
                    )


#tumor_ranks_50 <- c( "TCGA-ACC" =  20, "TCGA-BLCA" = 250, "TCGA-BRCA" = 250, "TCGA-CESC" = 250, "TCGA-CHOL" =  10, 
#                    "TCGA-COAD" = 250, "TCGA-DLBC" =  10, "TCGA-ESCA" =  30, "TCGA-HNSC" = 250, "TCGA-KICH" =  10, 
#                    "TCGA-KIRC" = 200, "TCGA-KIRP" = 250, "TCGA-LAML" =  30,  "TCGA-LGG" = 250, "TCGA-LIHC" = 250,
#                    "TCGA-LUAD" = 250, "TCGA-LUSC" = 250, "TCGA-MESO" =  20,   "TCGA-OV" = 250, "TCGA-PAAD" =  50,
#                    "TCGA-PCPG" =  50, "TCGA-PRAD" = 250, "TCGA-READ" =  30, "TCGA-SARC" = 250, "TCGA-SKCM" =  20,
#                    "TCGA-STAD" = 250, "TCGA-TGCT" =  20, "TCGA-THCA" = 250, "TCGA-THYM" =  20, "TCGA-UCEC" = 250,
#                    "TCGA-UCS"  =  10,  "TCGA-UVM" =  10)

#tumor_ranks_100 <- c( "TCGA-ACC" =  20, "TCGA-BLCA" = 250, "TCGA-BRCA" = 150, "TCGA-CESC" = 250, "TCGA-CHOL" =  10, 
#                     "TCGA-COAD" =  50, "TCGA-DLBC" =  10, "TCGA-ESCA" =  40, "TCGA-HNSC" = 250, "TCGA-KICH" =  10,
#                     "TCGA-KIRC" = 250, "TCGA-KIRP" = 250, "TCGA-LAML" =  20,  "TCGA-LGG" = 250, "TCGA-LIHC" = 250, 
#                     "TCGA-LUAD" = 250, "TCGA-LUSC" = 250, "TCGA-MESO" =  10,   "TCGA-OV" = 250, "TCGA-PAAD" =  30, 
#                     "TCGA-PCPG" =  40, "TCGA-PRAD" = 200, "TCGA-READ" =  30, "TCGA-SARC" = 250, "TCGA-SKCM" =  20,
#                     "TCGA-STAD" = 250, "TCGA-TGCT" =  20, "TCGA-THCA" = 250, "TCGA-THYM" =  10, "TCGA-UCEC" = 250, 
#                     "TCGA-UCS"  =  10,  "TCGA-UVM" =  10)

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

tumor_path <- "G:/My Drive/data_201903"
result_path <- sprintf("/results_%d", filter_ratio)

heatmap_mirna_max_size <- 10
heatmap_mrna_max_size <- 50

construct_pathway_dictionary <- function(gene_list) {
  pathways <- union(load_pathways("h.all"), union(load_pathways("c2.cp.kegg"), load_pathways("c2.cp.pid")))
  pathway_dictionary <- matrix(0, length(pathways), length(gene_list))
  rownames(pathway_dictionary) <- sort(unlist(sapply(pathways, function(x) x$name)))
  colnames(pathway_dictionary) <- gene_list
  for (pathway in pathways) {
    pathway_dictionary[pathway$name, intersect(pathway$symbols, gene_list)] <- 1
  }
  return(pathway_dictionary)
}

load_association_databases <- function() {
  mir2cancer <- read.csv("data/mircancer/miRCancerDecember2017.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  mir2cancer_table <- unique(cbind("miRNA" = tolower(mir2cancer$mirId),
                                   "Cancer(s)" = mir2cancer$Cancer,
                                   "Database" = "miR2Cancer"))
  mirna_names <- unique(mir2cancer_table[,"miRNA"])
  mir2cancer_table <- unique(cbind("miRNA" = mirna_names,
                                   "Cancer(s)" = sapply(mirna_names, FUN = function(mirna_name) {paste0(sort(mir2cancer_table[which(mir2cancer_table[, "miRNA"] == mirna_name), "Cancer(s)"]), collapse = ", ")}, USE.NAMES = FALSE),
                                   "Database" = "miR2Cancer"))
  mir2cancer_table <- mir2cancer_table[order(mir2cancer_table[, "miRNA"]),]
  return(mir2cancer_table)
}

load_pathways <- function(family) {
  symbols_lines <- read.table(sprintf("data/msigdb/%s.v6.1.symbols.gmt", family), header = FALSE, sep = ",", stringsAsFactor = FALSE)
  pathways <- vector("list", dim(symbols_lines)[1])
  for (i in 1:dim(symbols_lines)[1]) {
    symbols_entries <- strsplit(symbols_lines[i, 1], "\t")
    pathways[[i]]$name <- symbols_entries[[1]][1]
    pathways[[i]]$link <- symbols_entries[[1]][2]
    pathways[[i]]$symbols <- sort(symbols_entries[[1]][-2:-1])
  }
  return(pathways)
}

association_table <- load_association_databases()

dir.create(sprintf("./data/cache_%d", filter_ratio), showWarnings = FALSE)
for (selected_cohort in cohort_names) {
  dir.create(sprintf("./data/cache_%d/%s", filter_ratio, selected_cohort), showWarnings = FALSE)
  
  load(sprintf("%s/%s.RData", tumor_path, selected_cohort))
  common_patients <- intersect(rownames(TCGA$mirna), rownames(TCGA$mrna))
  TCGA$mirna <- TCGA$mirna[common_patients,]
  TCGA$mirna <- TCGA$mirna[,which(colMeans(TCGA$mirna > 0) >= (filter_ratio / 100))]
  TCGA$mirna <- scale(log2(TCGA$mirna + 1))
  TCGA$mrna <- TCGA$mrna[common_patients,]
  TCGA$mrna <- TCGA$mrna[,which(colMeans(TCGA$mrna > 0) >= (filter_ratio / 100))]
  TCGA$mrna <- scale(log2(TCGA$mrna + 1))
  TCGA$clinical <- TCGA$clinical[common_patients,]
  rownames(TCGA$clinical) <- common_patients
  TCGA$mutation <- NULL
  tcga_cohort <- TCGA
  save(tcga_cohort, file = sprintf("./data/cache_%d/%s/%s.RData", filter_ratio, selected_cohort, selected_cohort))

  survival_matrix <- matrix(NA, nrow(tcga_cohort$clinical), 2)
  rownames(survival_matrix) <- rownames(tcga_cohort$clinical)
  colnames(survival_matrix) <- c("time", "status")
  alive_patients <- which(tcga_cohort$clinical[,"vital_status"] == "Alive")
  dead_patients <- which(tcga_cohort$clinical[,"vital_status"] == "Dead")
  survival_matrix[alive_patients, "status"] <- 0
  survival_matrix[dead_patients, "status"] <- 1
  survival_matrix[alive_patients, "time"] <- sapply(alive_patients, function(patient) {max(tcga_cohort$clinical[patient, "days_to_last_followup"], tcga_cohort$clinical[patient, "days_to_last_known_alive"], na.rm = TRUE)})
  survival_matrix[dead_patients, "time"] <- tcga_cohort$clinical[dead_patients, "days_to_death"]
  survival_matrix[which(is.infinite(survival_matrix[,"time"]) == TRUE), "time"] <- NA
  save(survival_matrix, file = sprintf("./data/cache_%d/%s/%s_survival_matrix.RData", filter_ratio, selected_cohort, selected_cohort))

  #pathway_dictionary <- construct_pathway_dictionary(colnames(summary$Wy))
  pathway_dictionary <- construct_pathway_dictionary(colnames(tcga_cohort$mrna))
  save(pathway_dictionary, file = sprintf("./data/cache_%d/%s/%s_pathway_dictionary.RData", filter_ratio, selected_cohort, selected_cohort))

  load(sprintf("./%s/%s_summary.RData", result_path, selected_cohort))
  result_summary <- summary
  
  module_weights_table <- data.frame("Module" = 1:length(result_summary$weights_normalized),
                                     "Weight" = result_summary$weights_normalized, 
                                     check.names = FALSE, row.names = NULL, stringsAsFactors = FALSE)
  save(module_weights_table, file = sprintf("./data/cache_%d/%s/%s_module_weights_table.RData", filter_ratio, selected_cohort, selected_cohort))
  
  for (selected_module in 1:tumor_ranks[selected_cohort]) {
    mirna_weights_table <- data.frame("miRNA" = rownames(result_summary$Wx_normalized),
                                      "Weight" = result_summary$Wx_normalized[,selected_module],
                                      "Selected" = "No",
                                      check.names = FALSE, row.names = NULL, stringsAsFactors = FALSE)
    mirna_weights_table[abs(mirna_weights_table[,"Weight"]) > +2 / sqrt(nrow(result_summary$Wx_normalized)), "Selected"] <- "Yes"
    mirna_weights_table <- mirna_weights_table[order(abs(mirna_weights_table[,"Weight"]), decreasing = TRUE),]
    save(mirna_weights_table, file = sprintf("./data/cache_%d/%s/%s_mirna_weights_table_module_%d.RData", filter_ratio, selected_cohort, selected_cohort, selected_module))

    mrna_weights_table <- data.frame("mRNA" = colnames(result_summary$Wy_normalized), 
                                     "Weight" = result_summary$Wy_normalized[selected_module,],
                                     "Selected" = "No",
                                     check.names = FALSE, row.names = NULL, stringsAsFactors = FALSE)
    mrna_weights_table[abs(mrna_weights_table[,"Weight"]) > +2 / sqrt(ncol(result_summary$Wy_normalized)), "Selected"] <- "Yes"
    mrna_weights_table <- mrna_weights_table[order(abs(mrna_weights_table[,"Weight"]), decreasing = TRUE),]
    save(mrna_weights_table, file = sprintf("./data/cache_%d/%s/%s_mrna_weights_table_module_%d.RData", filter_ratio, selected_cohort, selected_cohort, selected_module))
 #   for (selected_module in 1:tumor_ranks[selected_cohort]) { 
    selected_mirna_names <- mirna_weights_table[which(mirna_weights_table[,"Selected"] == "Yes"), "miRNA"]
    selected_mrna_names <- mrna_weights_table[which(mrna_weights_table[,"Selected"] == "Yes"), "mRNA"]

    success_in_sample <- rowSums(pathway_dictionary[,selected_mrna_names])
    success_in_background <- rowSums(pathway_dictionary)
    failure_in_background <- ncol(pathway_dictionary) - success_in_background
    sample_size <- length(selected_mrna_names)
    p_values <- phyper(q = success_in_sample - 1, m = success_in_background, n = failure_in_background, k = sample_size, lower.tail= FALSE)
    enrichment_result_table <- data.frame("Pathway/gene set" = rownames(pathway_dictionary),
                                          "Overlap ratio" = sprintf("%s / %s", success_in_sample, success_in_background),
                                          "p-value" = p_values, 
                                          check.names = FALSE, row.names = NULL, stringsAsFactors = FALSE)
    save(enrichment_result_table, file = sprintf("./data/cache_%d/%s/%s_enrichment_result_table_module_%d.RData", filter_ratio, selected_cohort, selected_cohort, selected_module))
#}
    heatmap_mirna_count <- min(length(selected_mirna_names), heatmap_mirna_max_size)
    heatmap_mirna_weights <- mirna_weights_table[1:heatmap_mirna_count, "Weight"]
    names(heatmap_mirna_weights) <- mirna_weights_table[1:heatmap_mirna_count, "miRNA"]
    heatma_mirna_weights <- heatmap_mirna_weights / max(abs(heatmap_mirna_weights))
    save(heatmap_mirna_weights, file = sprintf("./data/cache_%d/%s/%s_heatmap_mirna_weights_module_%d.RData", filter_ratio, selected_cohort, selected_cohort, selected_module))
    heatmap_mrna_count <- min(length(selected_mrna_names), heatmap_mrna_max_size)
    heatmap_mrna_weights <- mrna_weights_table[1:heatmap_mrna_count, "Weight"]
    names(heatmap_mrna_weights) <- mrna_weights_table[1:heatmap_mrna_count, "mRNA"]
    heatmap_mrna_weights <- heatmap_mrna_weights / max(abs(heatmap_mrna_weights))
    save(heatmap_mrna_weights, file = sprintf("./data/cache_%d/%s/%s_heatmap_mrna_weights_module_%d.RData", filter_ratio, selected_cohort, selected_cohort, selected_module))
    
    mirna_expression <- tcga_cohort$mirna[,selected_mirna_names,drop = FALSE]
    mrna_expression <- tcga_cohort$mrna[,selected_mrna_names,drop = FALSE]
    
    patient_clustering <- kmeans(mrna_expression, centers = 2, nstart = 100)$cluster
    save(patient_clustering, file = sprintf("./data/cache_%d/%s/%s_patient_clustering_module_%d.RData", filter_ratio, selected_cohort, selected_cohort, selected_module))
    
    patient_order <- order(patient_clustering)
    
    heatmap_expression <- rbind(t(mirna_expression[patient_order, 1:heatmap_mirna_count, drop = FALSE]), t(mrna_expression[patient_order, 1:heatmap_mrna_count, drop = FALSE]))
    heatmap_expression[heatmap_expression > +3] <- +3
    heatmap_expression[heatmap_expression < -3] <- -3
    save(heatmap_expression, file = sprintf("./data/cache_%d/%s/%s_heatmap_expression_module_%d.RData", filter_ratio, selected_cohort, selected_cohort, selected_module))
    
    association_result_table <- association_table[match(intersect(selected_mirna_names, association_table[,"miRNA"]), association_table[,"miRNA"], nomatch = NULL),,drop = FALSE]
    save(association_result_table, file = sprintf("./data/cache_%d/%s/%s_association_result_table_module_%d.RData", filter_ratio, selected_cohort, selected_cohort, selected_module))
  }
}

pan_cancer_enrichment_result <- NULL
for (selected_cohort in cohort_names) {
  for (selected_module in 1:min(10, tumor_ranks[selected_cohort])) {
    load(file = sprintf("./data/cache_%d/%s/%s_enrichment_result_table_module_%d.RData", filter_ratio, selected_cohort, selected_cohort, selected_module))
    if (is.null(pan_cancer_enrichment_result) == TRUE) {
      pan_cancer_enrichment_result <- enrichment_result_table[,"p-value",drop=FALSE]
    } else {
      pan_cancer_enrichment_result <- cbind(pan_cancer_enrichment_result, enrichment_result_table[,"p-value",drop=FALSE])
    }
    colnames(pan_cancer_enrichment_result)[ncol(pan_cancer_enrichment_result)] <- sprintf("%s%02d", selected_cohort, selected_module)
  }
}
rownames(pan_cancer_enrichment_result) <- enrichment_result_table[, "Pathway/gene set"]
save(pan_cancer_enrichment_result, file = sprintf("./data/cache_%d/pan_cancer_enrichment_result.RData", filter_ratio))
