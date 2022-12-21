shinyServer(function(input, output, session) {
  selected_cohort <<- NULL
  tcga_cohort <<- NULL
  survival_matrix <<- NULL
  module_weights_table <<- NULL

  selected_module <<- NULL
  mirna_weights_table <<- NULL
  mrna_weights_table <<- NULL
  enrichment_result_table <<- NULL
  heatmap_expression <<- NULL
  heatmap_mirna_weights <<- NULL
  heatmap_mrna_weights <<- NULL
  patient_clustering <<- NULL
  association_result_table <<- NULL
  
  selected_gene_set <<- NULL
  selected_cut_off <<- NULL
  pan_cancer_enrichment_result <<- NULL
  pan_cancer_enrichment_result_thresholded <<- NULL

  output$tumor_count <- renderValueBox({
    load_cohort(input$cohort_name)
    load_module(1)
    updateSelectInput(session = session, inputId = "regulatory_module_index", label = "Regulatory module:", choices = 1:tumor_ranks[selected_cohort], selected = 1)
    valueBox(value = nrow(tcga_cohort$mrna), subtitle = "primary tumors", icon = icon("users"), color = "light-blue")
  })

  output$mirna_count <- renderValueBox({
    load_cohort(input$cohort_name)
    valueBox(value = ncol(tcga_cohort$mirna), subtitle = "miRNAs", icon = icon("sign-out-alt"), color = "purple")
  })

  output$mrna_count <- renderValueBox({
    load_cohort(input$cohort_name)
    valueBox(value = ncol(tcga_cohort$mrna), subtitle = "mRNAs", icon = icon("sign-in-alt"), color = "purple")
  })

  output$module_count <- renderValueBox({
    load_cohort(input$cohort_name)
    valueBox(value = tumor_ranks[selected_cohort], subtitle = "regulatory modules", icon = icon(name = "sort-amount-down"), color = "olive")
  })

  draw_module_weights <- function() {
    load_cohort(input$cohort_name)
    par(mar = c(4, 5, 1, 1), oma = c(0, 0, 0, 0))
    barplot(height = module_weights_table[,"Weight"] / sum(module_weights_table[,"Weight"]) * 100, 
            xlab = "Module", ylab = "Weight (%)", cex.lab = 1.25, cex.axis = 1.25, cex.names = 1.25,
            las = 1, names.arg = 1:length(module_weights_table[,"Weight"]), col = "#4292C6")
    abline(h = 100 / length(module_weights_table[,"Weight"]), col = "#EF3B2C", lwd = 2, lty = 1)
  }

  output$module_weights_plot <- renderPlot({
    draw_module_weights()
  })

  output$download_module_weights_plot <- downloadHandler(
    filename = function() {
      sprintf("TCGA-%s_module_weights.pdf", selected_cohort)
    },
    content = function(file) {
      pdf(file = file, width = 8.5 * 1.25, height = 4 * 1.25)
      draw_module_weights()
      dev.off()
    }
  )

  output$download_module_weights_table <- downloadHandler(
    filename = function() {
      sprintf("TCGA-%s_module_weights.csv", selected_cohort)
    },
    content = function(file) {
      write.csv(module_weights_table, file = file, row.names = FALSE)
    }
  )

  draw_rna_weights <- function() {
    load_cohort(input$cohort_name)
    load_module(as.numeric(input$regulatory_module_index))

    layout(matrix(c(1, 2), 1, 2))
    par(mar = c(5, 5, 1, 1), oma = c(0, 0, 0, 0))

    hist(mirna_weights_table[,"Weight"], freq = FALSE, breaks = 100,
         main = "", xlab = "miRNA weight", las = 1,
         cex.lab = 1.25, cex.axis = 1.25, col = "#4292C6")
    abline(v = +2 / sqrt(nrow(mirna_weights_table)), col = "#EF3B2C", lwd = 2, lty = 1)
    abline(v = -2 / sqrt(nrow(mirna_weights_table)), col = "#EF3B2C", lwd = 2, lty = 1)

    hist(mrna_weights_table[,"Weight"], freq = FALSE, breaks = 100,
         main = "", xlab = "mRNA weight", las = 1,
         cex.lab = 1.25, cex.axis = 1.25, col = "#4292C6")
    abline(v = +2 / sqrt(nrow(mrna_weights_table)), col = "#EF3B2C", lwd = 2, lty = 1)
    abline(v = -2 / sqrt(nrow(mrna_weights_table)), col = "#EF3B2C", lwd = 2, lty = 1)
  }

  output$rna_weights_plot <- renderPlot({
    draw_rna_weights()
  })

  output$download_rna_weights_plot <- downloadHandler(
    filename = function() {
      sprintf("TCGA-%s_rna_weights_module_%d.pdf", selected_cohort, selected_module)
    },
    content = function(file) {
      pdf(file = file, width = 8.5 * 1.25, height = 4 * 1.25)
      draw_rna_weights()
      dev.off()
    }
  )

  output$download_mirna_weights_table <- downloadHandler(
    filename = function() {
      sprintf("TCGA-%s_mirna_weights_module_%d.csv", selected_cohort, selected_module)
    },
    content = function(file) {
      write.csv(mirna_weights_table, file = file, row.names = FALSE)
    }
  )

  output$download_mrna_weights_table <- downloadHandler(
    filename = function() {
      sprintf("TCGA-%s_mrna_weights_module_%d.csv", selected_cohort, selected_module)
    },
    content = function(file) {
      write.csv(mrna_weights_table, file = file, row.names = FALSE)
    }
  )

  output$enrichment_results_table <- renderDataTable({
    load_cohort(input$cohort_name)
    load_module(as.numeric(input$regulatory_module_index))

    enrichment_result_table <- enrichment_result_table[which(sapply(enrichment_result_table[,"Pathway/gene set"], FUN = function (x) {strsplit(x, fixed = TRUE, split = "_")[[1]][1]}) %in% input$cancer_specific_gene_sets),]
    enrichment_result_table[,"Pathway/gene set"] <- sprintf("<a href='http://www.broadinstitute.org/gsea/msigdb/cards/%s.html' target='_blank'>%s</a>", enrichment_result_table[,"Pathway/gene set"], enrichment_result_table[,"Pathway/gene set"])
    enrichment_result_table[,"FDR q-value"] <- p.adjust(p = enrichment_result_table[,"p-value"], method = "fdr")
    enrichment_result_table <- enrichment_result_table[order(enrichment_result_table[,"FDR q-value"], enrichment_result_table[,"p-value"], decreasing = FALSE),]
    if (input$show_significant_results == TRUE) {
      enrichment_result_table <- enrichment_result_table[which(enrichment_result_table[,"FDR q-value"] < 0.05),]
    }
    enrichment_result_table[,"p-value"] <- sprintf("%g", enrichment_result_table[,"p-value"])
    enrichment_result_table[,"FDR q-value"] <- sprintf("%g", enrichment_result_table[,"FDR q-value"])
    return(enrichment_result_table)
  }, options = list(searching = FALSE, ordering = FALSE, paging = TRUE, pageLength = 10, lengthChange = FALSE), escape = FALSE)

  output$download_enrichment_results_table <- downloadHandler(
    filename = function() {
      sprintf("TCGA-%s_enrichment_results_module_%d.csv", selected_cohort, selected_module)
    },
    content = function(file) {
      enrichment_result_table <- enrichment_result_table[which(sapply(enrichment_result_table[,"Pathway/gene set"], FUN = function (x) {strsplit(x, fixed = TRUE, split = "_")[[1]][1]}) %in% input$cancer_specific_gene_sets),]
      enrichment_result_table[,"FDR q-value"] <- p.adjust(p = enrichment_result_table[,"p-value"], method = "fdr")
      enrichment_result_table <- enrichment_result_table[order(enrichment_result_table[,"FDR q-value"], enrichment_result_table[,"p-value"], decreasing = FALSE),]
      write.csv(enrichment_result_table, file = file, row.names = FALSE)
    }
  )
  
  draw_expression_heatmap <- function() {
    load_cohort(input$cohort_name)
    load_module(as.numeric(input$regulatory_module_index))

    annotation_frame <- as.data.frame(patient_clustering)
    colnames(annotation_frame) <- c("Group")

    pheatmap(heatmap_expression, cluster_rows = FALSE, cluster_cols = FALSE,
             border_color = NA, gaps_row = c(0, 0, rep(length(heatmap_mirna_weights), 4)), show_colnames = FALSE,
             annotation_col = annotation_frame, annotation_colors = list(Group = survival_colors), annotation_legend = FALSE,
             cellheight = heatmap_cell_height, fontsize = heatmap_cell_height, fontsize_row = c(45 * abs(heatmap_mirna_weights), 14 * abs(heatmap_mrna_weights)))
  }

  output$dynamic_expression_heatmap_plot <- renderUI({
    load_cohort(input$cohort_name)
    load_module(as.numeric(input$regulatory_module_index))
    
    plotOutput(outputId = "expression_heatmap_plot", width = 850, height = heatmap_cell_height * (5 + length(heatmap_mirna_weights) + length(heatmap_mrna_weights)))
  })
  
  output$expression_heatmap_plot <- renderPlot({
    draw_expression_heatmap()
  })

  output$download_expression_heatmap_plot <- downloadHandler(
    filename = function() {
      sprintf("TCGA-%s_expression_heatmap_module_%d.pdf", selected_cohort, selected_module)
    },
    content = function(file) {
      pdf(file = file, width = 8.5 * 1.25, height = heatmap_cell_height * (5 + length(heatmap_mirna_weights) + length(heatmap_mrna_weights)) / 750 * 8.5 * 1.25)
      draw_expression_heatmap()
      dev.off()
    }
  )
  
  draw_survival_analysis <- function() {
    survival_matrix <- cbind(survival_matrix, cluster = patient_clustering)
    survival_matrix <- survival_matrix[which(is.na(survival_matrix[,"status"]) == FALSE & survival_matrix[,"time"] > 0),]
    surv.fit <- survfit(Surv(time, status) ~ cluster, data = as.data.frame(survival_matrix))
    surv.diff <- survdiff(Surv(time, status) ~ cluster, data = as.data.frame(survival_matrix))
    p.val <- 1 - pchisq(surv.diff$chisq, length(surv.diff$n) - 1)

    par(mar = c(4.5, 4.5, 2, 1), oma = c(0, 0, 0, 0))
    plot(surv.fit, lty = c(1, 1), cex = 2, mark = c(15, 19), lwd = 3, main = sprintf("p-value = %.4f", p.val), las = 1, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, xlab = "Days", ylab = "Survival ratio", col = survival_colors) 
    legend(x = "topright", legend = c(sprintf("Group #1 (%d patients)", sum(survival_matrix[,"cluster"] == 1)), sprintf("Group #2 (%d patients)", sum(survival_matrix[,"cluster"]  == 2))), cex = 1.5, pt.cex = 1.5, col = survival_colors, lty = c(1, 1), lwd = 3, pch = c(15, 19))
  }
  
  output$survival_analysis_plot <- renderPlot({
    load_cohort(input$cohort_name)
    load_module(as.numeric(input$regulatory_module_index))
    
    draw_survival_analysis()
  })
  
  output$download_survival_analysis_plot <- downloadHandler(
    filename = function() {
      sprintf("TCGA-%s_survival_analysis_module_%d.pdf", selected_cohort, selected_module)
    },
    content = function(file) {
      pdf(file = file, width = 8.5 * 1.25, height = 6 * 1.25)
      draw_survival_analysis()
      dev.off()
    }
  )
  
  output$disease_associations_table <- renderDataTable({
    load_cohort(input$cohort_name)
    load_module(as.numeric(input$regulatory_module_index))
    
    association_result_table[,"Database"] <- sprintf("<a href='http://mircancer.ecu.edu' target='_blank'>%s</a>", association_result_table[,"Database"])
    return(association_result_table)
  }, options = list(searching = FALSE, ordering = FALSE, paging = TRUE, pageLength = 10, lengthChange = FALSE, autoWidth = TRUE, columnDefs = list(list(targets = 0, width = '15%'), list(targets = 2, width = '15%'))), escape = FALSE)
  
  output$download_disease_associations_table <- downloadHandler(
    filename = function() {
      sprintf("TCGA-%s_disease_associations_module_%d.csv", selected_cohort, selected_module)
    },
    content = function(file) {
      write.csv(association_result_table, file = file, row.names = FALSE)
    }
  )
  
  draw_pan_cancer_enrichment <- function() {
    rownames(pan_cancer_enrichment_result_thresholded) <- sub(pattern = sprintf("%s_", input$pan_cancer_gene_set), replacement = "", x = rownames(pan_cancer_enrichment_result_thresholded))
    ph <- pheatmap(pan_cancer_enrichment_result_thresholded, cluster_rows = TRUE, cluster_cols = TRUE, 
                   border_color = NA, show_colnames = FALSE, color = pan_cancer_colors, legend = FALSE,#, cutree_cols = 3, cutree_rows = 3,
                   cellheight = heatmap_cell_height, fontsize = heatmap_cell_height, fontsize_row = 14, clustering_method = input$clustering_method)
    print(ph)
    print(names(which(cutree(ph$tree_row, 3) == 2)))
    print(names(which(cutree(ph$tree_col, 3) == 2)))
    return(ph)
  }

  output$dynamic_pan_cancer_enrichment_plot <- renderUI({
    load_pan_cancer_results(input$pan_cancer_gene_set, as.numeric(input$fdr_cut_off))

    plotOutput(outputId = "pan_cancer_enrichment_plot", width = 850, height = heatmap_cell_height * (6 + nrow(pan_cancer_enrichment_result_thresholded)))
  })
  
  output$pan_cancer_enrichment_plot <- renderPlot({
    load_pan_cancer_results(input$pan_cancer_gene_set, as.numeric(input$fdr_cut_off))

    draw_pan_cancer_enrichment()
  })
  
  output$download_pan_cancer_enrichment_plot <- downloadHandler(
    filename = function() {
      sprintf("pan_cancer_enrichment_gene_set_%s_clustering_method_%s_cut_off_%g.pdf", selected_gene_set, input$clustering_method, selected_cut_off)
    },
    content = function(file) {
      pdf(file = file, width = 8.5 * 1.25, height = heatmap_cell_height * (6 + nrow(pan_cancer_enrichment_result_thresholded)) / 750 * 8.5 * 1.25)
      draw_pan_cancer_enrichment()
      dev.off()
    }
  )
  
  output$download_pan_cancer_enrichment_table <- downloadHandler(
    filename = function() {
      sprintf("pan_cancer_enrichment_gene_set_%s_clustering_method_%s_cut_off_%g.csv", selected_gene_set, input$clustering_method, selected_cut_off)
    },
    content = function(file) {
      write.csv(pan_cancer_enrichment_result, file = file, row.names = FALSE)
    }
  )
})