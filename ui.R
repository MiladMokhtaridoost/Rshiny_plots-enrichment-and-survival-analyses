shinyUI(
  dashboardPage(
    skin = "blue",
    title = "RNAExplorer: Bayesian inference of cancer-specific miRNA–mRNA regulatory modules",
    header = dashboardHeader(title = "RNAExplorer", titleWidth = 300),
    sidebar = dashboardSidebar(width = 300,
                               sidebarMenu(
                                 id = "analysis_type",
                                 menuItem(
                                   text = "Cancer-specific results",
                                   tabName = "cancer_specific",
                                   selected = TRUE,
                                   startExpanded = TRUE
                                 ),
                                 menuItem(text = "Pan-cancer results",
                                          tabName = "pan_cancer")
                               )),
    body = dashboardBody(
      tabItems(
        tabItem(tabName = "cancer_specific",
                fluidRow(
                  column(
                    width = 9,
                    box(
                      title = selectInput(
                        inputId = "cohort_name",
                        label = "Cohort name:",
                        choices = cohort_names,
                        selectize = TRUE,
                        width = 850
                      ),
                      width = NULL,
                      collapsible = FALSE,
                      solidHeader = TRUE,
                      status = "danger",
                      valueBoxOutput(outputId = "tumor_count", width = 3),
                      valueBoxOutput(outputId = "mirna_count", width = 3),
                      valueBoxOutput(outputId = "mrna_count", width = 3),
                      valueBoxOutput(outputId = "module_count", width = 3)
                    ),
                    box(
                      title = "Module weights",
                      width = NULL,
                      collapsible = TRUE,
                      solidHeader = TRUE,
                      status = "success",
                      downloadButton(outputId = "download_module_weights_plot", label = "Download plot"),
                      downloadButton(outputId = "download_module_weights_table", label = "Download module weights"),
                      HTML("<br>"),
                      plotOutput(
                        outputId = "module_weights_plot",
                        width = 850,
                        height = 400
                      )
                    ),
                    box(
                      title = "Analysis parameters",
                      width = NULL,
                      collapsible = FALSE,
                      solidHeader = TRUE,
                      status = "success",
                      fluidRow(
                        column(
                          width = 3,
                          selectInput(
                            inputId = "regulatory_module_index",
                            label = "Regulatory module:",
                            choices = 1:10,
                            selected = 1
                          )
                        ),
                        column(
                          width = 4,
                          selectizeInput(
                            inputId = "cancer_specific_gene_sets",
                            label = "Pathway/gene set collection(s):",
                            choices = c("HALLMARK", "KEGG", "PID"),
                            selected = c("HALLMARK"),
                            multiple = TRUE
                          ),
                          checkboxInput(
                            inputId = "show_significant_results",
                            label = "Show significant results only",
                            value = TRUE
                          )
                        )
                      )
                    ),
                    box(
                      title = "miRNA and mRNA weights",
                      width = NULL,
                      collapsible = TRUE,
                      solidHeader = TRUE,
                      status = "success",
                      downloadButton(outputId = "download_rna_weights_plot", label = "Download plot"),
                      downloadButton(outputId = "download_mirna_weights_table", label = "Download miRNA weights"),
                      downloadButton(outputId = "download_mrna_weights_table", label = "Download mRNA weights"),
                      HTML("<br>"),
                      plotOutput(
                        outputId = "rna_weights_plot",
                        width = 850,
                        height = 400
                      )
                    ),
                    box(
                      title = "Gene set enrichment analysis",
                      width = NULL,
                      collapsible = TRUE,
                      solidHeader = TRUE,
                      status = "success",
                      downloadButton(outputId = "download_enrichment_results_table", label = "Download enrichment results"),
                      HTML("<br>"),
                      dataTableOutput(outputId = "enrichment_results_table")
                    ),
                    box(
                      title = "miRNA and mRNA expression profiles",
                      width = NULL,
                      collapsible = TRUE,
                      solidHeader = TRUE,
                      status = "success",
                      downloadButton(outputId = "download_expression_heatmap_plot", label = "Download plot"),
                      HTML("<br>"),
                      uiOutput(outputId = "dynamic_expression_heatmap_plot")
                    ),
                    box(
                      title = "Survival analysis",
                      width = NULL,
                      collapsible = TRUE,
                      solidHeader = TRUE,
                      status = "success",
                      downloadButton(outputId = "download_survival_analysis_plot", label = "Download plot"),
                      HTML("<br>"),
                      plotOutput(
                        outputId = "survival_analysis_plot",
                        width = 850,
                        height = 600
                      )
                    ),
                    box(
                      title = "miRNA–disease associations",
                      width = NULL,
                      collapsible = TRUE,
                      solidHeader = TRUE,
                      status = "success",
                      downloadButton(outputId = "download_disease_associations_table", label = "Download disease associations"),
                      HTML("<br>"),
                      dataTableOutput(outputId = "disease_associations_table")
                    )
                  )
                ))
        ,
        tabItem(tabName = "pan_cancer",
                fluidRow(column(
                  width = 9,
                  box(
                    title = "Analysis parameters",
                    width = NULL,
                    collapsible = FALSE,
                    solidHeader = TRUE,
                    status = "success",
                    fluidRow(
                      column(
                        width = 3,
                        selectInput(
                          inputId = "pan_cancer_gene_set",
                          label = "Pathway/gene set collection:",
                          choices = c("HALLMARK", "KEGG", "PID"),
                          selected = "HALLMARK"
                        )
                      ),
                      column(
                        width = 3,
                        selectInput(
                          inputId = "clustering_method",
                          label = "Clustering method:",
                          choices = c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"),
                          selected = "complete"
                        )
                      ),
                      column(
                        width = 3,
                        selectInput(
                          inputId = "fdr_cut_off",
                          label = "FDR q-value cut-off:",
                          choices = c("0.1", "0.05", "0.01", "0.005", "0.001", "0.0005", "0.0001", "0.00005", "0.00001"),
                          selected = "0.001"
                        )
                      )
                    )
                  ),
                  box(
                    title = "Pan-cancer gene set enrichment analysis",
                    width = NULL,
                    collapsible = FALSE,
                    solidHeader = TRUE,
                    status = "success",
                    downloadButton(outputId = "download_pan_cancer_enrichment_plot", label = "Download plot"),
                    downloadButton(outputId = "download_pan_cancer_enrichment_table", label = "Download FDR q-values"),
                    HTML("<br>"),
                    uiOutput(outputId = "dynamic_pan_cancer_enrichment_plot")
                  )
                )))
      ),
      tags$head(
        tags$style(type = "text/css", ".container {max-width: 1400px; width: 1400px;}"),
        tags$style(type = "text/css", ".col-sm-9 {max-width: 900px; width: 900px;}")
      )
    )
  )
)