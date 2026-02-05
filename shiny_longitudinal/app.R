### Longitudinal Analysis of Gavage Data - Shiny App

## Startup code
# load libraries
library(shiny)
library(bslib)
library(tidyverse)
library(ggplot2)
library(shinyWidgets)
library(DT)
library(vegan)
library(gratia)
library(pheatmap)

# setwd
setwd("/Users/kristinvandenham/kmvanden/RStudio/microbiome_analysis/shiny_longitudinal")

# load precomputed RDS files
metadata <- readRDS("metadata.rds")
alpha_diversity <- readRDS("alpha_diversity.rds")
beta_diversity <- readRDS("beta_diversity.rds")
differential_abundance <- readRDS("differential_abundance.rds")


# types of diversity metrics
alpha_diversity_metrics <- names(alpha_diversity$models)
beta_diversity_metrics <- names(beta_diversity)

## UI
ui <- fluidPage(theme = bs_theme(version = 5, bootswatch = "flatly", primary = "#2C3E50"), # responsive Bootstrap 5 layout
                titlePanel("Longitudinal Microbiome Analysis"),
                tags$style(HTML(".tab-content-spacing {margin-top: 20px;}")),
                sidebarLayout(
                  sidebarPanel(pickerInput(inputId = "analysis_type", label = "Analysis Type",
                                           choices = c("Alpha Diversity", "Beta Diversity", "Differential Abundance"), selected = "Alpha Diversity" ),
                               
                               # global gavage picker (plots only)
                               pickerInput(inputId = "gavage_plot",
                                           label = "Groups to Display",
                                           choices = unique(metadata$gavage),
                                           selected = unique(metadata$gavage),
                                           multiple = TRUE),
                               helpText("Note: Gavage selection only affects plots. Statistical models are fit using all groups."),
                               
                               ### Alpha diversity dropdowns
                               conditionalPanel(condition = "input.analysis_type == 'Alpha Diversity'",
                                                pickerInput(inputId = "alpha_metric", label = "Alpha Diversity Metric",
                                                             choices = alpha_diversity_metrics, selected = "Observed"), # alpha diversity metric dropdown
                                                 pickerInput(inputId = "model_type", label = "Model Type",
                                                             choices = c("LMM", "GAMM"), selected = "LMM")),
                               
                               ### Beta diversity dropdowns
                               conditionalPanel(condition = "input.analysis_type == 'Beta Diversity'",
                                                pickerInput(inputId = "beta_metric", label = "Beta Diversity Metric",
                                                            choices = beta_diversity_metrics, selected = "bray_curtis")),
                               
                               ### Differential abundance dropdowns
                               conditionalPanel(condition = "input.analysis_type == 'Differential Abundance'",
                                                
                                                # overview type
                                                pickerInput(inputId = "diff_abun_view", label = "Overview Type",
                                                            choices = c("Community Overview", "Individual Taxon"), selected = "Community Overview"),
                                                
                                                # community overview
                                                conditionalPanel(condition = "input.diff_abun_view == 'Community Overview'",
                                                                 
                                                                 # choose method: top N or manual
                                                                 radioButtons(inputId = "comm_taxa_method", label = "Select taxa by:",
                                                                              choices = c("Top N Most Abundant Taxa" = "top_n", "Manual Selection" = "manual"),
                                                                              selected = "top_n", inline = TRUE),
                                                                 
                                                                 # top N input (default in top 10)
                                                                 conditionalPanel(condition = "input.comm_taxa_method == 'top_n'",
                                                                                  numericInput(inputId = "top_n_taxa", label = "Select Number of Taxa to Display",
                                                                                               value = 10, min = 1, max = 25, step = 1)),
                                                                 
                                                                 # manual taxon selection
                                                                 conditionalPanel(condition = "input.comm_taxa_method == 'manual'",
                                                                                  pickerInput(inputId = "comm_taxa_select", label = "Select Taxa to Display",
                                                                                              choices = sort(unique(differential_abundance$gamm$data$taxon)), 
                                                                                              selected = sort(unique(differential_abundance$gamm$data$taxon))[1:10], 
                                                                                              multiple = TRUE, 
                                                                                              options = list(`live-search` = TRUE, `max-options` = 25, `actions-box` = TRUE)))),
                                                
                                                # individual taxon
                                                conditionalPanel(condition = "input.diff_abun_view == 'Individual Taxon'",
                                                                 pickerInput(inputId = "taxon_select", label = "Choose Taxon",
                                                                             choices = sort(unique(differential_abundance$gamm$data$taxon)), # TO UPDATE LATER (WHEN MODELS FINISHED)
                                                                             options = list(`live-search` = TRUE)),
                                                                 pickerInput(inputId = "model_type_indiv", label = "Model Type",
                                                                             choices = c("GAMM"), selected = "GAMM")))),
                  
                  mainPanel(
                    
                    ### Alpha diversity - LMM tabs 
                    conditionalPanel(condition = "input.analysis_type == 'Alpha Diversity' && input.model_type == 'LMM'",
                                     tabsetPanel(tabPanel("Alpha Diversity Metric Over Time Plot", 
                                                          div(class = "tab-content-spacing",
                                                          h4("Alpha Diversity Metric Over Time Plot", 
                                                          actionButton("info_alpha_plot_lmm", label = NULL, icon = icon("question-circle"), 
                                                                       style = "font-size:12px;", class ="btn btn-info btn-sm")),
                                                              plotOutput("alpha_plot_lmm", height = 400))),
                                                 tabPanel("ANOVA (LMM)", 
                                                          div(class = "tab-content-spacing",
                                                          h4("ANOVA (LMM)", 
                                                          actionButton("info_anova_table", label = NULL, icon = icon("question-circle"), 
                                                                       style = "font-size:12px;", class ="btn btn-info btn-sm")),
                                                              tableOutput("anova_table"))),
                                                 tabPanel("Pairwise Comparisons (LMM)", 
                                                          div(class = "tab-content-spacing",
                                                          h4("Pairwise Comparisons (LMM)", 
                                                          actionButton("info_emmeans_table", label = NULL, icon = icon("question-circle"), 
                                                                       style = "font-size:12px;", class ="btn btn-info btn-sm")),
                                                              DTOutput("emmeans_table"))),
                                                 tabPanel("Model Summary (LMM)", 
                                                          div(class = "tab-content-spacing",
                                                          h4("Linear Mixed Model Summary", 
                                                          actionButton("info_lmm_model_summary", label = NULL, icon = icon("question-circle"), 
                                                                       style = "font-size:12px;", class ="btn btn-info btn-sm")),
                                                              verbatimTextOutput("lmm_model_summary"))),
                                                 tabPanel("Model Diagnostics Plots (LMM)",
                                                          div(class = "tab-content-spacing",
                                                          h4("Model Diagnostics Plots (LMM)", 
                                                          actionButton("info_lmm_resid_plot", label = NULL, icon = icon("question-circle"), 
                                                                       style = "font-size:12px;", class ="btn btn-info btn-sm")),
                                                          plotOutput("lmm_resid_plot", height = 400), plotOutput("lmm_qq_plot", height = 400), plotOutput("lmm_resp_plot", height = 400))))),
                    
                    ### Alpha diversity - GAMM tabs 
                    conditionalPanel(condition = "input.analysis_type == 'Alpha Diversity' && input.model_type == 'GAMM'",
                                     tabsetPanel(tabPanel("Alpha Diversity Metric Over Time Plot",
                                                          div(class = "tab-content-spacing",
                                                           h4("Alpha Diversity Metric Over Time Plot", 
                                                          actionButton("info_alpha_plot_gamm", label = NULL, icon = icon("question-circle"), 
                                                                       style = "font-size:12px;", class ="btn btn-info btn-sm")),
                                                          plotOutput("alpha_plot_gamm", height = 400))),
                                                 tabPanel("Smooth Plot (GAMM)", 
                                                          div(class = "tab-content-spacing",
                                                          h4("Smooth Plot (GAMM)", 
                                                          actionButton("info_gamm_smooth_plot", label = NULL, icon = icon("question-circle"), 
                                                                       style = "font-size:12px;", class ="btn btn-info btn-sm")),
                                                              plotOutput("gamm_smooth_plot", height = 400))),
                                                 tabPanel("Model Summary (GAMM)", 
                                                          div(class = "tab-content-spacing",
                                                          h4("Generalized Additive Mixed Model Summary", 
                                                          actionButton("info_gamm_model_summary", label = NULL, icon = icon("question-circle"), 
                                                                      style = "font-size:12px;", class ="btn btn-info btn-sm")),
                                                              verbatimTextOutput("gamm_model_summary"))),
                                                 tabPanel("Smooth Differences (GAMM)", 
                                                          div(class = "tab-content-spacing",
                                                          h4("Smooth Differences Plot (GAMM)", 
                                                          actionButton("info_gamm_smooth_diff_plot", label = NULL, icon = icon("question-circle"), 
                                                                       style = "font-size:12px;", class ="btn btn-info btn-sm")),
                                                          plotOutput("gamm_smooth_diff_plot", height = 800))),
                                                 tabPanel("Model Diagnostics Plots (GAMM)", 
                                                          div(class = "tab-content-spacing",
                                                          h4("Model Diagnostics Plots (GAMM)", 
                                                          actionButton("info_gamm_resid_plot", label = NULL, icon = icon("question-circle"), 
                                                                       style = "font-size:12px;", class ="btn btn-info btn-sm")),
                                                              plotOutput("gamm_resid_plot", height = 400), plotOutput("gamm_qq_plot", height = 400), plotOutput("gamm_resp_plot", height = 400))),
                                                 tabPanel("Concurvity (GAMM)", 
                                                          div(class = "tab-content-spacing",
                                                          h4("Concurvity (GAMM)", 
                                                          actionButton("info_gamm_concurvity", label = NULL, icon = icon("question-circle"), 
                                                                       style = "font-size:12px;", class ="btn btn-info btn-sm")),
                                                              tableOutput("gamm_concurvity"))))),
                    
                    ### Beta diversity tabs
                    conditionalPanel(condition = "input.analysis_type == 'Beta Diversity'",
                                     tabsetPanel(tabPanel("PCoA Plot",
                                                          div(class = "tab-content-spacing",
                                                          h4("Principal Coordinate Analysis Plot", 
                                                          actionButton("info_pcoa_plot", label = NULL, icon = icon("question-circle"), 
                                                                       style = "font-size:12px;", class ="btn btn-info btn-sm")),
                                                          plotOutput("pcoa_plot", height = 400))),
                                                 tabPanel("NMDS Plot",
                                                          div(class = "tab-content-spacing",
                                                          h4("Non-metric Multidimensional Scaling Plot", 
                                                          actionButton("info_nmds_plot", label = NULL, icon = icon("question-circle"), 
                                                                       style = "font-size:12px;", class ="btn btn-info btn-sm")),
                                                          plotOutput("nmds_plot", height = 400))),
                                                 tabPanel("PCA Plot",
                                                          div(class = "tab-content-spacing",
                                                          h4("Principal Component Analysis Plot", 
                                                          actionButton("info_pca_plot", label = NULL, icon = icon("question-circle"), 
                                                                       style = "font-size:12px;", class ="btn btn-info btn-sm")),
                                                          plotOutput("pca_plot", height = 400))),
                                                 tabPanel("Overall PERMANOVA",
                                                          div(class = "tab-content-spacing",
                                                          h4("Overall PERMANOVA", 
                                                          actionButton("info_permanova_overall", label = NULL, icon = icon("question-circle"), 
                                                                       style = "font-size:12px;", class ="btn btn-info btn-sm")),
                                                          tableOutput("permanova_overall"))),
                                                 tabPanel("Time-resolved PERMANOVA", 
                                                          div(class = "tab-content-spacing",
                                                          h4("Time-resolved PERMANOVA", 
                                                          actionButton("info_permanova_time", label = NULL, icon = icon("question-circle"), 
                                                                       style = "font-size:12px;", class ="btn btn-info btn-sm")),
                                                          DTOutput("permanova_time"))),
                                                 tabPanel("Overall Dispersion", 
                                                          div(class = "tab-content-spacing",
                                                          h4("Overall Dispersion", 
                                                          actionButton("info_dispersion_overall", label = NULL, icon = icon("question-circle"), 
                                                                       style = "font-size:12px;", class ="btn btn-info btn-sm")),
                                                          tableOutput("dispersion_overall"))),
                                                 tabPanel("Time-resolved Dispersion", 
                                                          div(class = "tab-content-spacing",
                                                          h4("Time-resolved Dispersion", 
                                                          actionButton("info_dispersion_time", label = NULL, icon = icon("question-circle"), 
                                                                       style = "font-size:12px;", class ="btn btn-info btn-sm")),
                                                          tableOutput("dispersion_time"))),
                                                 tabPanel("Constrained Ordination Plot",
                                                          div(class = "tab-content-spacing",
                                                          h4("Constrained Ordination (dbRDA/RDA) Plot", 
                                                          actionButton("info_con_ord_plot", label = NULL, icon = icon("question-circle"), 
                                                                       style = "font-size:12px;", class ="btn btn-info btn-sm")),
                                                          plotOutput("con_ord_plot", height = 400))),
                                                 tabPanel("Constrained Ordination Permutation Tests",
                                                          div(class = "tab-content-spacing",
                                                          h4("Constrained Ordination Permutation Tests", 
                                                          actionButton("info_con_ord_overall", label = NULL, icon = icon("question-circle"), 
                                                                       style = "font-size:12px;", class ="btn btn-info btn-sm")),
                                                          wellPanel(h5("Overall Permutation Test"), tableOutput("con_ord_overall")),
                                                          wellPanel(h5("Permutation Test by Term"), tableOutput("con_ord_terms")),
                                                          wellPanel(h5("Permutation Test by Axis"), tableOutput("con_ord_axis")))))),
                    
                    
                    ### Differential Abundance - Community Overview
                    conditionalPanel(condition = "input.analysis_type == 'Differential Abundance' && input.diff_abun_view == 'Community Overview'",
                                     tabsetPanel(tabPanel("Heatmap",
                                                          div(class = "tab-content-spacing",
                                                          h4("Heatmap", 
                                                          actionButton("info_rel_abun_heatmap", label = NULL, icon = icon("question-circle"), 
                                                                       style = "font-size:12px;", class ="btn btn-info btn-sm")),
                                                          plotOutput("rel_abun_heatmap", height = 400))),
                                                tabPanel("Stacked Barplot",
                                                         div(class = "tab-content-spacing",
                                                         h4("Stacked Barplot", 
                                                         actionButton("info_rel_abun_barplot", label = NULL, icon = icon("question-circle"), 
                                                                      style = "font-size:12px;", class ="btn btn-info btn-sm")),
                                                         plotOutput("rel_abun_barplot", height = 400))))),
                                     
                    
                    ### Differential Abundance - GAMM tabs 
                    conditionalPanel(condition = "input.analysis_type == 'Differential Abundance' && input.diff_abun_view == 'Individual Taxon' && input.model_type_indiv == 'GAMM'",
                                     tabsetPanel(tabPanel("Relative Abundance Over Time Plot",
                                                          div(class = "tab-content-spacing",
                                                          h4("Relative Abundance Over Time Plot", 
                                                          actionButton("info_gamm_indiv_taxon_rel_abun_plot", label = NULL, icon = icon("question-circle"), 
                                                                       style = "font-size:12px;", class ="btn btn-info btn-sm")),
                                                          plotOutput("gamm_indiv_taxon_rel_abun_plot", height = 400))),
                                                 tabPanel("Smooth Plot (GAMM)", 
                                                          div(class = "tab-content-spacing",
                                                          h4("Smooth Plot (GAMM)", 
                                                          actionButton("info_gamm_indiv_taxon_smooth_plot", label = NULL, icon = icon("question-circle"), 
                                                                       style = "font-size:12px;", class ="btn btn-info btn-sm")),
                                                          plotOutput("gamm_indiv_taxon_smooth_plot", height = 400))),
                                                 tabPanel("Model Summary (GAMM)",
                                                          div(class = "tab-content-spacing",
                                                          h4("Model Summary (GAMM)", 
                                                          actionButton("info_gamm_indiv_taxon_model_summary", label = NULL, icon = icon("question-circle"), 
                                                                       style = "font-size:12px;", class ="btn btn-info btn-sm")),
                                                          wellPanel(h5("Parametric Coefficients"), tableOutput("gamm_model_para_summary")),
                                                          wellPanel(h5("Smooth Terms"), tableOutput("gamm_model_smooth_summary")))),
                                                 tabPanel("Smooth Differences (GAMM)", 
                                                          div(class = "tab-content-spacing",
                                                          h4("Smooth Differences Plot (GAMM)", 
                                                          actionButton("info_gamm_indiv_taxon_smooth_diff_plot", label = NULL, icon = icon("question-circle"), 
                                                                       style = "font-size:12px;", class ="btn btn-info btn-sm")),
                                                          plotOutput("gamm_indiv_taxon_smooth_diff_plot", height = 800))),
                                                 tabPanel("Model Diagnostics Plots (GAMM)", 
                                                          div(class = "tab-content-spacing",
                                                          h4("Model Diagnostics Plots (GAMM)", 
                                                          actionButton("info_gamm_indiv_taxon_resid_plot", label = NULL, icon = icon("question-circle"), 
                                                                       style = "font-size:12px;", class ="btn btn-info btn-sm")),
                                                          plotOutput("gamm_indiv_taxon_resid_plot", height = 400), plotOutput("gamm_indiv_taxon_qq_plot", height = 400), plotOutput("gamm_indiv_taxon_resp_plot", height = 400))),
                                                 tabPanel("Concurvity (GAMM)", 
                                                          div(class = "tab-content-spacing",
                                                          h4("Concurvity (GAMM)", 
                                                          actionButton("info_gamm_indiv_taxon_concurvity", label = NULL, icon = icon("question-circle"), 
                                                                       style = "font-size:12px;", class ="btn btn-info btn-sm")),
                                                          tableOutput("gamm_indiv_taxon_concurvity"))))),
                  )
                )
)

## Server
server <- function(input, output, session) {
  
  ### define global color map for gavage group plots
  gavage_levels <- sort(unique(metadata$gavage))
  
  gavage_colors <- setNames(RColorBrewer::brewer.pal(n = length(gavage_levels), name = "Set1"),
                            gavage_levels)
  
  
  ### define global color map for gavage group plots
  pair_levels <- unique(differential_abundance$gamm$smooth_diffs$pair)
  
  pair_colors <- setNames(RColorBrewer::brewer.pal(n = length(pair_levels), name = "Set1"),
                          pair_levels)
  
  
  ######################################################
  ##########     ALPHA DIVERSITY ANALYSIS     ##########
  ######################################################
  
  ### filter for selected gavage groups for alpha diversity plots
  alpha_plot_df <- reactive({
    alpha_diversity$metrics %>%
      left_join(metadata, by = c("sample_name" = "sample_id")) %>%
      filter(gavage %in% input$gavage_plot)
  })
  
  ### Linear mixed effects model
  
  # plot alpha diversity metric over time (LMM)
  output$alpha_plot_lmm <- renderPlot({
    ggplot(alpha_plot_df(), aes(x = day, y = .data[[input$alpha_metric]], color = gavage)) +
      geom_jitter(width = 0.2, alpha = 0.5) + theme_minimal() +
      stat_summary(aes(group = gavage), fun = mean, geom = "line", linewidth = 0.8) +
      scale_color_manual(values = gavage_colors, drop = FALSE) +
      labs(y = input$alpha_metric, x = "Day", color = "Gavage")
  }) 
  
  # info pop-up for alpha diversity metric over time (LMM)
  observeEvent(input$info_alpha_plot_lmm, {
    showModal(
      modalDialog(
        title = "Alpha Diversity Metric Over Time",
        p(strong("Observed Richness:"), "Measures the total number of unique taxa."),
        p(strong("Shannon Index:"), "Measures how hard it is to predict the identity of an individual randomly drawn from the community. 
          Increasing richness and/or evenness increases the uncertainty of this prediction and thus also increases the Shannon index. 
          The Shannon index is sensitive to rare taxa since they increase richness"),
        p(strong("Simpson Index:"), "Measures the probability that two randomly drawn individuals are from different taxa. 
          Thus, increasing evenness decreases the Simpson index. Abundant taxa have a large effect on the measure, whereas rare taxa contribute very little."),
        p(strong("Chao1:"), "Estimates total species richness. Chao1 infers the number of unseen species by using the observed counts of singletons and doubletons 
          and the assumption that individuals are randomly and independently sampled from the community."),
        hr(),
        p(strong("Linear Mixed Model (LMM):"), "Models both fixed effects (e.g., time and treatment) and random effects (e.g., subject id). LMMs can only model linear 
          relationships between predictors and the response."),
        p(strong("Generalized Additive Mixed Model (GAMM):"), "Models both fixed effects (e.g., time and treatment) and random effects (e.g., subject id). 
          GAMMs allow non-linear relationships between predictors and the response via smooth functions (splines)."),
        easyClose = TRUE,
        footer = modalButton("Close")
      )
    )
  })
  
  
  # LMM ANOVA table
  output$anova_table <- renderTable({
    df <- as.data.frame(alpha_diversity$models[[input$alpha_metric]]$lmm$anova) 
    
    # format table
    df_formatted <- df %>% 
      tibble::rownames_to_column("Model Term") %>%
      mutate(`Sum Sq` = format(round(`Sum Sq`, 4), nsmall = 4),
             `Mean Sq` = format(round(`Mean Sq`, 4), nsmall = 4),
             NumDF = format(round(NumDF, 0), nsmall = 0),
             DenDF = format(round(DenDF, 0), nsmall = 0),
             `F value` = format(round(`F value`, 4), nsmall = 4),
             `Pr(>F)` = format(round(`Pr(>F)`, 5), nsmall = 5))
    
    df_formatted
  })
  
  # info pop-up for the ANOVA table
  observeEvent(input$info_anova_table, {
    showModal(
      modalDialog(
        title = "How to Read the ANOVA Table",
        p("The ANOVA tests whether each term in the model explain variation in the response after accounting 
          for all other terms in the model (i.e., if a term is removed from the full model, does the model ger significanlty worse?)."),
        hr(),
        p(strong("Sum Sq (Sum of squares):"), "Total variance explained by that term, after accounting for all other fixed effects."),
        p(strong("Mean Sq (Sum Sq divided by NumDF):"), "Variance explained per degree of freedom."),
        p(strong("NumDf (Numerator degrees of freedom):"), "How many parameters are being tested for that term."),
        p(strong("DenDf (Denominator degrees of freedom):"), "Effective residual degrees of freedom available for testing fixed effects, 
          accounting for random effects and unbalanced data using Satterhwaite's approximation."),
        p(strong("F:"), "Parametric F-statistic (ratio of signal to noise)."),
        p(strong("Pr(>F):"), "P-value."),
        easyClose = TRUE,
        footer = modalButton("Close")
      )
    )
  })
  
  
  # LMM pairwise comparisons (estimated marginal means table)
  output$emmeans_table <- renderDT({
    df <- as.data.frame(alpha_diversity$models[[input$alpha_metric]]$lmm$emmeans)
    
    # format table
    df_formatted <- df %>%
      mutate(estimate = format(round(estimate, 2), nsmall = 2),
             SE = format(round(SE, 2), nsmall = 2),
             df = format(round(df, 0), nsmall = 0),
             t.ratio = format(round(t.ratio, 3), nsmall = 3),
             p.value = format(round(p.value, 4), nsmall = 4))
    
    df_formatted
  })
  
  # info pop-up for the LMM pairwise comparisons (estimated marginal means table)
  observeEvent(input$info_emmeans_table, {
    showModal(
      modalDialog(
        title = "How to Read the Pairwise Comparisons Table",
        p("Pairwise comparisons of gavage groups at each timepoint using model-based estimated marginal means."),
        p(strong("day_factor:"), "The timepoint at which the comparison was made."),
        p(strong("SE (Standard error):"), "Uncertainty in the estimate difference."),
        p(strong("estimate:"), "Estimated difference in the alpha diversity metric between the two groups based on the mixed model, not on raw means."),
        p(strong("SE (Standard error):"), "Uncertainty in the estimate difference."),
        p(strong("df (Degrees of freedom):"), "Effective degrees of freedom used in the t-test."),
        p(strong("t.ratio:"), "How many standard errord the contrast is away from zero. The larger this value is, the stronger the evidence is. 
          2 = borderline, >3 = strong and > 4 = very strong"),
        p(strong("p.value:"), "Family-wise error controlled (by timepoint) p-values."),
        footer = modalButton("Close")
      )
    )
  })
  
  
  # LMM summary 
  output$lmm_model_summary <- renderPrint({
    summary(alpha_diversity$models[[input$alpha_metric]]$lmm$model)
  })
  
  # info pop-up for the LMM summary
  observeEvent(input$info_lmm_model_summary, {
    showModal(
      modalDialog(
        title = "How to Interpret the LMM Summary",
        p(strong("Scaled residuals:"), "Shows the distribution of the residuals Ideally, residuals should be centered areound zero, with a roughly symmetric spread."),
        p(strong("Random effects:"), "Variance of random effects (effects due to grouping or clustering, that account for correlation within groups)."),
        p(strong("Fixed effects:"), "Whether there is a significant difference between intercept, main effect of gavage, main effect of day and interactions (modification 
          of gavage effect acros days) and baseline."),
        footer = modalButton("Close")
      )
    )
  })
  
  
  # plot LMM residuals versus fitted
  output$lmm_resid_plot <- renderPlot({
    model <- alpha_diversity$models[[input$alpha_metric]]$lmm$model
    plot(fitted(model), resid(model),
         xlab = "Fitted values", ylab = "Residuals",
         main = "LMM Residuals vs Fitted")
    abline(h = 0, lty = 2)
  })
  
  # plot LMM normal Q-Q
  output$lmm_qq_plot <- renderPlot({
    model <- alpha_diversity$models[[input$alpha_metric]]$lmm$model
    qqnorm(resid(model), main = "LMM Normal Q-Q")
    qqline(resid(model))
  })
  
  
  # plot LMM response plot
  output$lmm_resp_plot <- renderPlot({
    model <- alpha_diversity$models[[input$alpha_metric]]$lmm$model
    plot(fitted(model), alpha_diversity$metrics[[input$alpha_metric]],
         ylab = "Response", xlab = "Fitted Values",
         main = "LMM Response vs Fitted Values")
    abline(a = 0, b = 1, lty = 2)
  })
  
  
  # info pop-up for LMM residuals vs fitted and Q-Q plots
  observeEvent(input$info_lmm_resid_plot, {
    showModal(
      modalDialog(
        title = "How to Interpret the Residuals vs Fitted Values Plot, the Normal Q-Q Plot and the Response vs Fitted Values Plot",
        p("The Residuals versus Fitted Values plot is used to check whether the statistical model is appropriate for the data by checking whether the remaining errors (residuals) are random 
          or if they show structure the model failed to capture. The points should be randomly scattered around zero (i.e., no clear pattern or trend)."),
        p(strong("Curved pattern:"), "Model failed to capture nonlinear strucutre in the mean (e.g., missing interaction)."),
        p(strong("Curved clusters or bands:"), "Group-specific trends or correlation structure is not adequately modeled."),
        p(strong("Widening/narrowing:"), "Non-constant variance (heteroscedasticity)."),
        p(strong("Extreme outliers:"), "Influential points or possible errors."),
        hr(),
        p("The Normal Q-Q plot (quantile-quantile plot) is used to compare the distribution of the residuals with a theoretical distribution (often normal). 
          Each point represents how one residual compares to what would be expected under normality. LMMs and Gaussian GAMMs assume that residuals are normally distributed. 
          If residuals are not normal, p-values and confidence intervals may be unreliable."),
        p(strong("Points lie on a straight 45° line:"), "Residuals are approximately normally distributed."),
        p(strong("S-shaped curve:"), "Less extreme values (concave up then down) or more extreme values (concave down then up) than normal."),
        p(strong("Points curve away at ends:"), "Indicates skewed residuals."),
        p(strong("Large deviations at the ends:"), "Potential outliers."),
        hr(),
        p("The Response versus Fitted Values plot compares the observed response values to the values predicted by the model. This plot helps you see how well the model 
          captures the overall trends in the data and whether there are systematic deviations."),
        p(strong("Points lie roughly along the y = x line:"), "The model fits the data well; predicted values match observed values."),
        p(strong("Points systematically above or below the line:"), "The model under- or over-predicts in certain ranges."),
        p(strong("Nonlinear patterns or curves:"), "The model may be missing key nonlinear effects or interactions."),
        p(strong("Large scatter or wide spread around the line:"), "High residual variance; model may not explain much of the variation."),
        p(strong("Clusters or gaps:"), "Indicates group-specific effects or unmodeled structure in the data."),
        footer = modalButton("Close")
      )
    )
  })
  
  
  ### Generalized additive mixed model
  
  # plot alpha diversity metric over time (GAMM)
  output$alpha_plot_gamm <- renderPlot({
    ggplot(alpha_plot_df(), aes(x = day, y = .data[[input$alpha_metric]], color = gavage)) +
      geom_jitter(width = 0.2, alpha = 0.5) + theme_minimal() +
      stat_summary(aes(group = gavage), fun = mean, geom = "line", linewidth = 0.8) +
      scale_color_manual(values = gavage_colors, drop = FALSE) +
      labs(y = input$alpha_metric, x = "Day", color = "Gavage")
  }) 
  
  # info pop-up for alpha diversity metric over time (gamm)
  observeEvent(input$info_alpha_plot_gamm, {
    showModal(
      modalDialog(
        title = "Alpha Diversity Metric Over Time",
        p(strong("Observed Richness:"), "Measures the total number of unique taxa."),
        p(strong("Shannon Index:"), "Measures how hard it is to predict the identity of an individual randomly drawn from the community. 
          Increasing richness and/or evenness increases the uncertainty of this prediction and thus also increases the Shannon index. 
          The Shannon index is sensitive to rare taxa since they increase richness"),
        p(strong("Simpson Index:"), "Measures the probability that two randomly drawn individuals are from different taxa. 
          Thus, increasing evenness decreases the Simpson index. Abundant taxa have a large effect on the measure, whereas rare taxa contribute very little."),
        p(strong("Chao1:"), "Estimates total species richness. Chao1 infers the number of unseen species by using the observed counts of singletons and doubletons 
          and the assumption that individuals are randomly and independently sampled from the community."),
        hr(),
        p(strong("Linear Mixed Model (LMM):"), "Models both fixed effects (e.g., time and treatment) and random effects (e.g., subject id). LMMs can only model linear 
          relationships between predictors and the response."),
        p(strong("Generalized Additive Mixed Model (GAMM):"), "Models both fixed effects (e.g., time and treatment) and random effects (e.g., subject id). 
          GAMMs allow non-linear relationships between predictors and the response via smooth functions (splines)."),
        easyClose = TRUE,
        footer = modalButton("Close")
      )
    )
  })
  
  
  # GAMM smooth plots
  output$gamm_smooth_plot <- renderPlot({
    smooth_df <- alpha_diversity$models[[input$alpha_metric]]$gamm$smooth_df
    ggplot() + theme_minimal() +
      geom_point(data = alpha_plot_df(), aes(x = day, y = .data[[input$alpha_metric]], color = gavage), alpha = 0.4) +
      geom_line(data = smooth_df, aes(x = day, y = fit, color = gavage), linewidth = 1) +
      geom_ribbon(data = smooth_df, aes(x = day, ymin = lower, ymax = upper, fill = gavage), alpha = 0.2, color = NA) +
      scale_color_manual(values = gavage_colors, drop = FALSE) +
      scale_fill_manual(values = gavage_colors, drop = FALSE) +
      labs(y = input$alpha_metric, x = "Day", color = "Gavage", fill = "Gavage")
  })
  
  # info pop-up for smooths plot
  observeEvent(input$info_gamm_smooth_plot, {
    showModal(
      modalDialog(
        title = "GAMM-estimated Smooths Over Time",
        p("Generalized Additive Mixed Models model both fixed effects (e.g., time and treatment) and random effects (e.g., subject id). 
          GAMMs allow non-linear relationships between predictors and the response via smooth functions (splines)."),
        p(strong("GAMM-estimated smooth trajectories:"), "Estimate the relationship between a predictor and a response without assuming a fixed parametric form. 
          Smooths are penalized and therefore only allow curvature if the data strongly supports it."),
        footer = modalButton("Close")
      )
    )
  })
  
  
  # plot GAMM residuals plot
  output$gamm_resid_plot <- renderPlot({
    model <- alpha_diversity$models[[input$alpha_metric]]$gamm$model$gam
    plot(fitted(model), resid(model),
        ylab = "Residuals", xlab = "Fitted Values",
         main = "GAMM Residuals vs Fitted Values Plot")
    abline(h = 0, lty = 2)
  })
  
  # plot GAMM normal Q-Q plot
  output$gamm_qq_plot <- renderPlot({
    res <- alpha_diversity$models[[input$alpha_metric]]$gamm$residuals
    qqnorm(res, 
           ylab = "Deviance Residuals", xlab = "Theoretical Quantiles",
           main = "GAMM Normal Q-Q Plot")
    qqline(res)
  })
  
  
  # plot GAMM response plot
  output$gamm_resp_plot <- renderPlot({
    model <- alpha_diversity$models[[input$alpha_metric]]$gamm$model$gam
    plot(fitted(model), alpha_diversity$metrics[[input$alpha_metric]],
         ylab = "Response", xlab = "Fitted Values",
         main = "GAMM Response vs Fitted Values Plot")
    abline(a = 0, b = 1, lty = 2)
  })
    
    
  # info pop-up for GAM residuals vs fitted plot, normal Q-Q plot and response vs fitted plot
  observeEvent(input$info_gamm_resid_plot, {
    showModal(
      modalDialog(
        title = "How to Interpret the Residuals vs Fitted Values Plot, the Normal Q-Q Plot and the Response vs Fitted Values Plot",
        p("The Residuals versus Fitted Values plot is used to check whether the statistical model is appropriate for the data by checking whether the remaining errors (residuals) are random 
          or if they show structure the model failed to capture. The points should be randomly scattered around zero (i.e., no clear pattern or trend)."),
        p(strong("Curved pattern:"), "Model failed to capture nonlinear strucutre in the mean (e.g., missing interaction)."),
        p(strong("Curved clusters or bands:"), "Group-specific trends or correlation structure is not adequately modeled."),
        p(strong("Widening/narrowing:"), "Non-constant variance (heteroscedasticity)."),
        p(strong("Extreme outliers:"), "Influential points or possible errors."),
        hr(),
        p("The Normal Q-Q plot (quantile-quantile plot) is used to compare the distribution of the residuals with a theoretical distribution (often normal). 
          Each point represents how one residual compares to what would be expected under normality. LMMs and Gaussian GAMMs assume that residuals are normally distributed. 
          If residuals are not normal, p-values and confidence intervals may be unreliable."),
        p(strong("Points lie on a straight 45° line:"), "Residuals are approximately normally distributed."),
        p(strong("S-shaped curve:"), "Less extreme values (concave up then down) or more extreme values (concave down then up) than normal."),
        p(strong("Points curve away at ends:"), "Indicates skewed residuals."),
        p(strong("Large deviations at the ends:"), "Potential outliers."),
        hr(),
        p("The Response versus Fitted Values plot compares the observed response values to the values predicted by the model. This plot helps you see how well the model 
          captures the overall trends in the data and whether there are systematic deviations."),
        p(strong("Points lie roughly along the y = x line:"), "The model fits the data well; predicted values match observed values."),
        p(strong("Points systematically above or below the line:"), "The model under- or over-predicts in certain ranges."),
        p(strong("Nonlinear patterns or curves:"), "The model may be missing key nonlinear effects or interactions."),
        p(strong("Large scatter or wide spread around the line:"), "High residual variance; model may not explain much of the variation."),
        p(strong("Clusters or gaps:"), "Indicates group-specific effects or unmodeled structure in the data."),
        footer = modalButton("Close")
      )
    )
  })
  

  # GAMM summary
  output$gamm_model_summary <- renderPrint({
    alpha_diversity$models[[input$alpha_metric]]$gamm$summary
  })
  
  # info pop-up for the GAMM summary
  observeEvent(input$info_gamm_model_summary, {
    showModal(
      modalDialog(
        title = "How to Interpret the GAMM Summary",
        p(strong("Family:"), "Distribution used when modelling."),
        p(strong("Link function:"), "Whether the data was transformed. Identity = no transformation"),
        p(strong("Parametric coefficients:"), "Represent the linear (fixed) effects in the model. Each group term is the linear difference from the reference group averaged over time."),
        p(strong("Smooth terms:"), "Indicates whether adding nonlinearity significantly improves the model for a given group."),
        p(strong("edf (effective degrees of freedom:"), "How wiggly the smooth is. > 1 = some nonlinearity"),
        p(strong("R-sq. (adj):"), "Proportion of the variance explained by the model."),
        p(strong("Scale est."), "In Gaussian models, the scale parameter is an estimate of the residual variance (i.e., the amount of variation in your response variable that your model cannot explain)."),
        footer = modalButton("Close")
      )
    )
  })
  
  
  # GAMM smooth difference plots 
  output$gamm_smooth_diff_plot <- renderPlot({ 
    smooth_diff_plot_df <- alpha_diversity$models[[input$alpha_metric]]$gamm$smooth_diff 
    draw(smooth_diff_plot_df) & 
      scale_x_continuous(name = "Day")
  })
  
  
  # info pop-up for the GAMM smooth differences plots
  observeEvent(input$info_gamm_smooth_diff_plot, {
    showModal(
      modalDialog(
        title = "How to Interpret the GAMM Smooth Differences Plots",
        p("These plots show the pairwise differences between the estimated smooth trajectories of the selected taxon of the selected alpha 
        diversity metrics for each gavage group (i.e., whether and when one group differs significantly from another over time.)"),
        p(strong("Positive values:"), "The first group in the comparison has higher CLR abundance than the second group at that timepoint."),
        p(strong("Negative values:"), "The first group in the comparison has lower CLR abundance than the second group at that timepoint."),
        p(strong("Shaded regions:"), "Indicate 95% confidence intervals for the differences. If the shaded region does not overlap zero, then 
          the difference is statistically significant at that timepoint."),
        footer = modalButton("Close")
      )
    )
  })
  
 
  # GAMM concurvity
  output$gamm_concurvity <- renderTable({
    df <- as.data.frame(alpha_diversity$models[[input$alpha_metric]]$gamm$concurvity)
    
    # format table
    df_formatted <- df %>%
      tibble::rownames_to_column("Measure") %>%
      mutate(para = format(round(para, 2), nsmall = 2),
             `s(day_c):gavageG_1DMD` = formatC(`s(day_c):gavageG_1DMD`, format = "e", digits = 4),
             `s(day_c):gavageG_4DMD` = formatC(`s(day_c):gavageG_4DMD`, format = "e", digits = 4),
             `s(day_c):gavageG_4W7C` = formatC(`s(day_c):gavageG_4W7C`, format = "e", digits = 4),
             `s(day_c):gavageG_4WMD` = formatC(`s(day_c):gavageG_4WMD`, format = "e", digits = 4))
    
    df_formatted
  })
  
  # info pop-up for concurvity table
  observeEvent(input$info_gamm_concurvity, {
    showModal(
      modalDialog(
        title = "How to Read the Concurvity Table",
        p("How much redundancy or linear dependence exists in the model (i.e., how much one term can explain another term). 
        If two terms are highly dependendent, it is hard for the model to separate their effects. 
        0 = no redundancy, > 0.5 = moderate redundancy, and > 0.9 = very high redundancy."),
        p(strong("Measure:"), "worst = maximum concurvity for the term across all other terms, observed = concurvity calculated 
          directly from the fitted model, and estimate = smoothed/adjusted estimate of concurvity."),
        p(strong("para:"), "The parametric coefficients (see Model Summary (GAMM))."),
        p(strong("s(day):gavage...:"), "The smooth terms/group-specific nonlinear trends over time (see Model Summary (GAMM))."),
        easyClose = TRUE,
        footer = modalButton("Close")
      )
    )
  })
  
  
  ######################################################
  ##########      BETA DIVERSITY ANALYSIS     ##########
  ######################################################
  
  ### filter for selected gavage groups for beta diversity plots
  beta_plot_meta <- reactive({
    metadata %>% filter(gavage %in% input$gavage_plot) 
  })
  
  
  # reactive data.frame for PCoA plot
  beta_pcoa_df <- reactive({
    metric <- input$beta_metric
    ord <- beta_diversity[[metric]]$ordination$pcoa
    df <- as.data.frame(ord$vectors) # extract coordinates
    df$sample_id <- rownames(df)
    
    # filter df to only include sample_ids that are in the filtered metadata
    df <- df %>% filter(sample_id %in% beta_plot_meta()$sample_id)
    
    # merge with filtered metadata
    df <- left_join(df, beta_plot_meta(), by = "sample_id")
    
    # compute variance explained
    var_explained <- ord$values$Relative_eig * 100
    
    # compute centroids per gavage x day
    df_centroids <- df %>%
      group_by(gavage, day_factor) %>%
      summarize(mean_Axis1 = mean(Axis.1),
                mean_Axis2 = mean(Axis.2),
                .groups = "drop") %>%
      arrange(gavage, as.numeric(as.character(day_factor)))
    
    list(df = df, df_centroids = df_centroids, var_explained = var_explained)
  })
  
  # PCoA plot
  output$pcoa_plot <- renderPlot({
    pcoa_data <- beta_pcoa_df()
    
    ggplot(pcoa_data$df, aes(x = Axis.1, y = Axis.2)) +
      geom_path(data = pcoa_data$df_centroids,
                aes(x = mean_Axis1, y = mean_Axis2, group = gavage, color = gavage),
                arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.5) +
      geom_point(aes(x = Axis.1, y = Axis.2, color = gavage, size = day_factor), 
                 shape = 21, stroke = 0.6, fill = NA, alpha = 1) +
      geom_point(aes(x = Axis.1, y = Axis.2, fill = gavage, size = day_factor), 
                 shape = 21, stroke = 0, alpha = 0.35) +
      stat_ellipse(aes(color = gavage, group = gavage), type = "norm", level = 0.95, linewidth = 0.5) +
      # ellipses assume multivariate normality and are thus only provided for visualization
      scale_color_manual(values = gavage_colors, drop = TRUE, name = "Gavage") +
      scale_fill_manual(values = gavage_colors, drop = TRUE, name = "Gavage") +
      scale_size_manual(values = c("7" = 1, "14" = 2, "21" = 3, "28" = 4), name = "Day") +
      xlab(paste0("PC1 (", round(pcoa_data$var_explained[1], 1), "%)")) +
      ylab(paste0("PC2 (", round(pcoa_data$var_explained[2], 1), "%)")) +
      theme_minimal()
  })
  
  # info pop-up for PCoA Plot
  observeEvent(input$info_pcoa_plot, {
    showModal(
      modalDialog(
        title = "Principal Coordinate Analysis Plot",
        p("PCoA is an unconstrained ordination method that represents samples in a low-dimensional space while 
          best preserving the pairwise distances or dissimilarities between samples."), 
        p("Each axis (principal coordinate) represents the amount of variation in the distance matrix explained by that axis."),
        p("PCoA is appropriate for both Euclidean and non-Euclidean distance metrics and is often the preferred 
          ordination method when distances have a meaningful metric interpretation."),
        p(strong("Bray-Curtis distance:"), "Measures the compositional dissimilarity between two samples based on the abundances of taxa present in at least one of the samples. 
          It is computed as the weighted sum of absolute differences, where weights are the abundances, thus the metric is dominated by abundant taxa."),
        p(strong("Canberra distance:"), "Measures compositional dissimilarity between two samples based on the abundances of taxa present in at least one of the samples. 
          It sums the ratios of the absolute differences to the sums of the abundances for each taxon, thereby giving relatively equal weight to rare and abundant taxa 
          (i.e., rare taxa contribute proportionally more than they would in Bray-Curtis)."),
        p(strong("Jaccard distance:"), "Measures the dissimilarity between two samples based on presence/absence of taxa, ignoring their abundances. 
          It quantifies how many taxa are shared versus unique to each sample."),
        p(strong("Aitchison distance:"), "The Euclidean (straight-line) distance between two samples after centered log ratio (CLR) transformation. CLR transformation maps 
        compositional microbiome data from the simplex into Euclidean space by removing the constant-sum constraint, making Euclidean distance geometrically meaningful."),
        easyClose = TRUE,
        footer = modalButton("Close")
      )
    )
  })
  
  
  # reactive data.frame for NMDS plot
  beta_nmds_df <- reactive({
    metric <- input$beta_metric
    ord <- beta_diversity[[metric]]$ordination$nmds
    df <- as.data.frame(ord$points) # extract coordinates
    df$sample_id <- rownames(df)
    
    # filter df to only include sample_ids that are in the filtered metadata
    df <- df %>% filter(sample_id %in% beta_plot_meta()$sample_id)
    
    # merge with filtered metadata
    df <- left_join(df, beta_plot_meta(), by = "sample_id")
    
    # compute centroids per gavage x day
    df_centroids <- df %>%
      group_by(gavage, day_factor) %>%
      summarize(mean_MDS1 = mean(MDS1),
                mean_MDS2 = mean(MDS2),
                .groups = "drop") %>%
      arrange(gavage, as.numeric(as.character(day_factor)))
    
    list(df = df, df_centroids = df_centroids)
  })
  
  # NMDS plot
  output$nmds_plot <- renderPlot({
    
    # only allow non-Euclidean distance metrics to be used for NMDS plot
    metric <- input$beta_metric
    ord <- beta_diversity[[metric]]$ordination$nmds
    
    validate(
      need(!is.null(ord),
           "NMDS is only available for non-Euclidean distance metrics (Bray-Curtis, Jaccard and Canberra). Euclidean distances have closed-form ordination solutions (i.e., PCA), so applying an iterative, rank-based method like NMDS would unnecessarily discard metric and variance information."))
    
    # add stress tests subtitle
    stress_value  <- ord$stress 
    
    # define a simple interpretation
    stress_label <- if(stress_value < 0.05) {
      "excellent representation"
    } else if(stress_value < 0.1) {
      "good representation"
    } else if(stress_value < 0.2) {
      "fair representation"
    } else if(stress_value < 0.3) {
      "poor representation"
    } else {
      "unreliable representation"
    }
    
    stress_text <- paste0("stress = ", round(stress_value, 4), " (", stress_label, ")")
    
    nmds_data <- beta_nmds_df()
    
    ggplot(nmds_data$df, aes(x = MDS1, y = MDS2)) +
      geom_path(data = nmds_data$df_centroids,
                aes(x = mean_MDS1, y = mean_MDS2, group = gavage, color = gavage),
                arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.5) +
      geom_point(aes(x = MDS1, y = MDS2, color = gavage, size = day_factor), 
                 shape = 21, stroke = 0.6, fill = NA, alpha = 1) +
      geom_point(aes(x = MDS1, y = MDS2, fill = gavage, size = day_factor), 
                 shape = 21, stroke = 0, alpha = 0.35) +
      stat_ellipse(aes(color = gavage, group = gavage), type = "norm", level = 0.95, linewidth = 0.5) +
      # ellipses assume multivariate normality and are thus only provided for visualization
      scale_color_manual(values = gavage_colors, drop = TRUE, name = "Gavage") +
      scale_fill_manual(values = gavage_colors, drop = TRUE, name = "Gavage") +
      scale_size_manual(values = c("7" = 1, "14" = 2, "21" = 3, "28" = 4), name = "Day") +
      labs(subtitle = stress_text) +
      xlab("NMDS1") + ylab("NMDS2") + theme_minimal()
  })
  
  # info pop-up for NMDS Plot
  observeEvent(input$info_nmds_plot, {
    showModal(
      modalDialog(
        title = "Non-metric Multidimensional Scaling Plot",
        p("NMDS is an unconstrained ordination method that represents samples in a low-dimensional space by preserving 
          the rank order of pairwise dissimilarities rather than their absolute values."),
        p("Unlike PCoA, NMDS does not rely on eigen-decomposition and does not produce axes that explain a defined 
          proportion of variance. Instead, NMDS uses iterative optimization to minimize stress (i.e., a measure of mismatch 
          between observed dissimilarities and distances in the ordination space)."),
        p("NMDS is particularly well-suited for non-Euclidean distance metrics commonly used in microbiome studies 
          (e.g., Bray–Curtis, Jaccard, Canberra), where preserving rank relationships is often more meaningful 
          than preserving absolute distances (microbiome data is compositional and not absolute)."),
        p(strong("Bray-Curtis distance:"), "Measures the compositional dissimilarity between two samples based on the abundances of taxa present in at least one of the samples. 
          It is computed as the weighted sum of absolute differences, where weights are the abundances, thus the metric is dominated by abundant taxa."),
        p(strong("Canberra distance:"), "Measures compositional dissimilarity between two samples based on the abundances of taxa present in at least one of the samples. 
          It sums the ratios of the absolute differences to the sums of the abundances for each taxon, thereby giving relatively equal weight to rare and abundant taxa 
          (i.e., rare taxa contribute proportionally more than they would in Bray-Curtis)."),
        p(strong("Jaccard distance:"), "Measures the dissimilarity between two samples based on presence/absence of taxa, ignoring their abundances. 
          It quantifies how many taxa are shared versus unique to each sample."),
        easyClose = TRUE,
        footer = modalButton("Close")
      )
    )
  })
  
  
  # reactive data.frame for PCA plot
  beta_pca_df <- reactive({
    metric <- input$beta_metric
    ord <- beta_diversity[[metric]]$ordination$pca
    df <- as.data.frame(ord$x) # extract coordinates
    df$sample_id <- rownames(df)
    
    # filter df to only include sample_ids that are in the filtered metadata
    df <- df %>% filter(sample_id %in% beta_plot_meta()$sample_id)
    
    # merge with filtered metadata
    df <- left_join(df, beta_plot_meta(), by = "sample_id")
    
    # compute variance explained
    var_explained <- round(100 * summary(ord)$importance[2, 1:2], 1) # extract percentage of variance explained
    
    # compute centroids per gavage x day
    df_centroids <- df %>%
      group_by(gavage, day_factor) %>%
      summarize(mean_PC1 = mean(PC1),
                mean_PC2 = mean(PC2),
                .groups = "drop") %>%
      arrange(gavage, as.numeric(as.character(day_factor)))
    
    list(df = df, df_centroids = df_centroids, var_explained = var_explained)
  })
  
  # PCA plot
  output$pca_plot <- renderPlot({
    
    # only allow Euclidean distance metrics to be used for PCA plot
    metric <- input$beta_metric
    ord <- beta_diversity[[metric]]$ordination$pca

    validate(
      need(!is.null(ord),
            "PCA is only available for Euclidean distances (Aitchison). PCA relies on linear geometry and variance–covariance structure, which are only well-defined in Euclidean space. Non-Euclidean distance metrics (e.g., Bray–Curtis, Jaccard and Canberra) do not preserve this structure, so applying PCA would result in axes and variance estimates that are not mathematically meaningful."))
    
    pca_data <- beta_pca_df()
    
    ggplot(pca_data$df, aes(x = PC1, y = PC2)) +
      geom_path(data = pca_data$df_centroids,
                aes(x = mean_PC1, y = mean_PC2, group = gavage, color = gavage),
                arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.5) +
      geom_point(aes(x = PC1, y = PC2, color = gavage, size = day_factor), 
                 shape = 21, stroke = 0.6, fill = NA, alpha = 1) +
      geom_point(aes(x = PC1, y = PC2, fill = gavage, size = day_factor), 
                 shape = 21, stroke = 0, alpha = 0.35) +
      stat_ellipse(aes(color = gavage, group = gavage), type = "norm", level = 0.95, linewidth = 0.5) +
      # ellipses assume multivariate normality and are thus only provided for visualization
      scale_color_manual(values = gavage_colors, drop = TRUE, name = "Gavage") +
      scale_fill_manual(values = gavage_colors, drop = TRUE, name = "Gavage") +
      scale_size_manual(values = c("7" = 1, "14" = 2, "21" = 3, "28" = 4), name = "Day") +
      xlab(paste0("PC1 (", round(pca_data$var_explained[1], 1), "%)")) +
      ylab(paste0("PC2 (", round(pca_data$var_explained[2], 1), "%)")) +
      theme_minimal()
  })
  
  # info pop-up for PCA Plot
  observeEvent(input$info_pca_plot, {
    showModal(
      modalDialog(
        title = "Principal Component Analysis Plot",
        p("PCA is an unconstrained ordination method that represents samples in a low-dimensional space by 
          directly decomposing the multivariate data matrix, rather than a distance or dissimilarity matrix."),
        p("Each axis (principal component) represents the amount of variation in the original data explained by that axis."),
        p("PCA requires Euclidean geometry and is therefore only appropriate when the data lie in Euclidean space."),
        p(strong("Aitchison distance:"), "The Euclidean (straight-line) distance between two samples after centered-log ratio (CLR) transformation. CLR transformation maps 
        compositional microbiome data from the simplex into Euclidean space by removing the constant-sum constraint, making Euclidean distance geometrically meaningful. 
          Distances therefore reflect differences in log-ratios (i.e., relative fold changes) between components."),
        easyClose = TRUE,
        footer = modalButton("Close")
      )
    )
  })
  

  # overall PERMANOVA table
  output$permanova_overall <- renderTable({
    metric <- input$beta_metric
    df <- as.data.frame(beta_diversity[[metric]]$permanova$overall)
    
    # format table
    df_formatted <- df %>% 
      tibble::rownames_to_column("Model Term") %>%
      mutate(Df = format(round(Df, 0), nsmall = 0),
             SumOfSqs = format(round(SumOfSqs, 4), nsmall = 4),
             R2 = format(round(R2, 5), nsmall = 5),
             F = ifelse(is.na(F), "", format(round(F, 3), nsmall = 3)),
             `Pr(>F)` = ifelse(is.na(`Pr(>F)`), "", format(round(`Pr(>F)`, 3), nsmall = 3)))
    
    df_formatted
  })
  
  # info pop-up for overall PERMANOVA table
  observeEvent(input$info_permanova_overall, {
    showModal(
      modalDialog(
        title = "How to Read the PERMANOVA Table",
        p("PERMANOVA tests whether the centroids of groups in multivariate space are significantly farther apart than expected by chance. 
          Permutations are used to calculate p-values as the distribution of distances is usually non-normal."),
        p(strong("Df (Degrees of freedom):"), "Number of independent observations in a dataset that are available to estimate parameters or variability.
          In repeated-measures designs or stratified analysis, Df is reduced since observations are not fully independent."),
        p(strong("SumOfSqs (Sum of squares):"), "Total variance explained by the model (model) and not explained by the model (residual)."),
        p(strong("R2:"), "Proportion of variance explained."),
        p(strong("F:"), "F-statistic from permutation testing (ratio of signal to noise)."),
        p(strong("Pr(>F):"), "Permutation-based p-value."),
        easyClose = TRUE,
        footer = modalButton("Close")
      )
    )
  })
  
  
  # time-resolved PERMANOVA table
  output$permanova_time <- renderDT({
    metric <- input$beta_metric
    df <- as.data.frame(beta_diversity[[metric]]$permanova$by_time)
    
    # format table
    df_formatted <- df %>%
      mutate(SumsOfSqs = format(round(SumsOfSqs, 4), nsmall = 4),
             F.Model = format(round(F.Model, 4), nsmall = 4),
             R2 = format(round(R2, 4), nsmall = 4),
             p.value = format(round(p.value, 4), nsmall = 4),
             p.adjusted = format(round(p.adjusted, 4), nsmall = 4),
             p.adjusted_global = format(round(p.adjusted_global, 4), nsmall = 4))
    
    df_formatted
  })
  
  # info pop-up for time-resolved PERMANOVA table
  observeEvent(input$info_permanova_time, {
    showModal(
      modalDialog(
        title = "How to Read the Time-resolved PERMANOVA Table",
        p("PERMANOVA tests whether the centroids of groups are significantly farther apart in the multivariate space than expected by chance. 
          Permutations are used to calculate p-values as the distribution of distances is usually non-normal."),
        p(strong("pairs:"), "The two treatment groups being compared."),
        p(strong("SumsOfSqs (Sum of squares):"), "Total variance explained by the difference between the two groups"),
        p(strong("F.model:"), "F-statistic from permutation testing (ratio of signal to noise)."),
        p(strong("R2:"), "Proportion of variation, at that timepoint, explained by the group."),
        p(strong("p.value:"), "Raw p-value from the permutation test."),
        p(strong("p.adjusted:"), "The p-value adjusted for multiple testing within the set of pairwise comparisons at the indicated timepoint."),
        p(strong("Time:"), "The timepoint at which the comparison was made."),
        p(strong("p.adjusted_global:"), "The p-value adjusted for multiple testing across all time points and comparisons."),
        easyClose = TRUE,
        footer = modalButton("Close")
      )
    )
  })
  
  
  # overall dispersion
  output$dispersion_overall <- renderTable({
    metric <- input$beta_metric
    df <- as.data.frame(beta_diversity[[metric]]$dispersion$overall$tab)
    
    # format table
    df_formatted <- df %>%
      tibble::rownames_to_column("Term") %>%
      mutate(Df = format(round(Df, 0), nsmall = 0),
             `Sum Sq` = format(round(`Sum Sq`, 6), nsmall = 6),
             `Mean Sq` = format(round(`Mean Sq`, 7), nsmall = 7),
             F = ifelse(is.na(F), "", format(round(F, 4), nsmall = 4)),
             N.Perm = ifelse(is.na(N.Perm), "", format(round(N.Perm, 0), nsmall = 0)),
             `Pr(>F)` = ifelse(is.na(`Pr(>F)`), "", format(round(`Pr(>F)`, 3), nsmall = 3)))
    
    df_formatted
  })
  
  # info pop-up for overall dispersion table
  observeEvent(input$info_dispersion_overall, {
    showModal(
      modalDialog(
        title = "How to Read the Multivariate Dispersion Table",
        p("PERMANOVA assumes homogeneity of dispersion. If dispersions are significantly different, significant PERMANOVA results could be 
          due to differences in variance and not differences in group centroids."),
        p(strong("Df (Degrees of freedom):"), "Number of independent observations in a dataset that are available to estimate parameters or variability.
          In repeated-measures designs or stratified analysis, Df is reduced since observations are not fully independent."),
        p(strong("Sum Sq (Sum of squares):"), "Total variance explained by the group differences in dispersion (Groups) and not explained by the group differences in dispersion (Residuals)."),
        p(strong("Mean Sq (Sum Sq divided by Df):"), "Average squared distance per degree of freedom"),
        p(strong("F:"), "F-statistic from permutation testing (ratio of signal to noise)."),
        p(strong("N.Perm:"), "Number of permutations used to calculate the p-value."),
        p(strong("Pr(>F):"), "Permutation-based p-value."),
        easyClose = TRUE,
        footer = modalButton("Close")
      )
    )
  })
  
  
  # time-resolved dispersion
  output$dispersion_time <- renderTable({
    metric <- input$beta_metric
    
    # extract and format dispersion table from each timepoint
    df_list <- lapply(beta_diversity[[metric]]$dispersion$by_time, function(x) {
      x$tab %>%
        
        # format table
        as.data.frame() %>%
        tibble::rownames_to_column("Term") %>%
        mutate(Df = format(round(Df, 0), nsmall = 0),
               `Sum Sq` = format(round(`Sum Sq`, 6), nsmall = 6),
               `Mean Sq` = format(round(`Mean Sq`, 7), nsmall = 7),
               F = ifelse(is.na(F), "", format(round(F, 4), nsmall = 4)),
               N.Perm = ifelse(is.na(N.Perm), "", format(round(N.Perm, 0), nsmall = 0)),
               `Pr(>F)` = ifelse(is.na(`Pr(>F)`), "", format(round(`Pr(>F)`, 3), nsmall = 3)))
    })
    
    # add a Time column for timepoints
    time_points <- sort(unique(metadata$day))
    stopifnot(length(time_points) == length(df_list))
    
    for (i in seq_along(df_list)) {
      df_list[[i]]$Time <- time_points[i]
    }
    
    combined <- bind_rows(df_list)
    combined

  })
  
  # info pop-up for time-resolved dispersion table
  observeEvent(input$info_dispersion_time, {
    showModal(
      modalDialog(
        title = "How to Read the Time-resolved Multivariate Dispersion Table",
        p("PERMANOVA assumes homogeneity of dispersion. If dispersions are significantly different, significant PERMANOVA results could be 
          due to differences in variance and not differences in group centroids."),
        p(strong("Df (Degrees of freedom):"), "Number of independent observations in a dataset that are available to estimate parameters or variability.
          In repeated-measures designs or stratified analysis, Df is reduced since observations are not fully independent."),
        p(strong("Sum Sq (Sum of squares):"), "Total variance explained by the group differences in dispersion (Groups) and not explained by the group differences in dispersion (Residuals)."),
        p(strong("Mean Sq (Sum Sq divided by Df):"), "Average squared distance per degree of freedom"),
        p(strong("F:"), "F-statistic from permutation testing (ratio of signal to noise)."),
        p(strong("N.Perm:"), "Number of permutations used to calculate the p-value."),
        p(strong("Pr(>F):"), "Permutation-based p-value."),
        p(strong("Time"), "The timepoint at which the comparison was made."),
        easyClose = TRUE,
        footer = modalButton("Close")
      )
    )
  })
  
  
  # reactive data.frame for constrained ordination plot
  beta_con_ord_df <- reactive({
    metric <- input$beta_metric
    con <- beta_diversity[[metric]]$ordination$con_ord
    cap_df <- as.data.frame(scores(con, display = "sites")) # get CAP coordinates
    cap_df$sample_id <- rownames(cap_df)

    # filter cap_df to only include sample_ids that are in the filtered metadata
    cap_df <- cap_df %>% filter(sample_id %in% beta_plot_meta()$sample_id)
    
    # merge with filtered metadata
    cap_df <- left_join(cap_df, beta_plot_meta(), by = "sample_id")
    
    # compute variance explained
    var_explained <- con$CCA$eig / sum(con$tot.chi) * 100
    
    # compute centroids
    axis_names <- colnames(cap_df)[1:2] # column names depend on method (capscale or rda)
    
    cap_df <- cap_df %>%
      rename(x = !!axis_names[1], 
             y = !!axis_names[2]) 
    
    df_centroids <- cap_df %>%
      group_by(gavage, day_factor) %>%
      summarize(mean_x = mean(x),
                mean_y = mean(y),
                .groups = "drop") %>%
      arrange(gavage, as.numeric(as.character(day_factor)))
    
    list(df = cap_df, df_centroids = df_centroids, var_explained = var_explained)
  })
  
  # constrained ordination plot
  output$con_ord_plot <- renderPlot({
    con_ord_data <- beta_con_ord_df()

    ggplot(con_ord_data$df, aes(x = x, y = y)) +
      geom_path(data = con_ord_data$df_centroids,
                aes(x = mean_x, y = mean_y, group = gavage, color = gavage),
                arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.5) +
      geom_point(aes(color = gavage, size = day_factor),
                 shape = 21, stroke = 0.6, fill = NA, alpha = 1) +
      geom_point(aes(fill = gavage, size = day_factor), 
                 shape = 21, stroke = 0, alpha = 0.35) +
      stat_ellipse(aes(color = gavage, group = gavage), type = "norm", level = 0.95, linewidth = 0.5) +
      # ellipses assume multivariate normality and are thus only provided for visualization
      scale_color_manual(values = gavage_colors, drop = TRUE, name = "Gavage") +
      scale_fill_manual(values = gavage_colors, drop = TRUE, name = "Gavage") +
      scale_size_manual(values = c("7" = 1, "14" = 2, "21" = 3, "28" = 4), name = "Day") +
      xlab(paste0("CAP1 (", round(con_ord_data$var_explained[1], 1), "%)")) +
      ylab(paste0("CAP2 (", round(con_ord_data$var_explained[2], 1), "%)")) +
      theme_minimal()
  })
  
  # info pop-up for constrained ordination plot
  observeEvent(input$info_con_ord_plot, {
    showModal(
      modalDialog(
        title = "Constrained Ordination Plot",
        p("Constrained ordination represents samples in a low-dimensional space while focusing on variation explained by one or more experimental variables (e.g., gavage and day)."),
        p("Unlike unconstrained methods (e.g., PCoA, NMDS or PCA), the ordination axes are restricted to linear combinations of the experimental variables (the percentage of variance on each axis reflects 
          the proportion of variation explained by the experimental variables). Thus constrained ordination is particularly useful for assessing how much variation in microbiome composition can be explained by these variables"), 
        p("Constrained ordiantion can be used with both Euclidean and non-Euclidean distance metrics."),
        hr(),
        p(strong("Bray-Curtis distance:"), "Measures the compositional dissimilarity between two samples based on the abundances of taxa present in at least one of the samples. 
          It is computed as the weighted sum of absolute differences, where weights are the abundances, thus the metric is dominated by abundant taxa."),
        p(strong("Canberra distance:"), "Measures compositional dissimilarity between two samples based on the abundances of taxa present in at least one of the samples. 
          It sums the ratios of the absolute differences to the sums of the abundances for each taxon, thereby giving relatively equal weight to rare and abundant taxa 
          (i.e., rare taxa contribute proportionally more than they would in Bray-Curtis)."),
        p(strong("Jaccard distance:"), "Measures the dissimilarity between two samples based on presence/absence of taxa, ignoring their abundances. 
          It quantifies how many taxa are shared versus unique to each sample."),
        p(strong("Aitchison distance:"), "The Euclidean (straight-line) distance between two samples after centered log ratio (CLR) transformation. CLR transformation maps 
        compositional microbiome data from the simplex into Euclidean space by removing the constant-sum constraint, making Euclidean distance geometrically meaningful."),
        easyClose = TRUE,
        footer = modalButton("Close")
      )
    )
  })
  

  # constrained ordination - overall permutation test
  output$con_ord_overall <- renderTable({
    metric <- input$beta_metric
    df <- as.data.frame(beta_diversity[[metric]]$anova_perm$overall)
    
    # format table
    df_formatted <- df %>% 
      tibble::rownames_to_column("Model Term") %>%
      mutate(Df = format(round(Df, 0), nsmall = 0),
             F = ifelse(is.na(F), "", format(round(F, 3), nsmall = 3)),
             `Pr(>F)` = ifelse(is.na(`Pr(>F)`), "", format(round(`Pr(>F)`, 3), nsmall = 3)))
    
    # ANOVA with capscale has SumOfSqs and ANOVA with rda has Variance
    if ("SumOfSqs" %in% colnames(df_formatted)) {
      df_formatted <- df_formatted %>%
        mutate(SumOfSqs = format(round(SumOfSqs, 4), nsmall = 4))
    } else if ("Variance" %in% colnames(df_formatted)) {
      df_formatted <- df_formatted %>%
        mutate(Variance = format(round(Variance, 4), nsmall = 4))
    }
    
    df_formatted
  })
  
  # constrained ordination - permutation tests by terms
  output$con_ord_terms <- renderTable({
    metric <- input$beta_metric
    df <- as.data.frame(beta_diversity[[metric]]$anova_perm$terms) 
    
    # format table
    df_formatted <- df %>% 
      tibble::rownames_to_column("Model Term") %>%
      mutate(Df = format(round(Df, 0), nsmall = 0),
             F = ifelse(is.na(F), "", format(round(F, 3), nsmall = 3)),
             `Pr(>F)` = ifelse(is.na(`Pr(>F)`), "", format(round(`Pr(>F)`, 3), nsmall = 3)))
    
    # ANOVA with capscale has SumOfSqs and ANOVA with rda has Variance
    if ("SumOfSqs" %in% colnames(df_formatted)) {
      df_formatted <- df_formatted %>%
        mutate(SumOfSqs = format(round(SumOfSqs, 4), nsmall = 4))
    } else if ("Variance" %in% colnames(df_formatted)) {
      df_formatted <- df_formatted %>%
        mutate(Variance = format(round(Variance, 4), nsmall = 4))
    }
    
    df_formatted
  })
  
  # constrained ordination - permutation tests by axis
  output$con_ord_axis <- renderTable({
    metric <- input$beta_metric
    df <- as.data.frame(beta_diversity[[metric]]$anova_perm$axis)
    
    # format table
    df_formatted <- df %>% 
      tibble::rownames_to_column("Model Term") %>%
      mutate(Df = format(round(Df, 0), nsmall = 0),
             F = ifelse(is.na(F), "", format(round(F, 3), nsmall = 3)),
             `Pr(>F)` = ifelse(is.na(`Pr(>F)`), "", format(round(`Pr(>F)`, 3), nsmall = 3)))
    
    # ANOVA with capscale has SumOfSqs and ANOVA with rda has Variance
    if ("SumOfSqs" %in% colnames(df_formatted)) {
      df_formatted <- df_formatted %>%
        mutate(SumOfSqs = format(round(SumOfSqs, 4), nsmall = 4))
    } else if ("Variance" %in% colnames(df_formatted)) {
      df_formatted <- df_formatted %>%
        mutate(Variance = format(round(Variance, 4), nsmall = 4))
    }
    
    df_formatted
  })
  
  # info pop-up for constrained ordination permutation test tables
  observeEvent(input$info_con_ord_overall, {
    showModal(
      modalDialog(
        title = "How to Read the Constrained Ordination ANOVA Tables",
        p("A constrained ordination ANOVA tests how much of the variation in community composition is explained by the constrained variables and whether 
        this explained variation is greater than expected by chance. Permutations are used to calculate p-values as the distribution of distances is usually non-normal."),
        p(strong("Overall Permutation Test:"), "Tests whether the total variation in community structure explained jointly by gavage and day is greater than expected by chance."), 
        p(strong("Permutation Test by Term:"), "Tests whether each term (gavage and day) explains more variation than expected by chance after accounting for the other term(s)."), 
        p(strong("Permutation Test by Axis:"), "Tests whether each constrained ordination axis explains more of the variation attributable to the constraining variables than expected by chance."), 
        p(strong("Df (Degrees of freedom):"), "Number of independent observations in a dataset that are available to estimate parameters or variability.
          In repeated-measures designs or stratified analysis, Df is reduced since observations are not fully independent."),
        p(strong("SumOfSqs (Sum of squares)/Variance:"), "Total variance explained by the model (model) and not explained by the model (residual). Sum of squares (sum of squared distances to a generalized centroid) 
        is used with non-Euclidean distance matrices to approximate variance, since variance is not well-defined in non-Euclidean space (distances don't obey Euclidean rules/Pythagoras’ theorem)."),
        p(strong("F:"), "F-statistic from permutation testing (ratio of signal to noise)."),
        p(strong("Pr(>F):"), "Permutation-based p-value."),
        easyClose = TRUE,
        footer = modalButton("Close")
      )
    )
  })
  
  
  ##############################################################
  ##########      DIFFERENTIAL ABUNDANCE ANALYSIS     ##########
  ##############################################################
  
  # filter for top N taxa or manually selected taxa
  selected_taxa <- reactive({
    
    df <- differential_abundance$comm_over
    
    if (input$comm_taxa_method == "top_n") {
      
      # top N most abundant taxa 
      df_taxa <- df %>%
        group_by(OTU) %>%
        summarize(mean_abundance = mean(abundance), .groups = "drop") %>%
        arrange(desc(mean_abundance)) %>%
        slice_head(n = input$top_n_taxa) %>%
        pull(OTU)
      
    } else if (input$comm_taxa_method == "manual") {
      
      # manually selected taxa
      df_taxa <- input$comm_taxa_select
    }
    
    df_taxa
  })
  
  
  # format matrix and plot heatmap
  output$rel_abun_heatmap <- renderPlot({
    
    heat_mat <- differential_abundance$comm_over %>%
      filter(gavage %in% input$gavage_plot) %>%
      filter(OTU %in% selected_taxa()) %>%
      group_by(gavage_day, OTU) %>%
      summarise(abundance = mean(abundance), .groups = "drop") %>% 
      pivot_wider(names_from = gavage_day, values_from = abundance, values_fill = 0) %>%
      column_to_rownames("OTU") %>%
      as.matrix()
    
    pheatmap(heat_mat, scale = "row", cluster_cols = FALSE,
             clustering_distance_rows = "euclidean",
             clustering_method = "complete")
    
  })
  
  
  # data.frame for stacked barplot
  bar_df <- reactive({
    df <- differential_abundance$comm_over %>%
      filter(gavage %in% input$gavage_plot) %>%
      mutate(plot_taxa = ifelse(OTU %in% selected_taxa(), OTU, "Other")) 
    
    # set factor levels
    df$plot_taxa <- factor(df$plot_taxa, levels = c(selected_taxa(), "Other"))
    
    # group and summarise by abundance
    df %>%
      group_by(gavage, day, day_factor, gavage_day, plot_taxa) %>%
      summarise(abundance = sum(abundance), .groups = "drop")
  })
  
  # reactive colors for stacked barplot
  base_colors <- c("#862185", "#009E73", "#88CCEE", "#CC6677", "#D55E00", 
                   "#44AA99", "#332288", "#E69F00", "#0072B2", "#AA4499",
                   "#F0E442", "#117733", "#DB72FB", "#619CFF", "#882255")
  
  fill_colors <- reactive({
    taxa_levels <- levels(bar_df()$plot_taxa)
    
    # all levels except Other
    selected_levels <- taxa_levels[taxa_levels != "Other"]
    
    # assign base_colors sequentially (cycle if more taxa than colors)
    cols <- rep(base_colors, length.out = length(selected_levels))
    names(cols) <- selected_levels
    
    # Other is gray
    c(cols, Other = "#999999")
  })
  
  # plot stacked barplot
  output$rel_abun_barplot <- renderPlot({
    ggplot(bar_df(), aes(x = gavage_day, y = abundance, fill = plot_taxa)) +
      geom_col() + theme_minimal() +
      scale_fill_manual(values = fill_colors(), name = "Taxa") +
      labs(x = "Gavage Group and Day", y = "Mean Relative Abundance") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  })
  
 
  ### relative abundance over time
  
  # filter for selected gavage groups for relative abundance plot
  rel_abun_plot_df <- reactive({
    differential_abundance$rel_abun %>%
      filter(gavage %in% input$gavage_plot) # gavage already in differential_abundance$rel_abun
  })
  
  # relative abundance over time plot
  output$gamm_indiv_taxon_rel_abun_plot <- renderPlot({
    
    plot_df <- rel_abun_plot_df() %>% filter(taxon == input$taxon_select)
    ggplot(plot_df, aes(x = day, y = mean_abundance, color = gavage, group = gavage, fill = gavage)) +
      geom_line() + geom_point() + theme_minimal() +
      geom_ribbon(aes(ymin = mean_abundance - sem, ymax = mean_abundance + sem), alpha = 0.2, color = NA) +
      scale_color_manual(values = gavage_colors, drop = FALSE) +
      scale_fill_manual(values = gavage_colors, drop = FALSE) +
      labs(x = "Day", y = "Mean Relative Abundance", color = "Gavage", fill = "Gavage")
  })
  
 
  # info pop-up for relative abundance over time plot
  observeEvent(input$info_gamm_indiv_taxon_rel_abun_plot, {
    showModal(
      modalDialog(
        title = "Relative Abundance of Individual Taxon Over Time",
        p("Plot of the mean relative abundance of the selected taxon over time for each gavage group. 
          Shaded ribbons represent the standard error of the mean (SEM)."),
        p(strong(":"), "."),
        footer = modalButton("Close")
      )
    )
  })
 
  
  ### generalized additive mixed model
  
  # filter for selected gavage groups for differential abundance smooth plot
  diff_abun_plot_df <- reactive({
    differential_abundance$gamm$data %>%
      filter(gavage %in% input$gavage_plot) # gavage already in differential_abundance$gamm$data
  })
  
  # GAMM smooth plots
  output$gamm_indiv_taxon_smooth_plot <- renderPlot({
    
    # check if model fit without error
    model_obj <- differential_abundance$gamm$models[[input$taxon_select]]
    
    validate(
      need(model_obj$success,
           paste0("GAMM failed to fit for ", input$taxon_select)))
    
    plot_df <- diff_abun_plot_df() %>% filter(taxon == input$taxon_select)
    smooth_df <- differential_abundance$gamm$models[[input$taxon_select]]$smooth_pred
    ggplot() + theme_minimal() +
      geom_point(data = plot_df, aes(x = day, y = abundance, color = gavage), alpha = 0.4) +
      geom_line(data = smooth_df, aes(x = day, y = fit, color = gavage), linewidth = 1) +
      geom_ribbon(data = smooth_df, aes(x = day, ymin = lower, ymax = upper, fill = gavage), alpha = 0.2, color = NA) +
      scale_color_manual(values = gavage_colors, drop = FALSE) +
      scale_fill_manual(values = gavage_colors, drop = FALSE) +
      labs(y = input$taxon_select, x = "Day", color = "Gavage", fill = "Gavage")
  })
  
  # info pop-up for smooths plot
  observeEvent(input$info_gamm_indiv_taxon_smooth_plot, {
    showModal(
      modalDialog(
        title = "GAMM-estimated Smooths Over Time",
        p("Plot of the CLR-transformed abundance of the selected taxon over time for each sample. Lines represent the estimated smooth 
        trajectories from the GAMM (the temporal trend for each group on the CLR scale). Shaded ribbons represent 95% confidence intervals."),
        p("Generalized Additive Mixed Models model both fixed effects (e.g., time and treatment) and random effects (e.g., subject id). 
          GAMMs allow non-linear relationships between predictors and the response via smooth functions (splines)."),
        p(strong("GAMM-estimated smooth trajectories:"), "Estimate the relationship between a predictor and a response without assuming a fixed parametric form. 
          Smooths are penalized and therefore only allow curvature if the data strongly supports it."),
        footer = modalButton("Close")
      )
    )
  })

  
  # gamm model summary - parametric coefficients
  output$gamm_model_para_summary <- renderTable({
    
    # check if model fit without error
    model_obj <- differential_abundance$gamm$models[[input$taxon_select]]
    
    validate(
      need(model_obj$success,
           paste0("GAMM failed to fit for ", input$taxon_select)))
    
    df <- model_obj$summary_param
    
    # format table
    df_formatted <- df %>% 
      mutate(estimate = format(round(estimate, 4), nsmall = 4),
             std_error = format(round(std_error, 4), nsmall = 4),
             t_value = format(round(t_value, 4), nsmall = 4),
             p_value = formatC(p_value, format = "e", digits = 4),
             p_adj = formatC(p_adj, format = "e", digits = 4)) %>%
      select(term, estimate, std_error, t_value, p_value, p_adj)
    
    df_formatted
  })
  
  # gamm model summary - smooths
  output$gamm_model_smooth_summary <- renderTable({
    
    # check if model fit without error
    model_obj <- differential_abundance$gamm$models[[input$taxon_select]]
    
    validate(
      need(model_obj$success,
           paste0("GAMM failed to fit for ", input$taxon_select)))
    
    df <- model_obj$summary_smooth
    
    # format table
    df_formatted <- df %>% 
      mutate(edf = format(round(edf, 4), nsmall = 4),
             ref_df = format(round(ref_df, 4), nsmall = 4),
             F_value = format(round(F_value, 4), nsmall = 4),
             p_value = formatC(p_value, format = "e", digits = 4),
             p_adj = formatC(p_adj, format = "e", digits = 4)) %>%
      select(smooth, edf, ref_df, F_value, p_value, p_adj)
    
    df_formatted
  })
  
  # info pop-up for the GAMM summary (parametric coefficients and smooths)
  observeEvent(input$info_gamm_indiv_taxon_model_summary, {
    showModal(
      modalDialog(
        title = "How to Interpret the GAMM Summary",
        p(strong("Parametric coefficients:"), "Represent the linear (fixed) effects in the model. Each group term is the linear difference from the reference group averaged over time."),
        p(strong("estimate"), "Estimated effect size (difference in CLR-transformed abundance relative to the reference group)."),
        p(strong("std_error"), "Uncertainty of the estimated effect size."),
        p(strong("t_value"), "Test statistic for the parametric effect (estimate divided by the standard error)."),
        hr(),
        p(strong("Smooth terms:"), "Indicates whether adding nonlinearity significantly improves the model for a given group."),
        p(strong("edf (effective degrees of freedom:"), "How wiggly the smooth is. > 1 = some nonlinearity"),
        p(strong("ref_df"), "Reference degrees of freedom used for significance testing of the smooth."),
        p(strong("F_value:"), "Test statistic assessing whether the smooth explains significant variation over time."),
        hr(),
        p(strong("p_value/p_adj"), "Significance before and after mutliple-testing correction (BH)."),
        footer = modalButton("Close")
      )
    )
  })
  
  
  # filter for taxon for GAMM smooth differences plots
  smooth_diff_indiv_taxon_plot_df <- reactive({
    differential_abundance$gamm$smooth_diff %>%
      filter(taxon == input$taxon_select)
  })
  
  # plot GAMM pairwise smooth differences
  output$gamm_indiv_taxon_smooth_diff_plot <- renderPlot({
    draw(smooth_diff_indiv_taxon_plot_df()) &
      scale_x_continuous(name = "Day")
  })
  
  # info pop-up for the GAMM smooth differences plots
  observeEvent(input$info_gamm_indiv_taxon_smooth_diff_plot, {
    showModal(
      modalDialog(
        title = "How to Interpret the GAMM Smooth Differences Plots",
        p("These plots show the pairwise differences between the estimated smooth trajectories of the selected taxon for each gavage group. 
        (i.e., whether and when one group differs significantly from another over time on the CLR-transformed scale.)"),
        p(strong("Positive values:"), "The first group in the comparison has higher CLR abundance than the second group at that timepoint."),
        p(strong("Negative values:"), "The first group in the comparison has lower CLR abundance than the second group at that timepoint."),
        p(strong("Shaded regions:"), "Indicate 95% confidence intervals for the differences. If the shaded region does not overlap zero, then 
          the difference is statistically significant at that timepoint."),
        footer = modalButton("Close")
      )
    )
  })
  
  
  # plot GAMM residuals plot
  output$gamm_indiv_taxon_resid_plot <- renderPlot({
    
    # check if model fit without error
    model_obj <- differential_abundance$gamm$models[[input$taxon_select]]
    
    validate(
      need(model_obj$success,
           paste0("GAMM failed to fit for ", input$taxon_select)))
    
    model <- model_obj$model$gam
    plot(fitted(model), resid(model),
         ylab = "Residuals", xlab = "Fitted Values",
         main = paste0("GAMM Residuals vs Fitted Values Plot: ", input$taxon_select))
    abline(h = 0, lty = 2)
  })
  
  # plot GAMM normal Q-Q plot
  output$gamm_indiv_taxon_qq_plot <- renderPlot({
    
    # check if model fit without error
    model_obj <- differential_abundance$gamm$models[[input$taxon_select]]
    
    validate(
      need(model_obj$success,
           paste0("GAMM failed to fit for ", input$taxon_select)))
    
    res <- model_obj$residuals
    qqnorm(res, 
           ylab = "Deviance Residuals", xlab = "Theoretical Quantiles",
           main = paste0("GAMM Normal Q-Q Plot: ", input$taxon_select))
    qqline(res)
  })
  
  
  # plot GAMM response plot
  output$gamm_indiv_taxon_resp_plot <- renderPlot({
    
    # check if model fit without error
    model_obj <- differential_abundance$gamm$models[[input$taxon_select]]
    
    validate(
      need(model_obj$success,
           paste0("GAMM failed to fit for ", input$taxon_select)))
    
    model <- model_obj$model$gam
    plot_df <- diff_abun_plot_df() %>% filter(taxon == input$taxon_select)
    plot(fitted(model), plot_df$abundance,
         ylab = "CLR Abundance", xlab = "Fitted Values",
         main = paste0("GAMM Response vs Fitted Values Plot: ", input$taxon_select))
    abline(a = 0, b = 1, lty = 2)
  })
  
  # info pop-up for GAM residuals vs fitted plot, normal Q-Q plot and response vs fitted plot
  observeEvent(input$info_gamm_indiv_taxon_resid_plot, {
    showModal(
      modalDialog(
        title = "How to Interpret the Residuals vs Fitted Values Plot, the Normal Q-Q Plot and the Response vs Fitted Values Plot",
        p("The Residuals versus Fitted Values plot is used to check whether the statistical model is appropriate for the data by checking whether the remaining errors (residuals) are random 
          or if they show structure the model failed to capture. The points should be randomly scattered around zero (i.e., no clear pattern or trend)."),
        p(strong("Curved pattern:"), "Model failed to capture nonlinear strucutre in the mean (e.g., missing interaction)."),
        p(strong("Curved clusters or bands:"), "Group-specific trends or correlation structure is not adequately modeled."),
        p(strong("Widening/narrowing:"), "Non-constant variance (heteroscedasticity)."),
        p(strong("Extreme outliers:"), "Influential points or possible errors."),
        hr(),
        p("The Normal Q-Q plot (quantile-quantile plot) is used to compare the distribution of the residuals with a theoretical distribution (often normal). 
          Each point represents how one residual compares to what would be expected under normality. LMMs and Gaussian GAMMs assume that residuals are normally distributed. 
          If residuals are not normal, p-values and confidence intervals may be unreliable."),
        p(strong("Points lie on a straight 45° line:"), "Residuals are approximately normally distributed."),
        p(strong("S-shaped curve:"), "Less extreme values (concave up then down) or more extreme values (concave down then up) than normal."),
        p(strong("Points curve away at ends:"), "Indicates skewed residuals."),
        p(strong("Large deviations at the ends:"), "Potential outliers."),
        hr(),
        p("The Response versus Fitted Values plot compares the observed response values to the values predicted by the model. This plot helps you see how well the model 
          captures the overall trends in the data and whether there are systematic deviations."),
        p(strong("Points lie roughly along the y = x line:"), "The model fits the data well; predicted values match observed values."),
        p(strong("Points systematically above or below the line:"), "The model under- or over-predicts in certain ranges."),
        p(strong("Nonlinear patterns or curves:"), "The model may be missing key nonlinear effects or interactions."),
        p(strong("Large scatter or wide spread around the line:"), "High residual variance; model may not explain much of the variation."),
        p(strong("Clusters or gaps:"), "Indicates group-specific effects or unmodeled structure in the data."),
        footer = modalButton("Close")
      )
    )
  })
  
  
  # GAMM concurvity
  output$gamm_indiv_taxon_concurvity <- renderTable({
    
    # check if model fit without error
    model_obj <- differential_abundance$gamm$models[[input$taxon_select]]
    
    validate(
      need(model_obj$success,
           paste0("GAMM failed to fit for ", input$taxon_select)))
    
    df <- as.data.frame(model_obj$concurvity)
    
    # format table
    df_formatted <- df %>% 
      tibble::rownames_to_column("Measure") %>%
      mutate(para = format(round(para, 2), nsmall = 2),
             `s(day_c):gavageG_1DMD` = formatC(`s(day_c):gavageG_1DMD`, format = "e", digits = 4),
             `s(day_c):gavageG_4DMD` = formatC(`s(day_c):gavageG_4DMD`, format = "e", digits = 4),
             `s(day_c):gavageG_4W7C` = formatC(`s(day_c):gavageG_4W7C`, format = "e", digits = 4),
             `s(day_c):gavageG_4WMD` = formatC(`s(day_c):gavageG_4WMD`, format = "e", digits = 4))
    
    df_formatted
  })
  
  # info pop-up for concurvity table
  observeEvent(input$info_gamm_indiv_taxon_concurvity, {
    showModal(
      modalDialog(
        title = "How to Read the Concurvity Table",
        p("How much redundancy or linear dependence exists in the model (i.e., how much one term can explain another term). 
        If two terms are highly dependendent, it is hard for the model to separate their effects. 
        0 = no redundancy, > 0.5 = moderate redundancy, and > 0.9 = very high redundancy."),
        p(strong("Measure:"), "worst = maximum concurvity for the term across all other terms, observed = concurvity calculated 
          directly from the fitted model, and estimate = smoothed/adjusted estimate of concurvity."),
        p(strong("para:"), "The parametric coefficients (see Model Summary (GAMM))."),
        p(strong("s(day):gavage...:"), "The smooth terms/group-specific nonlinear trends over time (see Model Summary (GAMM))."),
        easyClose = TRUE,
        footer = modalButton("Close")
      )
    )
  })
  
}

# launch Shiny app
shinyApp(ui, server)

