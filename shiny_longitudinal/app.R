### Longitudinal Analysis of Gavage Data - Shiny App

## Startup code
# load libraries
library(shiny)
library(bslib)
library(tidyverse)
library(ggplot2)
library(shinyWidgets)
library(DT)

# setwd
setwd("/Users/kristinvandenham/kmvanden/RStudio/microbiome_analysis/shiny_longitudinal")

# load precomputed RDS files
metadata <- readRDS("metadata.rds")
alpha_diversity <- readRDS("alpha_diversity.rds")
beta_diversity <- readRDS("beta_diversity.rds")


# types of diversity metrics
alpha_diversity_metrics <- names(alpha_diversity$models)
beta_diversity_metrics <- names(beta_diversity)

## UI
ui <- fluidPage(theme = bs_theme(version = 5, bootswatch = "flatly", primary = "#2C3E50"), # responsive Bootstrap 5 layout
                titlePanel("Longitudinal Microbiome Diversity Analysis"),
                tags$style(HTML(".tab-content-spacing {margin-top: 20px;}")),
                sidebarLayout(
                  sidebarPanel(pickerInput(inputId = "diversity_type", label = "Diversity Type",
                                           choices = c("alpha", "beta"), selected = "Alpha" ),
                               
                               ### Alpha diversity dropdowns
                               conditionalPanel(condition = "input.diversity_type == 'alpha'",
                                                pickerInput(inputId = "alpha_metric", label = "Alpha Diversity Metric",
                                                             choices = alpha_diversity_metrics, selected = "Observed"), # alpha diversity metric dropdown
                                                 pickerInput(inputId = "gavage", label = "Treatment Group",
                                                             choices = unique(metadata$gavage),
                                                             selected = unique(metadata$gavage),
                                                             multiple = TRUE), # treatment (gavage) group multi-select dropdown
                                                 pickerInput(inputId = "model_type", label = "Model Type",
                                                             choices = c("LMM", "GAMM"), selected = "LMM")),
                               
                               ### Beta diversity dropdowns
                               conditionalPanel(condition = "input.diversity_type == 'beta'",
                                                pickerInput(inputId = "beta_metric", label = "Beta Diversity Metric",
                                                            choices = beta_diversity_metrics, selected = "bray_curtis"),
                                                pickerInput(inputId = "gavage", label = "Treatment Group",
                                                            choices = unique(metadata$gavage),
                                                            selected = unique(metadata$gavage),
                                                            multiple = TRUE))),
                              
                  mainPanel(
                    
                    ### Alpha diversity tabs
                    
                    # LMM tabs
                    conditionalPanel(condition = "input.diversity_type == 'alpha' && input.model_type == 'LMM'",
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
                                                 tabPanel("Residuals and QQ Plots (LMM)",
                                                          div(class = "tab-content-spacing",
                                                          h4("Residuals and QQ Plots (LMM)", 
                                                          actionButton("info_lmm_resid_plot", label = NULL, icon = icon("question-circle"), 
                                                                       style = "font-size:12px;", class ="btn btn-info btn-sm")),
                                                          plotOutput("lmm_resid_plot", height = 400), plotOutput("lmm_qq_plot", height = 400))))),
                    
                    # GAMM tabs
                    conditionalPanel(condition = "input.diversity_type == 'alpha' && input.model_type == 'GAMM'",
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
                                                 tabPanel("Residuals and QQ Plots (GAMM)", 
                                                          div(class = "tab-content-spacing",
                                                          h4("Residuals and QQ Plots (GAMM)", 
                                                          actionButton("info_gamm_resid_plot", label = NULL, icon = icon("question-circle"), 
                                                                       style = "font-size:12px;", class ="btn btn-info btn-sm")),
                                                              plotOutput("gamm_resid_plot", height = 400), plotOutput("gamm_qq_plot", height = 400))),
                                                 tabPanel("Concurvity (GAMM)", 
                                                          div(class = "tab-content-spacing",
                                                          h4("Concurvity (GAMM)", 
                                                          actionButton("info_gamm_concurvity", label = NULL, icon = icon("question-circle"), 
                                                                       style = "font-size:12px;", class ="btn btn-info btn-sm")),
                                                              tableOutput("gamm_concurvity"))))),
                    
                    ### Beta diversity tabs
                    conditionalPanel(condition = "input.diversity_type == 'beta'",
                                     tabsetPanel(tabPanel("PCoA Plot",
                                                          div(class = "tab-content-spacing",
                                                          h4("Principal Coordinates Analysis Plot", 
                                                          actionButton("info_pcoa_plot", label = NULL, icon = icon("question-circle"), 
                                                                        style = "font-size:12px;", class ="btn btn-info btn-sm")),
                                                          plotOutput("pcoa_plot", height = 400))),
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
                                                              tableOutput("dispersion_time")))))
                  )
                )
)

## Server
server <- function(input, output, session) {
  
  ### ALPHA DIVERSITY
  
  # filter for selected treatment groups
  filtered_df <- reactive({
    alpha_diversity$metrics %>%
      left_join(metadata, by = c("sample_name" = "sample_id")) %>%
      filter(gavage %in% input$gavage)
  })
  
  
  ### Linear mixed effects model
  
  # plot alpha diversity metric over time (LMM)
  output$alpha_plot_lmm <- renderPlot({
    ggplot(filtered_df(), aes(x = day, y = .data[[input$alpha_metric]], color = gavage)) +
      geom_jitter(width = 0.2, alpha = 0.5) + theme_minimal() +
      stat_summary(aes(group = gavage), fun = mean, geom = "line", linewidth = 0.8) +
      labs(y = input$alpha_metric, x = "Day")
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
  
  # plot LMM Q-Q
  output$lmm_qq_plot <- renderPlot({
    model <- alpha_diversity$models[[input$alpha_metric]]$lmm$model
    qqnorm(resid(model), main = "LMM Normal Q-Q")
    qqline(resid(model))
  })
  
  # info pop-up for LMM residuals vs fitted and Q-Q plots
  observeEvent(input$info_lmm_resid_plot, {
    showModal(
      modalDialog(
        title = "How to Interpret the Residuals and Q-Q Plots",
        p("The Residuals versus Fitted plot is used to check whether the statistical model is appropriate for the data, whether the remaining errors (residuals) are random 
          or if they show structure the model failed to capture. The points should be randomly scattered around zero (i.e., no clear pattern or trend)."),
        p(strong("Curved pattern:"), "Model failed to capture nonlinear strucutre in the mean (e.g., missing interaction)."),
        p(strong("Curved clusters or bands:"), "Group-specific trends or correlation structure is not adequately modeled."),
        p(strong("Widening/narrowing:"), "Non-constant variance (heteroscedasticity)."),
        p(strong("Extreme outliers:"), "Influential points or possible errors."),
        p("The Q-Q plot (quantile-quantile plot) is used to compare the distribution of the residuals with a theoretical distribution (often normal). 
          Each point represents how one residual compares to what would be expected under normality. LMMs and Gaussian GAMMs assume that residuals are normally distributed. 
          If residuals are not normal, p-values are confidence intervals may be unreliable."),
        p(strong("Points lie on a straight 45° line:"), "Residuals are approximately normally distributed."),
        p(strong("S-shaped curve:"), "Less extreme values (concave up then down) or more extreme values (concave down then up) than normal."),
        p(strong("Points curve away at ends:"), "Indicates skewed residuals."),
        p(strong("Large deviations at the ends:"), "Potential outliers."),
        footer = modalButton("Close")
      )
    )
  })
  
  
  ### Generalized additive mixed model
  
  # plot alpha diversity metric over time (GAMM)
  output$alpha_plot_gamm <- renderPlot({
    ggplot(filtered_df(), aes(x = day, y = .data[[input$alpha_metric]], color = gavage)) +
      geom_jitter(width = 0.2, alpha = 0.5) + theme_minimal() +
      stat_summary(aes(group = gavage), fun = mean, geom = "line", linewidth = 0.8) +
      labs(y = input$alpha_metric, x = "Day")
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
      geom_point(data = filtered_df(), aes(x = day, y = .data[[input$alpha_metric]], color = gavage), alpha = 0.4) +
      geom_line(data = smooth_df, aes(x = day, y = fit, color = gavage), linewidth = 1) +
      geom_ribbon(data = smooth_df, aes(day, ymin = lower, ymax = upper, fill = gavage), alpha = 0.2, color = NA) +
      labs(y = input$alpha_metric, x = "Day")
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
  
  
  # plot GAMM residuals
  output$gamm_resid_plot <- renderPlot({
    model <- alpha_diversity$models[[input$alpha_metric]]$gamm$model$gam
    plot(fitted(model), resid(model),
         xlab = "Fitted values", ylab = "Residuals",
         main = "GAMM Residuals vs Fitted")
    abline(h = 0, lty = 2)
  })
  
  # plot GAMM Q-Q
  output$gamm_qq_plot <- renderPlot({
    res <- alpha_diversity$models[[input$alpha_metric]]$gamm$residuals
    qqnorm(res, main = "GAMM Residuals Q-Q")
    qqline(res)
  })
  
  # info pop-up for GAM residuals vs fitted and Q-Q plots
  observeEvent(input$info_gamm_resid_plot, {
    showModal(
      modalDialog(
        title = "How to Interpret the Residuals and Q-Q Plots",
        p("The Residuals versus Fitted plot is used to check whether the statistical model is appropriate for the data, whether the remaining errors (residuals) are random 
          or if they show structure the model failed to capture. The points should be randomly scattered around zero (i.e., no clear pattern or trend)."),
        p(strong("Curved pattern:"), "Model failed to capture nonlinear strucutre in the mean (e.g., missing interaction)."),
        p(strong("Curved clusters or bands:"), "Group-specific trends or correlation structure is not adequately modeled."),
        p(strong("Widening/narrowing:"), "Non-constant variance (heteroscedasticity)."),
        p(strong("Extreme outliers:"), "Influential points or possible errors."),
        p("The Q-Q plot (quantile-quantile plot) is used to compare the distribution of the residuals with a theoretical distribution (often normal). 
          Each point represents how one residual compares to what would be expected under normality. LMMs and Gaussian GAMMs assume that residuals are normally distributed. 
          If residuals are not normal, p-values are confidence intervals may be unreliable."),
        p(strong("Points lie on a straight 45° line:"), "Residuals are approximately normally distributed."),
        p(strong("S-shaped curve:"), "Less extreme values (concave up then down) or more extreme values (concave down then up) than normal."),
        p(strong("Points curve away at ends:"), "Indicates skewed residuals."),
        p(strong("Large deviations at the ends:"), "Potential outliers."),
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
        p(strong("Parametric coefficients:"), "Represent the linear (fixed) effects in the model. Each group term is the linear difference from the reference group."),
        p(strong("Smooth terms:"), "Indicates whether adding nonlinearity significantly improves the model for a given group."),
        p(strong("edf (estimated degrees of freedom:"), "How wiggly the smooth is. > 1 = some nonlinearity"),
        p(strong("R-sq. (adj):"), "Proportion of the variance explained by the model."),
        p(strong("Scale est."), "In Gaussian models, the scale parameter is an estimate of the residual variance (i.e., the amount of variation in your response variable that your model cannot explain)."),
        footer = modalButton("Close")
      )
    )
  })
  
  
  # GAMM concurvity
  output$gamm_concurvity <- renderTable({
    df <- as.data.frame(alpha_diversity$models[[input$alpha_metric]]$gamm$concurvity)
    
    # format table
    df_formatted <- df %>%
      mutate(para = format(round(para, 2), nsmall = 2),
             `s(day):gavageG_1DMD` = formatC(`s(day):gavageG_1DMD`, format = "e", digits = 4),
             `s(day):gavageG_4DMD` = formatC(`s(day):gavageG_4DMD`, format = "e", digits = 4),
             `s(day):gavageG_4W7C` = formatC(`s(day):gavageG_4W7C`, format = "e", digits = 4),
             `s(day):gavageG_4WMD` = formatC(`s(day):gavageG_4WMD`, format = "e", digits = 4))
    
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
        p(strong("para:"), "The parametric coefficients (see Model Summary (GAMM))."),
        p(strong("s(day):gavage...:"), "The smooth terms/group-specific nonlinear trends over time (see Model Summary (GAMM))."),
        easyClose = TRUE,
        footer = modalButton("Close")
      )
    )
  })
  
  
  ### BETA DIVERSITY
  
  # filter for selected treatment groups
  beta_meta_filtered <- reactive({
    metadata %>% filter(gavage %in% input$gavage)
  })
  
  # reactive data.frame for PCoA plot
  beta_pcoa_df <- reactive({
    metric <- input$beta_metric
    ord <- beta_diversity[[metric]]$ordination$pcoa
    df <- as.data.frame(ord$vectors) # extract coordinates
    df$sample_id <- rownames(df)
    
    # merge with filtered metadata
    df <- left_join(df, beta_meta_filtered(), by = "sample_id")
    
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
      scale_fill_discrete(guide = "none") +
      scale_size_manual(values = c("7" = 1, "14" = 2, "21" = 3, "28" = 4)) +
      xlab(paste0("PC1 (", round(pcoa_data$var_explained[1], 1), "%)")) +
      ylab(paste0("PC2 (", round(pcoa_data$var_explained[2], 1), "%)")) +
      theme_minimal()
  })
  
  # info pop-up for PCoA Plot
  observeEvent(input$info_pcoa_plot, {
    showModal(
      modalDialog(
        title = "Principal coordinate analysis",
        p("PCoA plots represent samples in a low-dimensional Euclidean space with the goal of preserving the pairwise distances or dissimilarities as faithfully as possible. 
        On a PCoA plot, each axis represents the amount of variation in the distance matrix explained by that coordinate, making it interpretable in terms of explained variance."),
        p(strong("Bray-Curtis distance:"), "Measures the compositional dissimilarity between two samples based on the abundances of taxa present in at least one of the samples. 
          It is computed as the weighted sum of absolute differences, where weights are the abundances, thus the metric is dominated by abundant taxa."),
        p(strong("Canberra distance:"), "Measures compositional dissimilarity between two samples based on the abundances of taxa present in at least one of the samples. 
          It sums the ratios of the absolute differences to the sums of the abundances for each taxon, thereby giving relatively equal weight to rare and abundant taxa 
          (i.e., rare taxa contribute proportionally more than they would in Bray-Curtis)."),
        p(strong("Jaccard distance:"), "Measures the dissimilarity between two samples based on presence/absence of taxa, ignoring their abundances. 
          It quantifies how many taxa are shared versus unique to each sample."),
        p(strong("Euclidean distance:"), "The geometric (straight-line) distance between two samples in multidimensional space, where each dimension corresponds to a feature 
          and and the coordinate of a sample along each axis is determined by the abundance of that feature. Euclidean distance assumes absolute values are important, which
          can be misleading with microbiome data."),
        p(strong("Aitchison distance:"), "The Euclidean distance between samples after centered log-ratio (CLR) transformation (which projects the compositional data from the simplex into Euclidean space)."),
        easyClose = TRUE,
        footer = modalButton("Close")
      )
    )
  })

  
  # overall PERMANOVA table
  output$permanova_overall <- renderTable({
    metric <- input$beta_metric
    df <- as.data.frame(beta_diversity[[metric]]$permanova$overall)
    
    # match adonis2 formatting
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
        p("A PERMANOVA tests whether the centroids of groups are significantly different in multivariate space. 
          Permutations to calculate p-values as the distribution of distances is usually non-normal."),
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
        p("A PERMANOVA tests whether the centroids of groups are significantly different in multivariate space. 
          Permutations to calculate p-values as the distribution of distances is usually non-normal."),
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
  
}

# launch Shiny app
shinyApp(ui, server)

