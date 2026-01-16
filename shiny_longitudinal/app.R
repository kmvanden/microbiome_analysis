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
alpha_div_df <- readRDS("alpha_diversity_df.rds")
models <- readRDS("alpha_diversity_models.rds")

# define constants
metrics <- names(models)

## UI
ui <- fluidPage(theme = bs_theme(version = 5, bootswatch = "flatly", primary = "#2C3E50"), # responsive Bootstrap 5 layout
                titlePanel("Alpha Diversity Analysis"),
                sidebarLayout(
                  sidebarPanel(pickerInput(inputId = "metric", label = "Alpha Diversity Metric",
                                           choices = metrics, selected = "Observed"), # alpha diversity metric dropdown
                               pickerInput(inputId = "gavage", label = "Treatment Group",
                                           choices = unique(alpha_div_df$gavage),
                                           selected = unique(alpha_div_df$gavage),
                                           multiple = TRUE), # treatment (gavage) group multi-select dropdown
                               pickerInput(inputId = "model_type", label = "Model Type",
                                           choices = c("LMM", "GAMM"), selected = "LMM")),
                  mainPanel(
                    
                    # LMM tabs
                    conditionalPanel(condition = "input.model_type == 'LMM'",
                                     tabsetPanel(
                                       tabPanel("Alpha Diversity Metric Over Time Plot", plotOutput("alpha_plot_lmm", height = 400)),
                                       tabPanel("ANOVA (LMM)", tableOutput("anova_table")),
                                       tabPanel("Pairwise Comparisons (LMM)", DTOutput("emmeans_table")),
                                       tabPanel("Model Summary (LMM)", verbatimTextOutput("lmm_model_summary")),
                                       tabPanel("Residuals (LMM)", plotOutput("lmm_resid_plot", height = 400), plotOutput("lmm_qq_plot", height = 400)))),
                    
                    # GAMM tabs
                    conditionalPanel(condition = "input.model_type == 'GAMM'",
                                     tabsetPanel(
                                       tabPanel("Alpha Diversity Metric Over Time Plot", plotOutput("alpha_plot_gamm", height = 400)),
                                       tabPanel("Smooth Plot (GAMM)", plotOutput("gamm_smooth_plot", height = 400)),
                                       tabPanel("Model Summary (GAMM)", verbatimTextOutput("gamm_model_summary")),
                                       tabPanel("Residuals (GAMM)", plotOutput("gamm_resid_plot", height = 400), plotOutput("gamm_qq_plot", height = 400)),
                                       tabPanel("Concurvity (GAMM)", DTOutput("gamm_concurvity"))))
                  )))

## Server
server <- function(input, output, session) {
  filtered_df <- reactive({
    alpha_div_df %>%
      filter(gavage %in% input$gavage) # returns filtered data.frame
  })
  
  # plot alpha diversity metric over time (LMM)
  output$alpha_plot_lmm <- renderPlot({
    ggplot(filtered_df(), aes(x = day, y = .data[[input$metric]], color = gavage)) +
      geom_jitter(width = 0.2, alpha = 0.5) + theme_minimal() +
      stat_summary(aes(group = gavage), fun = mean, geom = "line", linewidth = 0.8) +
      labs(title = paste0(input$metric, " over time"), y = input$metric, x = "Day")
  }) 
  
  ### LMM outputs
  # LMM ANOVA table
  output$anova_table <- renderTable({
    aov_tab <- as.data.frame(models[[input$metric]]$lmm$anova)
    aov_tab %>%
      tibble::rownames_to_column("Model Term")
  }, digits = 4)
  
  # LMM estimated marginal means table
  output$emmeans_table <- renderDT({
    datatable(models[[input$metric]]$lmm$emmeans,
              options = list(pageLength = 10))
  })
  
  # LMM summary 
  output$lmm_model_summary <- renderPrint({
    summary(models[[input$metric]]$lmm$model)
  })
  
  # plot LMM residuals versus fitted
  output$lmm_resid_plot <- renderPlot({
    model <- models[[input$metric]]$lmm$model
    plot(fitted(model), resid(model),
         xlab = "Fitted values", ylab = "Residuals",
         main = paste0(input$metric, " - LMM Residuals vs Fitted"))
    abline(h = 0, lty = 2)
  })
  
  # plot LMM Q-Q
  output$lmm_qq_plot <- renderPlot({
    model <- models[[input$metric]]$lmm$model
    qqnorm(resid(model), main = paste(input$metric, "- LMM Normal Q-Q"))
    qqline(resid(model))
  })
  
  ### GAMM outputs
  
  # plot alpha diversity metric over time (GAMM)
  output$alpha_plot_gamm <- renderPlot({
    ggplot(filtered_df(), aes(x = day, y = .data[[input$metric]], color = gavage)) +
      geom_jitter(width = 0.2, alpha = 0.5) + theme_minimal() +
      stat_summary(aes(group = gavage), fun = mean, geom = "line", linewidth = 0.8) +
      labs(title = paste0(input$metric, " over time"), y = input$metric, x = "Day")
  }) 
  
  # GAMM smooth plots
  output$gamm_smooth_plot <- renderPlot({
    smooth_df <- models[[input$metric]]$gamm$smooth_df
    ggplot() + theme_minimal() +
      geom_point(data = filtered_df(), aes(x = day, y = .data[[input$metric]], color = gavage), alpha = 0.4) +
      geom_line(data = smooth_df, aes(x = day, y = fit, color = gavage), linewidth = 1) +
      geom_ribbon(data = smooth_df, aes(day, ymin = lower, ymax = upper, fill = gavage), alpha = 0.2, color = NA) +
      labs(title = paste0("GAMM Smooths for ", input$metric), y = input$metric, x = "Day")
  })
  
  # plot GAMM residuals
  output$gamm_resid_plot <- renderPlot({
    model <- models[[input$metric]]$gamm$model$gam
    plot(fitted(model), resid(model),
         xlab = "Fitted values", ylab = "Residuals",
         main = paste0(input$metric, " - GAMM Residuals vs Fitted"))
    abline(h = 0, lty = 2)
  })
    
  # plot GAMM Q-Q
  output$gamm_qq_plot <- renderPlot({
    res <- models[[input$metric]]$gamm$residuals
    qqnorm(res, main = paste0(input$metric, " - GAMM Residuals Q-Q"))
    qqline(res)
  })

  # GAMM summary
  output$gamm_model_summary <- renderPrint({
    models[[input$metric]]$gamm$summary
  })
  
  # GAMM concurvity
  output$gamm_concurvity <- renderDT({
    conc <- models[[input$metric]]$gamm$concurvity
    conc_df <- as.data.frame(conc)
    
    conc_df[] <- lapply(names(conc_df), function(col_name) {
      if (col_name == "para") {
        round(conc_df[[col_name]], 2) # round "para" column to 2 decimals
      } else {
        formatC(conc_df[[col_name]], format = "e", digits = 2)  # scientific notation for all other columns
      }
    })
    
    colnames(conc_df) <- colnames(conc) # preserve column names
    datatable(conc_df, options = list(pageLength = nrow(conc_df)))
  })
  
}

# launch Shiny app
shinyApp(ui, server)

